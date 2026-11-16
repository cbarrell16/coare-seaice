"""Solve every tile and blend them. Used by :func:`coare_seaice.ice_fluxes`.

Flow: open water from pycoare; each ice tile gets a surface temperature from its method
and fluxes from the ice surface-layer solver; the floe-edge form drag is added to the
momentum flux; fluxes are blended by area fraction.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from . import thermo
from .form_drag import lupkes_form_drag
from .most import DEFAULT_OPTIONS, V2_OPTIONS, bulk_richardson, ice_surface_layer
from .openwater import open_water
from .results import BlendedResult, FluxResult, TileResult
from .roughness import ICE_Z0_DEFAULT, resolve_ice_roughness
from .surface_temperature.methods import AirFreezing, SurfaceTemperature, TileContext
from .surface_temperature.seb import seb_terms, surface_albedo
from .tiles import tile_fractions

NITS = 10  # iterations of the ice surface-layer loop
LIGHT_WIND = 1.0  # m s-1; Grachev et al. (2007) excluded winds of 1 m s-1 and below
# Humidity check over ice, after Andreas et al. (2002, J. Geophys. Res. 107, 8033: "Near-
# surface water vapor over polar sea ice is always near ice saturation"; citation from
# the reference list of Andreas et al. 2010, QJRMS). The paper itself was not available
# to check, and the thresholds are choices of this package, not values from it.
RH_ICE_RANGE = (80.0, 105.0)  # % relative to ice; outside this the sample is flagged
RH_ICE_MIN_SIC = 0.9  # only where the ice concentration is at least this
BLEND_MODES = ("flux", "single_temperature")
_SCHEMES = {"a87": "blended_a87", "ratio": "simple", "prescribed": "prescribed"}


@dataclass
class Forcing:
    """Checked inputs, one value per sample (heights too)."""

    u: np.ndarray  # m s-1 at zu
    t: np.ndarray  # deg C at zt
    rh: np.ndarray  # % over liquid water at zq
    p: np.ndarray  # mb
    sst: np.ndarray  # deg C
    sic: np.ndarray  # 0-1
    zu: np.ndarray
    zt: np.ndarray
    zq: np.ndarray
    lat: np.ndarray
    lw_down: np.ndarray  # W m-2
    sw_down: np.ndarray  # W m-2
    zi: np.ndarray  # m
    rain: np.ndarray | None = None
    cp: np.ndarray | None = None
    sigH: np.ndarray | None = None
    jcool: int = 1
    salinity: float | None = None
    time: np.ndarray | None = None  # s

    @property
    def n(self):
        return self.u.size


def _per_sample(x, n):
    x = np.asarray(x, dtype=float)
    return np.full(n, float(x.reshape(-1)[0])) if x.size == 1 else x


def run_open_water(F, coare="3.5", zrf=10.0, sst=None, jcool=None):
    """pycoare for the open-water tile (reference height ``zrf``).

    Without the cool-skin correction (jcool = 0) the turbulent fluxes do not depend on
    the radiation (tested), but pycoare returns NaN if it is NaN. Missing radiation is
    then replaced by zero for this call only, so that data without radiation (aircraft
    runs) can be used; nothing that is returned depends on the replacement.
    """
    jcool = F.jcool if jcool is None else jcool
    rs, rl = F.sw_down, F.lw_down
    if jcool == 0:
        rs, rl = np.where(np.isfinite(rs), rs, 0.0), np.where(np.isfinite(rl), rl, 0.0)
    return open_water(
        F.u, F.t, F.rh, F.sst if sst is None else sst, p=F.p, rs=rs, rl=rl,
        zu=F.zu, zt=F.zt, zq=F.zq, lat=F.lat, zi=F.zi, rain=F.rain, cp=F.cp, sigH=F.sigH,
        jcool=jcool, zrf=zrf, version=coare,
        salinity=F.salinity,
    )  # fmt: skip


def _water_tile(tile, fraction, ow, F, sst):
    t_skin = sst - ow.dter * ow.jcool
    q_skin = ow.q_sea - ow.dqer * ow.jcool
    rib = bulk_richardson(
        F.t + thermo.LAPSE_DRY * F.zt, t_skin, ow.q_air, q_skin, F.u, F.zu,
        F.t + thermo.TDK, thermo.gravity(F.lat),
    )  # fmt: skip
    used = [F.u, F.t, F.rh, F.p, F.lat, sst]
    if ow.jcool:  # the cool-skin correction needs the radiation
        used += [F.lw_down, F.sw_down]
    missing = ~np.all([np.isfinite(x) for x in used], axis=0)
    return TileResult(
        name=tile.name, kind="water", fraction=fraction, tau=ow.tau, hsb=ow.hsb, hlb=ow.hlb,
        usr=ow.usr, t_surface=t_skin, q_surface=q_skin, z0=ow.zo, z0t=ow.zot, z0q=ow.zoq,
        cdn10=ow.cdn_rf, chn10=ow.chn_rf, cen10=ow.cen_rf, L=ow.obukhov_length, zeta=ow.zet,
        rib=rib, wind_eff=ow.ut,
        qc={"rib_above_critical": rib > DEFAULT_OPTIONS.rib_critical, "missing_input": missing},
        extra={"open_water": ow, "u10n": ow.u_n_rf},
    )  # fmt: skip


def _prescribed_water_tile(tile, fraction, F, ow, sst, options):
    """A water tile with given roughness lengths (no COARE, no cool skin).

    It uses COARE's stability functions and no light-wind terms, whatever the ice tiles
    use; the stability limits and flag thresholds are those of ``options``. Surface
    humidity is 0.98 x saturation over water, also for COARE 3.6 with a salinity.
    """
    n = F.n
    water_options = replace(options, stability="coare", gustiness=False, windless=False)
    fl = ice_surface_layer(
        u=F.u, t=F.t, Q_air=ow.q_air, t_ice=sst, z0_ice=_per_sample(tile.z0, n), P=F.p,
        grav=thermo.gravity(F.lat), zu=F.zu, zt=F.zt, zq=F.zq, rhoa=ow.rhoa, cpa=thermo.CPA,
        von=thermo.VON, fdg=thermo.FDG, ta=F.t + thermo.TDK, ice_scheme="prescribed",
        scalar_ratio=1.0, nits=NITS, options=water_options, z0T=tile.z0t, z0q=tile.z0q,
        surface="water",
    )  # fmt: skip
    u10n = fl["usr"] / thermo.VON * np.log(10.0 / fl["z0"])
    missing = ~np.all([np.isfinite(x) for x in (F.u, F.t, F.rh, F.p, F.lat, sst)], axis=0)
    return TileResult(
        name=tile.name, kind="water", fraction=fraction, tau=fl["tau"], hsb=fl["hsb"],
        hlb=fl["hlb"], usr=fl["usr"], t_surface=sst, q_surface=fl["Qs"], z0=fl["z0"],
        z0t=fl["z0T"], z0q=fl["z0q"], cdn10=fl["Cdn"] / 1e3, chn10=fl["Chn"] / 1e3,
        cen10=fl["Cen"] / 1e3, rstar=fl["R_star"], L=fl["L"], zeta=fl["zeta"], rib=fl["Rib"],
        wind_eff=fl["S"],
        qc={"rib_above_critical": fl["qc_rib"], "zeta_limited": fl["zeta_limited"],
            "not_converged": fl["not_converged"], "missing_input": missing},
        extra={"u10n": u10n},
    )  # fmt: skip


def _ice_tile(tile, fraction, F, ow, q_air, default_method, default_z0, options, v2_compat,
              t_fixed=None):  # fmt: skip
    """Surface temperature, fluxes and energy-balance terms for one ice tile."""
    n = F.n
    method = tile.t_surface if tile.t_surface is not None else default_method
    z0 = resolve_ice_roughness(default_z0 if tile.z0 is None else tile.z0)
    z0 = z0 if callable(z0) else _per_sample(z0, n)

    def setting(name):  # the method's own value wins over the tile's
        own = getattr(method, name, None)
        return getattr(tile, name) if own is None else own

    h_ice, h_snow = _per_sample(setting("h_ice"), n), _per_sample(setting("h_snow"), n)
    albedo = setting("albedo")
    albedo = surface_albedo(h_snow) if albedo is None else _per_sample(albedo, n)
    k_snow, t_base = setting("k_snow"), setting("t_base")

    args = dict(
        u=F.u, t=F.t, Q_air=q_air, z0_ice=z0, P=F.p, grav=thermo.gravity(F.lat),
        zu=F.zu, zt=F.zt, zq=F.zq, rhoa=ow.rhoa, cpa=thermo.CPA, von=thermo.VON,
        fdg=thermo.FDG, ta=F.t + thermo.TDK,
        ice_scheme=_SCHEMES[tile.scalar_scheme], scalar_ratio=tile.scalar_ratio, nits=NITS,
        options=options, z0T=tile.z0t, z0q=tile.z0q,
    )  # fmt: skip

    def flux(Ts):
        return ice_surface_layer(t_ice=Ts, **args)

    ctx = TileContext(
        t_air=F.t, zt=F.zt, p=F.p, q_air=q_air, lw_down=F.lw_down, sw_down=F.sw_down,
        sic=F.sic, flux=flux, h_ice=h_ice, h_snow=h_snow, k_snow=k_snow, albedo=albedo,
        t_base=t_base, time=F.time, v2_compat=v2_compat,
    )  # fmt: skip
    st = method.solve(ctx) if t_fixed is None else SurfaceTemperature(t_fixed)
    Ts = st.t_surface
    fl = flux(Ts)
    lw_net, sw_net, f_cond = seb_terms(Ts, F.lw_down, F.sw_down, albedo, ctx.conductance, t_base)
    balanced = getattr(method, "energy_balance", False) and t_fixed is None
    if balanced and st.f_cond is not None:
        f_cond = st.f_cond
    if balanced and st.seb_residual is not None:
        residual = st.seb_residual
    else:  # the imbalance implied by a surface temperature that was not solved for
        residual = f_cond + lw_net + sw_net - fl["hsb"] - fl["hlb"]

    used = [F.u, F.t, F.rh, F.p, F.lat, F.sic, Ts]
    if balanced:
        used += [F.lw_down, F.sw_down, h_ice, h_snow, albedo]
    filled = np.zeros(n, dtype=bool) if st.filled is None else st.filled
    return TileResult(
        name=tile.name, kind="ice", fraction=fraction, tau=fl["tau"], hsb=fl["hsb"],
        hlb=fl["hlb"], usr=fl["usr"], t_surface=Ts, q_surface=fl["Qs"], z0=fl["z0"],
        z0t=fl["z0T"],
        z0q=fl["z0q"], cdn10=fl["Cdn"] / 1e3, chn10=fl["Chn"] / 1e3, cen10=fl["Cen"] / 1e3,
        rstar=fl["R_star"], L=fl["L"], zeta=fl["zeta"], rib=fl["Rib"], wind_eff=fl["S"],
        lw_net=lw_net, sw_net=sw_net, f_cond=f_cond, seb_residual=residual,
        qc={
            "rib_above_critical": fl["qc_rib"], "zeta_limited": fl["zeta_limited"],
            "rstar_out_of_range": fl["qc_rstar"], "not_converged": fl["not_converged"],
            "t_surface_filled": filled,
            "missing_input": ~np.all([np.isfinite(x) for x in used], axis=0),
        },
        extra={"iter_dH": fl["iter_dH"], "iter_dLE": fl["iter_dLE"]},
    )  # fmt: skip


def compute(
    F,
    tiles,
    *,
    blend="flux",
    coare="3.5",
    t_ice=None,
    z0_ice=None,
    options=None,
    rh_reference="water",
    t_surface=None,
    open_water_result=None,
    v2_compat=False,
):
    """Solve the tiles for forcing ``F`` and blend them. Returns a FluxResult.

    blend="flux": each tile has its own surface temperature and stability, and the
    fluxes are averaged by area (a mosaic).

    blend="single_temperature": every tile is solved with one aggregate surface
    temperature. The tiles still differ in roughness, in surface humidity (0.98 x
    saturation over water for the water tile, saturation over ice for ice tiles, at that
    temperature; at -16.5 C the water value is about 15 % higher), in latent heat, in
    stability correction, and over water in COARE's gustiness and its stress convention.
    For sensible heat in neutral conditions the area-mean flux is the area mean of the
    tiles' neutral transfer coefficients times one temperature difference, the blending
    of Elvidge et al. (2021, Eq. 12; there at 10 m, here at the measurement height).
    Their offline fluxes also use one surface humidity and one latent heat, which this
    mode does not. The aggregate temperature is ``t_surface`` (deg C, for example from a
    radiometer); if None it is the area mean of the tile temperatures from a first pass
    in flux mode (as Andreas et al. 2010, Eq. 6.4). The open-water tile is run at that
    temperature in place of the SST and without the cool-skin correction, even below the
    freezing point of sea water. Ice tiles cap the surface at 0 C.
    """
    if blend not in BLEND_MODES:
        raise ValueError(f"blend must be one of {BLEND_MODES}")
    options = V2_OPTIONS if v2_compat else (options or DEFAULT_OPTIONS)
    default_method = AirFreezing() if t_ice is None else t_ice
    default_z0 = ICE_Z0_DEFAULT if z0_ice is None else z0_ice
    fractions = tile_fractions(tiles, F.sic)
    if sum(t.kind == "water" for t in tiles) != 1:
        raise ValueError("a tile set needs exactly one water tile")

    ow = open_water_result if open_water_result is not None else run_open_water(F, coare)
    if rh_reference == "water":
        q_air = ow.q_air
    elif rh_reference == "ice":  # v1 behaviour: RH taken as relative to ice
        em = F.rh / 100.0 * thermo.esat_ice(np.minimum(F.t, 0.0), F.p)
        q_air = 0.622 * em / (F.p - 0.378 * em)
    else:
        raise ValueError("rh_reference must be 'water' or 'ice'")

    def solve_all(ow, sst, t_fixed=None):
        out = {}
        for tile in tiles:
            f = fractions[tile.name]
            if tile.kind == "water" and tile.z0 is not None:
                out[tile.name] = _prescribed_water_tile(tile, f, F, ow, sst, options)
            elif tile.kind == "water":
                out[tile.name] = _water_tile(tile, f, ow, F, sst)
            else:
                out[tile.name] = _ice_tile(tile, f, F, ow, q_air, default_method, default_z0,
                                           options, v2_compat, t_fixed)  # fmt: skip
        return out

    if blend == "flux":
        if t_surface is not None:
            raise ValueError("t_surface is only used with blend='single_temperature'")
        results = solve_all(ow, F.sst)
    else:
        if t_surface is None:
            first = solve_all(ow, F.sst).values()
            t_agg = sum(r.fraction * r.t_surface for r in first)
        else:
            t_agg = _per_sample(t_surface, F.n)
        ow = run_open_water(F, coare, sst=t_agg, jcool=0)
        results = solve_all(ow, t_agg, t_fixed=t_agg)
    water = next(r for r in results.values() if r.kind == "water")

    # Floe-edge form drag (momentum only): neutral 10 m coefficient x U10N^2 over water.
    a_form = sum(fractions[t.name] for t in tiles if t.form_drag)
    cdn10_form = lupkes_form_drag(a_form, water.z0)
    tau_form = ow.rhoa * cdn10_form * water.extra["u10n"] ** 2

    def area_mean(name):
        return sum(r.fraction * getattr(r, name) for r in results.values())

    blended = BlendedResult(
        tau=area_mean("tau") + tau_form, hsb=area_mean("hsb"), hlb=area_mean("hlb"),
        tau_form=tau_form, t_surface=area_mean("t_surface") if blend == "flux" else t_agg,
        cdn10=area_mean("cdn10") + cdn10_form, chn10=area_mean("chn10"),
        cen10=area_mean("cen10"), cdn10_form=cdn10_form,
    )  # fmt: skip
    meta = dict(blend=blend, coare=coare, tiles=",".join(t.name for t in tiles),
                rh_reference=rh_reference, options=options)  # fmt: skip
    qc = _quality_flags(F, results, rh_reference)
    # whatever the cause, a NaN blended flux is flagged as missing input
    qc["missing_input"] |= ~(np.isfinite(blended.tau) & np.isfinite(blended.hsb)
                             & np.isfinite(blended.hlb))  # fmt: skip
    return FluxResult(blended=blended, tiles=results, qc=qc, meta=meta)


def relative_humidity_over_ice(t, p, rh):
    """Relative humidity with respect to ice (%) from that with respect to water (%).

    Below 0 C saturation over ice is lower than over water, so the value is higher:
    100 % over ice is about 91 % over water at -10 C and 82 % at -20 C. Above 0 C there
    is no saturation over ice at the air temperature, so the result is NaN (and the
    humidity flag is not set).
    """
    t = np.asarray(t, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = rh * thermo.esat_water(t, p) / thermo.esat_ice(np.minimum(t, 0.0), p)
    return np.where(t > 0.0, np.nan, out)


def _quality_flags(F, results, rh_reference="water"):
    """Flags for the blended result: a tile counts where its area fraction is not zero."""
    n = F.n
    rh_ice = F.rh if rh_reference == "ice" else relative_humidity_over_ice(F.t, F.p, F.rh)

    def any_tile(flag=None, test=None, kind=None):
        out = np.zeros(n, dtype=bool)
        for r in results.values():
            if kind and r.kind != kind:
                continue
            hit = r.qc.get(flag, np.zeros(n, dtype=bool)) if test is None else test(r)
            out |= hit & ~(r.fraction == 0)
        return out

    return {
        "missing_input": any_tile("missing_input") | ~np.isfinite(F.sic),
        "rib_above_critical": any_tile("rib_above_critical"),
        "zeta_limited": any_tile("zeta_limited"),
        "height_above_obukhov_length": any_tile(test=lambda r: r.zeta > 1.0),
        "light_wind": F.u < LIGHT_WIND,
        "rstar_out_of_range": any_tile("rstar_out_of_range"),
        "not_converged": any_tile("not_converged"),
        "t_surface_filled": any_tile("t_surface_filled"),
        "melting_surface": any_tile(test=lambda r: r.t_surface >= 0.0, kind="ice"),
        "humidity_far_from_ice_saturation": (F.sic >= RH_ICE_MIN_SIC)
        & ((rh_ice < RH_ICE_RANGE[0]) | (rh_ice > RH_ICE_RANGE[1])),
    }
