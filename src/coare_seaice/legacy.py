"""The v2 interface: ``coare35vn_ice(...)`` returning the 56-column array.

Same call signature, defaults and output columns as ``coare35vn_ice.py`` version 2.0.0
(Nov 2026). Open-water columns 0-37 come from ``pycoare`` (through
:mod:`coare_seaice.openwater`); ice columns 38-55 from the ported ice modules.

Open-water columns differ slightly from v2, which ran its own copy of the NOAA COARE 3.5
loop: pycoare uses WGS84 gravity, exact 1/3 and 2/3 exponents and constants in a few
places, where the NOAA code has 0.333 / 0.3333 / 0.6667. Ice-only columns do not depend
on these; the blended columns and form drag inherit the small differences through the
open-water fluxes, roughness ``zo`` and 10 m neutral wind.

Columns (heat fluxes positive upward):
    0-37   open water, as the original coare35vn (see ``OUTPUT_COLUMNS``)
    38-46  tau_ice, hsb_ice, hlb_ice, Cdn_ice, Chn_ice, Cen_ice (x1e3, 10 m neutral),
           R_star, z0T_ice, z0q_ice
    47-49  tau_blend = (1-A) tau_w + A tau_i + tau_form; hsb_blend; hlb_blend
    50-55  T_ice (C), F_cond (+up), LWnet_ice (+down), SWnet_ice (+down),
           SEB_res (residual or melt flux), tau_form
"""

from __future__ import annotations

import warnings

import numpy as np

from . import thermo
from .form_drag import lupkes_form_drag
from .most import DEFAULT_OPTIONS, V2_OPTIONS, ice_surface_layer
from .openwater import open_water
from .roughness import ICE_Z0_DEFAULT, SCALAR_RATIO
from .surface_temperature import simple
from .surface_temperature.column import IceColumn
from .surface_temperature.properties import K_ICE_FRESH, K_SNOW
from .surface_temperature.seb import (
    effective_conductance,
    seb_terms,
    solve_seb_steady,
    surface_albedo,
)

__all__ = ["OUTPUT_COLUMNS", "T_ICE_METHODS", "coare35vn_ice"]

ICE_T_DEFAULT = thermo.T_FREEZE_SW  # deg C, ice base temperature
T_ICE_METHODS = ("air_freezing", "air_capped", "empirical", "prescribed", "seb_steady", "diffusion")
ZRF_U, ZRF_TQ = 10.0, 2.0  # v2 reference heights: wind 10 m, temperature and humidity 2 m
NITS = 10  # iterations of the ice surface-layer loop (as the v2 open-water loop)

OUTPUT_COLUMNS = [
    "usr", "tau", "hsb", "hlb", "hbb", "hsbb", "hlwebb", "tsr", "qsr", "zot", "zoq",
    "Cd", "Ch", "Ce", "L", "zet", "dter", "dqer", "tkt", "Urf", "Trf", "Qrf", "RHrf",
    "UrfN", "Rnl", "Le", "rhoa", "UN", "U10", "U10N", "Cdn_10", "Chn_10", "Cen_10",
    "RF", "Evap", "Qs", "Q10", "RH10",
    "tau_ice", "hsb_ice", "hlb_ice", "Cdn_ice", "Chn_ice", "Cen_ice", "R_star",
    "z0T_ice", "z0q_ice", "tau_blend", "hsb_blend", "hlb_blend",
    "T_ice", "F_cond", "LWnet_ice", "SWnet_ice", "SEB_res", "tau_form",
]  # fmt: skip


def _as_float(x):
    return np.array(x, dtype=float)  # always a copy


def _per_sample(x, n, name):
    """Return x as a length-n float array; x must have size 1 or n."""
    x = _as_float(x)
    if x.size == 1:
        return x.reshape(()) * np.ones(n)
    if x.size != n:
        raise ValueError(f"coare35vn_ice: {name} array has wrong length")
    return x


def coare35vn_ice(
    u,
    t,
    rh,
    ts,
    P=1015,
    Rs=150,
    Rl=370,
    zu=18,
    zt=18,
    zq=18,
    lat=45,
    zi=600,
    rain=None,
    cp=None,
    sigH=None,
    jcool=1,
    ice_conc=None,
    t_ice=None,
    z0_ice=ICE_Z0_DEFAULT,
    ice_scheme="blended_a87",
    scalar_ratio=SCALAR_RATIO,
    t_ice_method="seb_steady",
    t_ice_fallback="seb_steady",
    h_ice=1.0,
    h_snow=0.2,
    k_snow=K_SNOW,
    albedo=None,
    t_base=ICE_T_DEFAULT,
    t_ice_offset=2.0,
    time_s=None,
    diffusion_kwargs=None,
    rh_reference="water",
    verbose=True,
    surface_layer_options=None,
    return_diagnostics=False,
    _v2_compat=False,
):
    """COARE 3.5 with the sea-ice extension, v2 interface. Returns an (N, 56) array.

    Open-water inputs: u (m/s at zu), t (C at zt), rh (% over liquid water at zq),
    ts (C), P (mb), Rs, Rl (downwelling SW, LW, W m-2; also used in the ice SEB), zu, zt,
    zq (m; scalars of any real number type, or one value per sample), lat, zi (m),
    rain (mm/h), cp, sigH, jcool.

    Ice inputs:
        ice_conc        sea-ice concentration 0-1 (None: open water only, 38 columns)
        t_ice_method    'air_freezing' | 'air_capped' | 'empirical' | 'prescribed' |
                        'seb_steady' (default) | 'diffusion'
        t_ice           prescribed T_ice (C); giving it selects 'prescribed'; NaNs are
                        filled by t_ice_fallback
        h_ice, h_snow   ice and snow thickness (m), scalar or array (SEB methods). For
                        'diffusion' the column is one slab of fixed thickness: arrays
                        (and an albedo array) are replaced by their medians everywhere in
                        that method, with a warning
        k_snow          snow conductivity (W m-1 K-1)
        albedo          surface albedo (default 0.85 if h_snow > 2 cm else 0.65)
        t_base          ice-base temperature (C)
        t_ice_offset    dT for 'empirical' (K)
        time_s          time (s), required for 'diffusion'
        diffusion_kwargs  passed to IceColumn (and 'sic_reset', 'max_gap_s' to its run())
        rh_reference    'water' (default) or 'ice' (v1: RH taken as relative to ice)
        z0_ice, ice_scheme ('blended_a87' | 'simple' | 'open_ocean'), scalar_ratio

    Additions to the v2 interface (keyword only in practice):
        surface_layer_options  a :class:`coare_seaice.most.SurfaceLayerOptions`
                        (stability limit and flag threshold); default as v2
        return_diagnostics  if True, return ``(array, diagnostics)``; ``diagnostics`` is a
                        dict of per-sample arrays for the ice fraction: ``zeta``, ``Rib``,
                        ``iter_dH``, ``iter_dLE`` and the flags ``zeta_limited``,
                        ``qc_rib``, ``qc_rstar`` and ``not_converged`` (see
                        :func:`coare_seaice.most.ice_surface_layer`); ``t_ice_filled``
                        (T_ice does not come from the requested method: a prescribed
                        value was missing and ``t_ice_fallback`` was used, or the
                        'diffusion' column was not running and the steady value was
                        used); ``missing_input`` (an input needed for the ice fluxes is
                        NaN, so they are NaN). Empty without ice. The threshold flags
                        (``qc_rib``, ``zeta_limited``, ``qc_rstar``) are False where
                        the inputs are missing: read them together with
                        ``missing_input`` (the fluxes are NaN there anyway).

    Missing data: NaN in gives NaN out for every output that depends on the input, and
    nothing is filled silently. ``t_ice_fallback`` is the one explicit fill; pass
    ``t_ice_fallback=None`` to leave missing prescribed temperatures as NaN.

    ``_v2_compat=True`` switches off every audit fix of release 0.3.0 that changes
    results. It exists only so that the tests can still compare the port with the stored
    v2 outputs; do not use it for science.
    """
    # ---- open-water inputs -------------------------------------------------------------
    u, t, rh, ts = _as_float(u), _as_float(t), _as_float(rh), _as_float(ts)
    N = u.size
    if not (t.size == N and rh.size == N and ts.size == N):
        raise ValueError("coare35vn_ice: u, t, rh, ts arrays of different length")
    P, Rl, Rs, zi, lat = (
        _per_sample(x, N, name)
        for x, name in ((P, "P"), (Rl, "Rl"), (Rs, "Rs"), (zi, "zi"), (lat, "lat"))
    )
    if rain is not None:
        rain = _as_float(rain)
        if rain.size != N:
            raise ValueError("coare35vn_ice: rain array has wrong length")
    if cp is not None:
        cp = _per_sample(cp, N, "cp")
    if sigH is not None:
        sigH = _per_sample(sigH, N, "sigH")
    if cp is not None and sigH is not None:
        print("Using seastate dependent parameterization")
    elif cp is not None:
        print("Using waveage dependent parameterization")
    jcool = 0 if jcool == 0 else 1
    # Heights: any real number type, scalar or one value per sample.
    zu_a, zt_a, zq_a = (_per_sample(z, N, name) for z, name in ((zu, "zu"), (zt, "zt"), (zq, "zq")))
    for z, name in ((zu_a, "zu"), (zt_a, "zt"), (zq_a, "zq")):
        if np.any(z <= 0):  # NaN heights pass and give NaN fluxes
            raise ValueError(f"coare35vn_ice: {name} must be positive")

    # ---- ice inputs ----------------------------------------------------------------------
    ice_active = ice_conc is not None and ice_scheme != "open_ocean"
    if ice_active:
        A_ice = _per_sample(ice_conc, N, "ice_conc")
        A_ice = np.where(np.isfinite(A_ice), np.clip(A_ice, 0.0, 1.0), np.nan)
        if t_ice_method not in T_ICE_METHODS:
            raise ValueError(f"coare35vn_ice: t_ice_method must be one of {T_ICE_METHODS}")
        if t_ice_method == "prescribed" and t_ice is None:
            raise ValueError("coare35vn_ice: t_ice_method='prescribed' needs t_ice")
        if t_ice_method == "diffusion" and time_s is None:
            raise ValueError("coare35vn_ice: t_ice_method='diffusion' needs time_s")
        if t_ice is not None and t_ice_method != "prescribed":
            t_ice_method = "prescribed"  # v1 behaviour: an explicit t_ice means 'prescribed'
        if rh_reference not in ("water", "ice"):
            raise ValueError("coare35vn_ice: rh_reference must be 'water' or 'ice'")
        if np.asarray(z0_ice).size != 1:
            raise ValueError("coare35vn_ice: z0_ice should be a scalar constant")
        z0_ice_val = float(z0_ice)
        h_ice_v = _per_sample(h_ice, N, "h_ice")
        h_snow_v = _per_sample(h_snow, N, "h_snow")
        alb_v = surface_albedo(h_snow_v) if albedo is None else _per_sample(albedo, N, "albedo")
        if verbose:
            print(
                f"Ice scheme: {ice_scheme!r}, T_ice: {t_ice_method!r}, "
                f"z0_ice = {z0_ice_val:.2e} m, h_ice = {np.nanmedian(h_ice_v):.2f} m, "
                f"h_snow = {np.nanmedian(h_snow_v):.2f} m, RH ref: {rh_reference}"
            )

    # ---- open water (pycoare) ------------------------------------------------------------
    ow_kw = dict(
        p=P, rs=Rs, rl=Rl, zu=zu_a, zt=zt_a, zq=zq_a, lat=lat, zi=zi,
        rain=rain, cp=cp, sigH=sigH, jcool=jcool, version="3.5",
    )  # fmt: skip
    ow_u = open_water(u, t, rh, ts, zrf=ZRF_U, **ow_kw)
    ow_tq = open_water(u, t, rh, ts, zrf=ZRF_TQ, **ow_kw)
    out_ow = _open_water_columns(ow_u, ow_tq, rain is not None)
    if not ice_active:
        out = np.column_stack(out_ow)
        return (out, {}) if return_diagnostics else out

    out_ice, diag = _ice_columns(
        u=u, t=t, rh=rh, P=P, Rs=Rs, Rl=Rl, grav=thermo.gravity(lat),
        zu=zu_a, zt=zt_a, zq=zq_a, q_air=ow_u.q_air, rhoa=ow_u.rhoa,
        zo_w=ow_u.zo, u10n_w=ow_u.u_n_rf, tau_w=ow_u.tau, hsb_w=ow_u.hsb, hlb_w=ow_u.hlb,
        A_ice=A_ice, z0_ice=z0_ice_val, ice_scheme=ice_scheme, scalar_ratio=scalar_ratio,
        t_ice_method=t_ice_method, t_ice=t_ice, t_ice_fallback=t_ice_fallback,
        h_ice=h_ice_v, h_snow=h_snow_v, k_snow=k_snow, albedo=alb_v, t_base=t_base,
        t_ice_offset=t_ice_offset, time_s=time_s, diffusion_kwargs=diffusion_kwargs,
        rh_reference=rh_reference, verbose=verbose,
        options=V2_OPTIONS if _v2_compat else (surface_layer_options or DEFAULT_OPTIONS),
        v2_compat=_v2_compat,
    )  # fmt: skip
    out = np.column_stack(out_ow + out_ice)
    return (out, diag) if return_diagnostics else out


def _open_water_columns(ow_u, ow_tq, rain_given):
    """The 38 v2 open-water columns from two pycoare runs (zrf = 10 m and 2 m)."""
    jc = ow_u.jcool
    dqer = ow_u.dqer * jc
    # v2 returns zeros for the rain heat flux when no rain is given; pycoare returns NaN.
    rf = ow_u.rf if rain_given else np.zeros(ow_u.usr.size)
    return [
        ow_u.usr, ow_u.tau, ow_u.hsb, ow_u.hlb, ow_u.hbb, ow_u.hsbb, ow_u.hlwebb,
        ow_u.tsr, ow_u.qsr * 1000, ow_u.zot, ow_u.zoq,
        ow_u.cd, ow_u.ch, ow_u.ce, ow_u.obukhov_length, ow_u.zet,
        ow_u.dter, dqer, ow_u.tkt,
        ow_u.u_rf, ow_tq.t_rf, ow_tq.q_rf * 1000, ow_tq.rh_rf, ow_u.u_n_rf,
        ow_u.rnl, ow_u.le, ow_u.rhoa, ow_u.u_n, ow_u.u_rf, ow_u.u_n_rf,
        1000 * ow_u.cdn_rf, 1000 * ow_u.chn_rf, 1000 * ow_u.cen_rf,
        rf, ow_u.evap, (ow_u.q_sea - dqer) * 1000, ow_u.q_rf * 1000, ow_u.rh_rf,
    ]  # fmt: skip


def _ice_columns(
    *,
    u,
    t,
    rh,
    P,
    Rs,
    Rl,
    grav,
    zu,
    zt,
    zq,
    q_air,
    rhoa,
    zo_w,
    u10n_w,
    tau_w,
    hsb_w,
    hlb_w,
    A_ice,
    z0_ice,
    ice_scheme,
    scalar_ratio,
    t_ice_method,
    t_ice,
    t_ice_fallback,
    h_ice,
    h_snow,
    k_snow,
    albedo,
    t_base,
    t_ice_offset,
    time_s,
    diffusion_kwargs,
    rh_reference,
    verbose,
    options=DEFAULT_OPTIONS,
    v2_compat=False,
):
    """The 18 v2 ice columns (38-55), given the open-water quantities they depend on.

    zo_w and u10n_w are the open-water roughness and 10 m neutral wind (form drag);
    tau_w, hsb_w, hlb_w the open-water fluxes (blending). All arrays have length N.
    Returns (columns, diagnostics).
    """
    N = u.size
    # Air humidity over the ice fraction (kg/kg)
    if rh_reference == "water":
        Q_ice_air = q_air
    else:  # v1 behaviour: RH taken as relative to ice
        es_i = thermo.esat_ice(np.minimum(t, 0.0), P)
        em_i = rh / 100.0 * es_i
        Q_ice_air = 0.622 * em_i / (P - 0.378 * em_i)

    args = dict(
        u=u, t=t, Q_air=Q_ice_air, z0_ice=z0_ice, P=P, grav=grav, zu=zu, zt=zt, zq=zq,
        rhoa=rhoa, cpa=thermo.CPA, von=thermo.VON, fdg=thermo.FDG, ta=t + thermo.TDK,
        ice_scheme=ice_scheme, scalar_ratio=scalar_ratio, nits=NITS, options=options,
    )  # fmt: skip
    G = effective_conductance(h_ice, h_snow, K_ICE_FRESH, k_snow)

    def _simple(method):
        if method == "air_freezing":
            return simple.air_freezing(t, ICE_T_DEFAULT)
        if method == "air_capped":
            return simple.air_capped(t, 0.0)
        if method == "empirical":
            return simple.empirical(t, offset=t_ice_offset)
        if method == "seb_steady":
            return solve_seb_steady(args, Rl, Rs, albedo, G, t_base, t)[0]
        raise ValueError(f"unsupported fallback {method!r}")

    seb_res = np.full(N, np.nan)
    F_cond_column = None
    t_ice_filled = np.zeros(N, dtype=bool)
    if t_ice_method == "prescribed":
        T_ice = _per_sample(t_ice, N, "t_ice")
        miss = ~np.isfinite(T_ice)
        if miss.any() and t_ice_fallback is not None:
            T_ice[miss] = _simple(t_ice_fallback)[miss]
            t_ice_filled = miss & np.isfinite(T_ice)
            if verbose:
                print(
                    f"  prescribed T_ice missing in {miss.sum()} of {N} steps -> {t_ice_fallback!r}"
                )
    elif t_ice_method == "seb_steady":
        T_ice, seb_res = solve_seb_steady(args, Rl, Rs, albedo, G, t_base, t)
    elif t_ice_method == "diffusion":
        # Samples with a missing thickness or albedo get no T_ice (set after the run).
        no_slab = ~(np.isfinite(h_ice) & np.isfinite(h_snow) & np.isfinite(albedo))
        if not v2_compat:
            varying = [
                name
                for name, x in (("h_ice", h_ice), ("h_snow", h_snow), ("albedo", albedo))
                if np.nanmax(x) - np.nanmin(x) > 1e-9
            ]
            if varying:
                warnings.warn(
                    f"t_ice_method='diffusion' models one slab: {', '.join(varying)} "
                    "vary along the series and are replaced by their medians",
                    stacklevel=3,
                )
            h_ice = np.full(N, np.nanmedian(h_ice))
            h_snow = np.full(N, np.nanmedian(h_snow))
            albedo = np.full(N, np.nanmedian(albedo))
            G = effective_conductance(h_ice, h_snow, K_ICE_FRESH, k_snow)
        T_ss, _ = solve_seb_steady(args, Rl, Rs, albedo, G, t_base, t)
        fl_ss = ice_surface_layer(t_ice=T_ss, **args)  # conductances at steady state
        dk = dict(diffusion_kwargs or {})
        run_kw = {k: dk.pop(k) for k in ("sic_reset", "max_gap_s") if k in dk}
        # One column is one slab of fixed thickness. Thickness arrays are reduced to their
        # medians (see the docstring); v2 did that for the column but not for T_ss above.
        col = IceColumn(
            h_ice=float(np.nanmedian(h_ice)),
            h_snow=float(np.nanmedian(h_snow)),
            t_base=t_base,
            k_snow=k_snow,
            v2_compat=v2_compat,
            **dk,
        )
        theta_a = t + thermo.LAPSE_DRY * zt
        res = col.run(
            np.asarray(time_s, dtype=float), theta_a, Q_ice_air, P, Rl, Rs,
            fl_ss["h_s"], fl_ss["h_l"], T_ss, sic=A_ice,
            albedo=float(np.nanmedian(albedo)), **run_kw,
        )  # fmt: skip
        # The column is not run out of ice (SIC below sic_reset) or across missing
        # forcing. There the steady-state T_ice is used and flagged (t_ice_filled); it is
        # NaN itself wherever an input it needs is missing.
        T_ice = np.where(np.isfinite(res["Ts"]), res["Ts"], T_ss)
        t_ice_filled = ~np.isfinite(res["Ts"]) & np.isfinite(T_ss)
        if not v2_compat:  # the medians must not stand in for a missing value
            T_ice = np.where(no_slab, np.nan, T_ice)
            t_ice_filled &= ~no_slab
        seb_res = res["melt"]
        if v2_compat:
            # v2 back-calculated F_cond through an effective conductance, which is wrong
            # (and can be enormous) when the surface is warmer than the base.
            G = np.where(
                np.isfinite(res["F_cond"]),
                res["F_cond"] / np.maximum(t_base - np.minimum(T_ice, 0.0), 1e-6),
                G,
            )
        else:
            F_cond_column = res["F_cond"]  # the column model's own conductive flux
    else:
        T_ice = _simple(t_ice_method)

    fl = ice_surface_layer(t_ice=T_ice, **args)
    lw_net, sw_net, F_cond = seb_terms(T_ice, Rl, Rs, albedo, G, t_base)
    if F_cond_column is not None:
        F_cond = np.where(np.isfinite(F_cond_column), F_cond_column, F_cond)
    if t_ice_method not in ("seb_steady", "diffusion"):
        seb_res = F_cond + lw_net + sw_net - fl["hsb"] - fl["hlb"]  # implied imbalance

    # Form drag over the MIZ (momentum only): neutral 10 m coefficient x U10N^2
    CD10nf = lupkes_form_drag(A_ice, zo_w)
    tau_form = rhoa * CD10nf * u10n_w**2

    A_w = 1.0 - A_ice
    tau_blend = A_w * tau_w + A_ice * fl["tau"] + tau_form
    hsb_blend = A_w * hsb_w + A_ice * fl["hsb"]
    hlb_blend = A_w * hlb_w + A_ice * fl["hlb"]

    diag = {
        k: fl[k]
        for k in ("zeta", "Rib", "zeta_limited", "qc_rib", "qc_rstar", "iter_dH", "iter_dLE",
                  "not_converged")
    }  # fmt: skip
    diag["t_ice_filled"] = t_ice_filled
    used = [u, t, rh, P, grav, A_ice, T_ice]
    if t_ice_method in ("seb_steady", "diffusion"):
        used += [Rl, Rs, h_ice, h_snow, albedo]  # T_ice is NaN where a slab value was
    diag["missing_input"] = ~np.all([np.isfinite(x) for x in used], axis=0)
    cols = [
        fl["tau"], fl["hsb"], fl["hlb"], fl["Cdn"], fl["Chn"], fl["Cen"],
        fl["R_star"], fl["z0T"], fl["z0q"],
        tau_blend, hsb_blend, hlb_blend,
        T_ice, F_cond, lw_net, sw_net, seb_res, tau_form,
    ]  # fmt: skip
    return cols, diag
