"""The v2 interface: ``coare35vn_ice(...)`` returning the 56-column array.

Same call signature, defaults and output columns as ``coare35vn_ice.py`` version 2.0.0
(Nov 2026). Open-water columns 0-37 come from ``pycoare`` (through
:mod:`coare_seaice.openwater`); ice columns 38-55 from :func:`coare_seaice.ice_fluxes`
with two tiles, of which this is a thin wrapper.

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

from types import SimpleNamespace

import numpy as np

from .api import ice_fluxes
from .openwater import open_water
from .results import _LEGACY_ICE
from .roughness import ICE_Z0_DEFAULT, SCALAR_RATIO
from .surface_temperature.methods import (
    SEB,
    AirCapped,
    AirFreezing,
    Column,
    Empirical,
    Prescribed,
)
from .surface_temperature.properties import K_SNOW
from .surface_temperature.seb import surface_albedo
from .thermo import T_FREEZE_SW
from .tiles import two_tile

__all__ = ["OUTPUT_COLUMNS", "T_ICE_METHODS", "coare35vn_ice"]

ICE_T_DEFAULT = T_FREEZE_SW  # deg C, ice base temperature
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
        u=u, t=t, rh=rh, ts=ts, P=P, Rs=Rs, Rl=Rl, lat=lat, zi=zi, zu=zu_a, zt=zt_a, zq=zq_a,
        ow=ow_u, A_ice=A_ice, z0_ice=z0_ice_val, ice_scheme=ice_scheme,
        scalar_ratio=scalar_ratio, t_ice_method=t_ice_method, t_ice=t_ice,
        t_ice_fallback=t_ice_fallback, h_ice=h_ice_v, h_snow=h_snow_v, k_snow=k_snow,
        albedo=alb_v, t_base=t_base, t_ice_offset=t_ice_offset, time_s=time_s,
        diffusion_kwargs=diffusion_kwargs, rh_reference=rh_reference, verbose=verbose,
        options=surface_layer_options, v2_compat=_v2_compat, jcool=jcool,
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


def _method(name, *, t_ice=None, fallback=None, offset=2.0, diffusion_kwargs=None, n=None):
    """The surface temperature method object for a v2 method name."""
    if name == "air_freezing":
        return AirFreezing(ICE_T_DEFAULT)
    if name == "air_capped":
        return AirCapped(0.0)
    if name == "empirical":
        return Empirical(offset=offset)
    if name == "seb_steady":
        return SEB()
    if name == "prescribed":
        fb = None if fallback is None else _method(fallback, offset=offset)
        return Prescribed(_per_sample(t_ice, n, "t_ice"), fallback=fb)
    if name == "diffusion":
        dk = dict(diffusion_kwargs or {})
        run_kw = {k: dk.pop(k) for k in ("sic_reset", "max_gap_s") if k in dk}
        return Column(column_kwargs=dk, **run_kw)
    raise ValueError(f"unsupported fallback {name!r}")


def open_water_stub(*, rhoa, q_air, zo, u_n_rf, tau, hsb, hlb):
    """Stand-in for the pycoare result, holding only what the ice code and the blend use.

    For tests that feed stored open-water values in; everything else is NaN.
    """
    nan = np.full(np.shape(tau), np.nan)
    return SimpleNamespace(
        rhoa=rhoa, q_air=q_air, zo=zo, u_n_rf=u_n_rf, tau=tau, hsb=hsb, hlb=hlb, usr=nan,
        zot=nan, zoq=nan, cdn_rf=nan, chn_rf=nan, cen_rf=nan, obukhov_length=nan, zet=nan,
        ut=nan, q_sea=nan, dter=0.0, dqer=0.0, jcool=0,
    )  # fmt: skip


def _ice_columns(
    *, u, t, rh, ts, P, Rs, Rl, lat, zi, zu, zt, zq, ow, A_ice, z0_ice, ice_scheme,
    scalar_ratio, t_ice_method, t_ice, t_ice_fallback, h_ice, h_snow, k_snow, albedo, t_base,
    t_ice_offset, time_s, diffusion_kwargs, rh_reference, verbose, options=None,
    v2_compat=False, jcool=1,
):  # fmt: skip
    """The 18 v2 ice columns (38-55) and the diagnostics, from :func:`ice_fluxes`.

    ``ow`` is the open-water result at a 10 m reference height (or an
    :func:`open_water_stub`). Returns (columns, diagnostics).
    """
    N = u.size
    method = _method(t_ice_method, t_ice=t_ice, fallback=t_ice_fallback, offset=t_ice_offset,
                     diffusion_kwargs=diffusion_kwargs, n=N)  # fmt: skip
    if t_ice_method == "prescribed" and verbose and t_ice_fallback is not None:
        miss = int((~np.isfinite(method.values)).sum())
        if miss:
            print(f"  prescribed T_ice missing in {miss} of {N} steps -> {t_ice_fallback!r}")
    tiles = two_tile(
        z0=z0_ice, scalar_scheme="a87" if ice_scheme == "blended_a87" else "ratio",
        scalar_ratio=scalar_ratio, h_ice=h_ice, h_snow=h_snow, k_snow=k_snow, albedo=albedo,
        t_base=t_base,
    )  # fmt: skip
    try:
        res = ice_fluxes(
            u, t, rh, P, ts, A_ice, zu=zu, zt=zt, zq=zq, lat=lat, lw_down=Rl, sw_down=Rs,
            tiles=tiles, t_ice=method, options=options, zi=zi, jcool=jcool, time=time_s,
            rh_reference=rh_reference, _open_water=ow, _v2_compat=v2_compat,
        )  # fmt: skip
    except ValueError as err:  # keep the v2 wording for a wrong-length input
        raise ValueError(f"coare35vn_ice: {err}") from err
    ice = res.ice
    diag = {
        "zeta": ice.zeta, "Rib": ice.rib, "zeta_limited": ice.qc["zeta_limited"],
        "qc_rib": ice.qc["rib_above_critical"], "qc_rstar": ice.qc["rstar_out_of_range"],
        "iter_dH": ice.extra["iter_dH"], "iter_dLE": ice.extra["iter_dLE"],
        "not_converged": ice.qc["not_converged"], "t_ice_filled": ice.qc["t_surface_filled"],
        "missing_input": ice.qc["missing_input"],
    }  # fmt: skip
    src = {"ice": ice, "blended": res.blended}
    cols = [getattr(src[s_], f) * scale for _, s_, f, scale in _LEGACY_ICE]
    return cols, diag
