"""The public entry point: :func:`ice_fluxes`."""

from __future__ import annotations

import numpy as np

from .engine import Forcing, compute
from .surface_temperature.methods import Prescribed, method_from_name
from .tiles import two_tile

__all__ = ["ice_fluxes"]


def _array(x, n, name):
    """x as a float array of length n (a copy); size 1 is repeated."""
    x = np.array(x, dtype=float)
    if x.size == 1:
        return np.full(n, float(x.reshape(-1)[0]))
    if x.ndim != 1 or x.size != n:
        raise ValueError(f"ice_fluxes: {name} must be a number or have one value per sample")
    return x


def ice_fluxes(
    u,
    t,
    rh,
    p,
    sst,
    sic,
    *,
    zu,
    zt,
    zq,
    lat,
    lw_down=np.nan,
    sw_down=np.nan,
    coare="3.5",
    tiles=None,
    t_ice=None,
    z0_ice=None,
    blend="flux",
    t_surface=None,
    options=None,
    zi=600.0,
    rain=None,
    cp=None,
    sigH=None,
    jcool=1,
    salinity=None,
    time=None,
    rh_reference="water",
    _open_water=None,
    _v2_compat=False,
):
    """Turbulent fluxes over open water, sea ice and the marginal ice zone.

    Inputs are numbers or 1-D arrays with one value per sample.

    u : wind speed (m s-1) at height zu, relative to the surface.
    t : air temperature (deg C) at zt.
    rh : relative humidity (%) with respect to liquid water, at zq.
    p : surface air pressure (mb).
    sst : sea surface temperature (deg C) for the open-water tile.
    sic : sea-ice concentration (0 to 1). NaN gives NaN blended fluxes.
    zu, zt, zq : measurement heights (m), a number or one value per sample.
    lat : latitude (deg), for gravity.
    lw_down, sw_down : downwelling longwave and shortwave (W m-2). Needed by the
        energy-balance surface temperature methods and by the cool-skin correction
        over water (``jcool=1``). There is no default value: if they are missing,
        whatever depends on them is NaN. Turbulent fluxes with a non-energy-balance
        method and ``jcool=0`` do not need them.
    coare : "3.5" (default) or "3.6", the COARE version for open water (pycoare).
    tiles : a list of :class:`coare_seaice.tiles.Tile`; default ``two_tile()``.
    t_ice : surface temperature method for ice tiles that do not set their own: a
        method object (``AirFreezing()``, ``SEB(...)``, ...), a name ("air_freezing",
        "seb", ...), or an array of prescribed temperatures (deg C). Default
        ``AirFreezing()``: min(T_air, -1.8 C), the simple reference.
    z0_ice : momentum roughness (m) for ice tiles that do not set their own; default
        5e-4 m. Ice roughness is regional and rarely measured: treat it as an input with
        an uncertainty range.
    blend : "flux" (default) or "single_temperature"; see
        :func:`coare_seaice.engine.compute`. ``t_surface`` is the aggregate surface
        temperature for the latter.
    options : :class:`coare_seaice.most.SurfaceLayerOptions` (stability limits, flag
        thresholds, light-wind terms) for the ice tiles.
    zi, rain, cp, sigH, jcool, salinity : passed to COARE for open water (boundary-layer
        depth m, rain rate mm h-1, wave phase speed m s-1, significant wave height m,
        cool-skin switch, salinity psu for COARE 3.6).
    time : time in seconds, needed only by the ``Column`` method.
    rh_reference : "water" (default) or "ice" (the v1 reading of the humidity input).

    Returns a :class:`coare_seaice.results.FluxResult`: ``res.blended.hsb``,
    ``res.tiles["thick_ice"].t_surface``, ``res.qc``, ``res.to_xarray()``. Heat fluxes
    are positive upward. With the defaults the result is the simple reference: two
    tiles, ice at min(T_air, -1.8 C), z0_ice = 5e-4 m, COARE stability functions.

    Missing data: NaN in gives NaN out for everything that depends on it; flags mark,
    and never change, a value.
    """
    u = np.array(u, dtype=float).reshape(-1)
    n = u.size
    named = dict(t=t, rh=rh, p=p, sst=sst, sic=sic, zu=zu, zt=zt, zq=zq, lat=lat,
                 lw_down=lw_down, sw_down=sw_down, zi=zi)  # fmt: skip
    a = {k: _array(v, n, k) for k, v in named.items()}
    for k in ("zu", "zt", "zq"):
        if np.any(a[k] <= 0):
            raise ValueError(f"ice_fluxes: {k} must be positive")
    bad = np.isfinite(a["sic"]) & ((a["sic"] < 0) | (a["sic"] > 1))
    if bad.any():
        raise ValueError("ice_fluxes: sic must be between 0 and 1 (a fraction, not per cent)")
    optional = {k: None if v is None else _array(v, n, k)
                for k, v in dict(rain=rain, cp=cp, sigH=sigH, time=time).items()}  # fmt: skip
    F = Forcing(u=u, **a, **optional, jcool=0 if jcool == 0 else 1, salinity=salinity)

    if isinstance(t_ice, str):
        t_ice = method_from_name(t_ice)
    elif t_ice is not None and not hasattr(t_ice, "solve"):
        t_ice = Prescribed(_array(t_ice, n, "t_ice"))
    return compute(
        F, two_tile() if tiles is None else list(tiles), blend=blend, coare=coare,
        t_ice=t_ice, z0_ice=z0_ice, options=options, rh_reference=rh_reference,
        t_surface=t_surface, open_water_result=_open_water, v2_compat=_v2_compat,
    )  # fmt: skip
