"""Steady surface energy balance (SEB) at the ice surface.

Sign conventions: radiative terms positive downward; conductive flux positive upward
(toward the surface); turbulent fluxes positive upward.
"""

from __future__ import annotations

import numpy as np

from ..most import ice_surface_layer
from ..radiation import EPS_SFC, SIGMA
from ..thermo import T_FREEZE_SW, TK0
from .properties import ALBEDO_ICE, ALBEDO_SNOW, K_ICE_FRESH, K_SNOW


def effective_conductance(h_ice, h_snow, k_ice=K_ICE_FRESH, k_snow=K_SNOW):
    """Conductance (W m-2 K-1) of a snow layer over an ice slab in series.

    G = 1 / (h_i / k_i + h_s / k_s),   h_i bounded below by 1 mm, h_s by 0.
    """
    h_ice = np.asarray(h_ice, dtype=float)
    h_snow = np.asarray(h_snow, dtype=float)
    return 1.0 / (np.maximum(h_ice, 1e-3) / k_ice + np.maximum(h_snow, 0.0) / k_snow)


def surface_albedo(h_snow):
    """Winter albedo: 0.85 (snow) if h_snow > 2 cm, else 0.65 (bare ice); NaN if h_snow is."""
    h_snow = np.asarray(h_snow, dtype=float)
    albedo = np.where(h_snow > 0.02, ALBEDO_SNOW, ALBEDO_ICE)
    return np.where(np.isnan(h_snow), np.nan, albedo)


def seb_terms(t_ice, Rl, Rs, albedo, G, t_base):
    """Radiative and conductive SEB terms at the ice surface (W m-2).

    Returns (LW_net, SW_net, F_cond) with T_s = min(t_ice, 0 C):
    LW_net = eps Rl - eps sigma T_s^4, SW_net = (1 - albedo) Rs, F_cond = G (t_base - T_s).
    """
    TsK = np.minimum(t_ice, 0.0) + TK0
    lw_net = EPS_SFC * Rl - EPS_SFC * SIGMA * TsK**4
    sw_net = (1.0 - albedo) * Rs
    F_cond = G * (t_base - np.minimum(t_ice, 0.0))
    return lw_net, sw_net, F_cond


def solve_seb_steady(flux_kwargs, Rl, Rs, albedo, G, t_base, t_air, n_bisect=40):
    """Solve the steady SEB for the ice surface temperature by vectorised bisection.

        f(Ts) = F_cond(Ts) + LW_net(Ts) + SW_net - H(Ts) - LE(Ts) = 0

    H and LE come from the full ice surface-layer solver,
    :func:`coare_seaice.most.ice_surface_layer`, called with ``flux_kwargs`` (all its
    arguments except ``t_ice``); ``flux_kwargs`` may instead be a function of Ts that
    returns a dict with ``hsb`` and ``hlb``. The search interval is [min(t_air, -1.8) - 40, 0] C.
    f decreases with Ts. If f(0) > 0 the surface is melting: Ts = 0 C and f(0) is
    returned as the melt flux.

    Returns (Ts, residual); Ts is NaN where f(0) is not finite.
    """
    N = t_air.size
    lo = np.minimum(t_air, T_FREEZE_SW) - 40.0
    hi = np.zeros(N)
    flux = flux_kwargs if callable(flux_kwargs) else None

    def f(Ts):
        fl = flux(Ts) if flux else ice_surface_layer(t_ice=Ts, **flux_kwargs)
        lw, sw, fc = seb_terms(Ts, Rl, Rs, albedo, G, t_base)
        return fc + lw + sw - fl["hsb"] - fl["hlb"]

    f_hi = f(hi)
    melt = f_hi > 0
    for _ in range(n_bisect):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        pos = fm > 0
        lo = np.where(pos, mid, lo)
        hi = np.where(pos, hi, mid)
    Ts = np.where(melt, 0.0, 0.5 * (lo + hi))
    res = np.where(melt, f_hi, f(Ts))
    bad = ~np.isfinite(f_hi)
    Ts[bad] = np.nan
    return Ts, res
