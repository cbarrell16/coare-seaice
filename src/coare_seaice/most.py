"""Monin-Obukhov surface-layer solver for the ice fraction.

An iteration for u*, theta* and q* over ice, independent of the open-water COARE loop,
using the COARE 3.5 stability functions and either the Andreas (1987) scalar roughness
('blended_a87') or a fixed ratio z0T = z0q = scalar_ratio z0 ('simple').
"""

from __future__ import annotations

import numpy as np

from .roughness import NU_AIR, a87_scalar_roughness
from .stability import psit_26, psiu_26
from .thermo import LAPSE_DRY, LS_ICE, qsat_ice

# AUDIT F1: hard-coded stability clip |zeta| <= 10 with no flag (to be changed in task F1).
ZETA_CLIP = 10.0


def ice_surface_layer(
    u,
    t,
    Q_air,
    t_ice,
    z0_ice,
    P,
    grav,
    zu,
    zt,
    zq,
    rhoa,
    cpa,
    von,
    fdg,
    ta,
    ice_scheme,
    scalar_ratio,
    nits=10,
):
    """Turbulent fluxes over the ice fraction (v2 ``_ice_fluxes``).

    u (m/s at zu), t (deg C at zt), Q_air (kg/kg at zq), t_ice (deg C, capped at 0 C),
    P (mb), grav (m s-2), zu/zt/zq (m, arrays), rhoa (kg m-3), cpa, von, fdg constants,
    ta = t + 273.16 (K). z0_ice is a scalar momentum roughness (m).

    Returns a dict of arrays: tau, hsb, hlb (fluxes, heat positive upward), Cdn, Chn,
    Cen (neutral 10 m coefficients x 1e3), R_star, z0T, z0q, h_s and h_l (bulk
    conductances: H = h_s (Ts - theta), LE = h_l (qs - q)), Qs (kg/kg) and L (m).
    """
    N = u.size
    z0_ice = z0_ice * np.ones(N)

    t_sfc = np.minimum(t_ice, 0.0)
    Qs_ice = qsat_ice(t_sfc, P)
    dt_ice = t_sfc - t - LAPSE_DRY * zt  # potential temperature difference
    dq_ice = Qs_ice - Q_air

    # AUDIT F10: no gustiness (unstable) or windless (stable) term over ice; only a 0.5 m/s
    # floor on the wind (to be changed in task F10).
    ut = np.maximum(u, 0.5)
    usr_i = np.maximum(ut * von / np.log(zu / z0_ice), 0.001)
    tsr_i = np.zeros(N)
    qsr_i = np.zeros(N)
    L_i = np.full(N, 1e6)

    def scalar_rough(usr):
        if ice_scheme == "blended_a87":
            return a87_scalar_roughness(z0_ice, usr)
        # AUDIT F2: constant air viscosity in R* (to be changed in task F2).
        return scalar_ratio * z0_ice, scalar_ratio * z0_ice, z0_ice * usr / NU_AIR

    # AUDIT F6: fixed number of iterations with no convergence diagnostic (task F6).
    for _ in range(nits):
        z0T_i, z0q_i, R_star = scalar_rough(usr_i)
        tvsr_i = tsr_i + 0.61 * ta * qsr_i
        zet_i = np.clip(
            von * grav * zu / ta * tvsr_i / np.maximum(usr_i**2, 1e-6), -ZETA_CLIP, ZETA_CLIP
        )
        safe = np.abs(zet_i) > 1e-6
        L_i = np.where(safe, zu / np.where(safe, zet_i, 1.0), 1e6)
        psi_m = psiu_26(zu / L_i)
        psi_t = psit_26(zt / L_i)
        psi_q = psit_26(zq / L_i)
        denom_m = np.maximum(np.log(zu / z0_ice) - psi_m, 0.1)
        denom_t = np.log(zt / np.maximum(z0T_i, 1e-9)) - psi_t
        denom_q = np.log(zq / np.maximum(z0q_i, 1e-9)) - psi_q
        denom_t = np.where(np.abs(denom_t) > 0.1, denom_t, 0.1)
        denom_q = np.where(np.abs(denom_q) > 0.1, denom_q, 0.1)
        usr_i = ut * von / denom_m
        tsr_i = -dt_ice * von * fdg / denom_t
        qsr_i = -dq_ice * von * fdg / denom_q

    # Final pass with the converged u* and L.
    # AUDIT F6: u* is not updated in this final pass (to be changed in task F6).
    z0T_i, z0q_i, R_star = scalar_rough(usr_i)
    denom_t = np.log(zt / np.maximum(z0T_i, 1e-9)) - psit_26(zt / L_i)
    denom_q = np.log(zq / np.maximum(z0q_i, 1e-9)) - psit_26(zq / L_i)
    denom_t = np.where(np.abs(denom_t) > 0.1, denom_t, 0.1)
    denom_q = np.where(np.abs(denom_q) > 0.1, denom_q, 0.1)
    Ct_i = von * fdg / denom_t * von / np.maximum(np.log(zu / z0_ice) - psiu_26(zu / L_i), 0.1)
    Cq_i = von * fdg / denom_q * von / np.maximum(np.log(zu / z0_ice) - psiu_26(zu / L_i), 0.1)
    tsr_i = -dt_ice * von * fdg / denom_t
    qsr_i = -dq_ice * von * fdg / denom_q

    tau_ice = rhoa * usr_i**2
    hsb_ice = -rhoa * cpa * usr_i * tsr_i
    hlb_ice = -rhoa * LS_ICE * usr_i * qsr_i

    z0T_s, z0q_s = np.maximum(z0T_i, 1e-9), np.maximum(z0q_i, 1e-9)
    lz = np.log(10.0 / z0_ice)
    return dict(
        tau=tau_ice,
        hsb=hsb_ice,
        hlb=hlb_ice,
        Cdn=1000 * von**2 / lz**2,
        Chn=1000 * von**2 * fdg / (lz * np.log(10.0 / z0T_s)),
        Cen=1000 * von**2 * fdg / (lz * np.log(10.0 / z0q_s)),
        R_star=R_star,
        z0T=z0T_i,
        z0q=z0q_i,
        h_s=rhoa * cpa * Ct_i * ut,
        h_l=rhoa * LS_ICE * Cq_i * ut,
        Qs=Qs_ice,
        L=L_i,
    )
