"""Monin-Obukhov surface-layer solver for the ice fraction.

An iteration for u*, theta* and q* over ice, independent of the open-water COARE loop,
using the COARE 3.5 stability functions and either the Andreas (1987) scalar roughness
('blended_a87') or a fixed ratio z0T = z0q = scalar_ratio z0 ('simple').
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .roughness import NU_AIR_V2, a87_scalar_roughness
from .stability import psit_26, psiu_26
from .thermo import LAPSE_DRY, LS_ICE, TDK, air_viscosity, qsat_ice

ZETA_MAX_DEFAULT = 10.0  # limit on |zeta| = |zu / L| in the iteration (the v2 value)
RIB_CRITICAL = 0.2  # Grachev et al. (2007): critical bulk Richardson number, Ri_B,cr ~ 0.2


@dataclass(frozen=True)
class SurfaceLayerOptions:
    """Numerical and physical options of :func:`ice_surface_layer`.

    zeta_max : the stability parameter zeta = zu / L is limited to [-zeta_max, zeta_max]
        while iterating. Without a limit the iteration has no solution in some very
        stable, light-wind cases (the turbulent flux collapses). The limit is a numerical
        device, not physics: where it acts the result is set by the limit, and the output
        ``zeta_limited`` is True. The default, 10, is the v2 value.
    rib_critical : samples with a bulk Richardson number above this are flagged
        (``qc_rib``). Grachev et al. (2007) find Ri_B,cr ~ 0.2 at SHEBA and state that
        their profile functions are valid for Ri_B < Ri_B,cr (text after their Eq. 8 and
        Eq. 13); zeta ~ 10 corresponds on average to Ri_B ~ 0.2 there (their Fig. 3).
    viscosity : "coare" (default) uses the temperature-dependent kinematic viscosity of
        air (:func:`coare_seaice.thermo.air_viscosity`, at the air temperature) in the
        roughness Reynolds number R* = z0 u* / nu; "constant" uses the v2 value 1.5e-5.
    """

    zeta_max: float = ZETA_MAX_DEFAULT
    rib_critical: float = RIB_CRITICAL
    viscosity: str = "coare"

    def __post_init__(self):
        if self.viscosity not in ("coare", "constant"):
            raise ValueError("viscosity must be 'coare' or 'constant'")


DEFAULT_OPTIONS = SurfaceLayerOptions()
# The v2 behaviour, kept so that the port can still be tested against the stored v2
# outputs: every audit fix that changes results is switched off here.
V2_OPTIONS = SurfaceLayerOptions(viscosity="constant")


def bulk_richardson(theta_air, t_sfc, q_air, q_sfc, u, z, ta, grav):
    """Bulk Richardson number between the surface and height z (Grachev et al. 2007, Eq. 8).

        Ri_B = (g z / theta_v) (d_theta + 0.61 theta_v d_q) / U^2,

    with d_theta and d_q the air-minus-surface differences in potential temperature (K)
    and specific humidity (kg/kg), so Ri_B > 0 is stable. theta_v is taken as
    ta (1 + 0.61 q_air) with ta in K; U is the measured wind at z = zu with no floor, so
    Ri_B is large (or infinite) in calm conditions. Temperature and humidity measured at
    other heights are used as they are.
    """
    thv = ta * (1.0 + 0.61 * q_air)
    with np.errstate(divide="ignore", invalid="ignore"):
        return grav * z / thv * ((theta_air - t_sfc) + 0.61 * thv * (q_air - q_sfc)) / u**2


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
    options=DEFAULT_OPTIONS,
):
    """Turbulent fluxes over the ice fraction (v2 ``_ice_fluxes``).

    u (m/s at zu), t (deg C at zt), Q_air (kg/kg at zq), t_ice (deg C, capped at 0 C),
    P (mb), grav (m s-2), zu/zt/zq (m, arrays), rhoa (kg m-3), cpa, von, fdg constants,
    ta = t + 273.16 (K). z0_ice is a scalar momentum roughness (m).

    Returns a dict of arrays: tau, hsb, hlb (fluxes, heat positive upward), Cdn, Chn,
    Cen (neutral 10 m coefficients x 1e3), R_star, z0T, z0q, h_s and h_l (bulk
    conductances: H = h_s (Ts - theta), LE = h_l (qs - q)), Qs (kg/kg) and L (m).

    Stability diagnostics and flags (see :class:`SurfaceLayerOptions`): zeta (= zu / L as
    used in the final pass), zeta_limited (True where the limit on |zeta| acted in the
    last iteration), Rib (bulk Richardson number) and qc_rib (True where
    Rib > options.rib_critical: outside the validated range of the stability functions).
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

    nu = air_viscosity(ta - TDK) if options.viscosity == "coare" else NU_AIR_V2

    def scalar_rough(usr):
        if ice_scheme == "blended_a87":
            return a87_scalar_roughness(z0_ice, usr, nu)
        return scalar_ratio * z0_ice, scalar_ratio * z0_ice, z0_ice * usr / nu

    # AUDIT F6: fixed number of iterations with no convergence diagnostic (task F6).
    for _ in range(nits):
        z0T_i, z0q_i, R_star = scalar_rough(usr_i)
        tvsr_i = tsr_i + 0.61 * ta * qsr_i
        zet_raw = von * grav * zu / ta * tvsr_i / np.maximum(usr_i**2, 1e-6)
        zet_i = np.clip(zet_raw, -options.zeta_max, options.zeta_max)
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
    rib = bulk_richardson(t + LAPSE_DRY * zt, t_sfc, Q_air, Qs_ice, u, zu, ta, grav)
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
        zeta=zet_i,
        zeta_limited=np.abs(zet_raw) > options.zeta_max,
        Rib=rib,
        qc_rib=rib > options.rib_critical,
    )
