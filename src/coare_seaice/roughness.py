"""Roughness lengths over sea ice.

Scalar (heat and moisture) roughness after the Andreas (1987) surface-renewal model:

    ln(z0s / z0) = b0 + b1 ln R* + b2 (ln R*)^2,    R* = z0 u* / nu,

with coefficients from Andreas (1987) Table 1 for three regimes of the roughness
Reynolds number R* (smooth, transitional, rough).

Momentum roughness of ice is regional and rarely measured: it is an input with an
uncertainty range, not a constant. :data:`ICE_ROUGHNESS` lists published values by name,
and :func:`z0_andreas2010_winter` is a formula that depends on the friction velocity.
"""

from __future__ import annotations

import numpy as np

VON = 0.4

ICE_Z0_DEFAULT = 5.0e-4  # m, MetUM operational z0 over ice (C_DN10i ~ 1.6e-3, Lupkes 2012)
SCALAR_RATIO = 0.2  # z0T/z0 = z0q/z0 for the 'simple' scheme (MetUM)

# Constant kinematic viscosity used by v2 (m2 s-1). It corresponds to about +19 C; the
# solver now uses the temperature-dependent value, :func:`coare_seaice.thermo.air_viscosity`.
NU_AIR_V2 = 1.5e-5

# Andreas (1987) Table 1, rows: smooth (R* <= 0.135), transitional (0.135 < R* < 2.5),
# rough (2.5 <= R* <= 1000); columns: b0, b1, b2.
A87_ZOT = np.array(
    [
        [1.250, 0.000, 0.000],
        [0.149, -0.550, 0.000],
        [0.317, -0.565, -0.183],
    ]
)
A87_ZOQ = np.array(
    [
        [1.610, 0.000, 0.000],
        [0.351, -0.628, 0.000],
        [0.396, -0.512, -0.180],
    ]
)
R_SMOOTH = 0.135
R_ROUGH = 2.5
R_MAX = 1000.0  # upper end of the range covered by the Andreas (1987) rough-regime fit


def roughness_reynolds(z0, u_star, nu):
    """Roughness Reynolds number R* = z0 u* / nu (no lower bound)."""
    return z0 * u_star / nu


def a87_scalar_roughness(z0_i, u_star_i, nu=NU_AIR_V2, guard=True):
    """Andreas (1987) scalar roughness lengths over ice.

    Parameters are arrays of the same shape: momentum roughness z0_i (m), friction
    velocity u_star_i (m/s) and kinematic viscosity of air nu (m2 s-1, array or number;
    use :func:`coare_seaice.thermo.air_viscosity`; the default is the v2 constant).
    Returns (z0T, z0q, R*), R* bounded below by 1e-6.

    The rough-regime polynomial is fitted for 2.5 <= R* <= 1000. With ``guard`` (default)
    it is not extrapolated: for R* > 1000 the ratios z0T/z0 and z0q/z0 are held at their
    R* = 1000 values (4.5e-6 and 8.0e-6). The returned R* is the actual value, so callers
    can flag R* > ``R_MAX``. ``guard=False`` extrapolates as v2 did.
    """
    R_star = np.maximum(z0_i * u_star_i / nu, 1e-6)
    lnR = np.log(R_star)
    bt = A87_ZOT[0, 0] + 0 * lnR
    bq = A87_ZOQ[0, 0] + 0 * lnR
    k = (R_star > R_SMOOTH) & (R_star < R_ROUGH)
    bt[k] = A87_ZOT[1, 0] + A87_ZOT[1, 1] * lnR[k]
    bq[k] = A87_ZOQ[1, 0] + A87_ZOQ[1, 1] * lnR[k]
    k = R_star >= R_ROUGH
    lnr = np.minimum(lnR[k], np.log(R_MAX)) if guard else lnR[k]
    bt[k] = A87_ZOT[2, 0] + A87_ZOT[2, 1] * lnr + A87_ZOT[2, 2] * lnr**2
    bq[k] = A87_ZOQ[2, 0] + A87_ZOQ[2, 1] * lnr + A87_ZOQ[2, 2] * lnr**2
    return z0_i * np.exp(bt), z0_i * np.exp(bq), R_star


def z0_from_cdn10(cdn10):
    """Roughness length (m) for a neutral 10 m drag coefficient: z0 = 10 exp(-k / sqrt(C))."""
    return 10.0 * np.exp(-VON / np.sqrt(np.asarray(cdn10, dtype=float)))


def cdn10_from_z0(z0):
    """Neutral 10 m drag coefficient for a roughness length (m): C = (k / ln(10 / z0))^2."""
    return (VON / np.log(10.0 / np.asarray(z0, dtype=float))) ** 2


def z0_andreas2010_winter(u_star, nu):
    """Roughness of winter sea ice from the friction velocity (Andreas et al. 2010, J.
    Hydrometeorol. 11, 87-104, Eq. 5.1):

        z0 = 0.135 nu / u* + 2.30e-4 tanh^3(13 u*)        (m, with mks units)

    Fitted by eye to SHEBA winter data (their Fig. 2), so SHEBA is not independent
    evidence for it. The first term is aerodynamically smooth flow; z0 levels off at
    2.3e-4 m (C_DN10 = 1.4e-3) for u* above about 0.2 m s-1.
    """
    u_star = np.maximum(np.asarray(u_star, dtype=float), 1e-4)
    return 0.135 * nu / u_star + 2.30e-4 * np.tanh(13.0 * u_star) ** 3


# Published neutral 10 m drag coefficients over (nearly) complete ice cover, by name.
# cdn10 is the median (or the stated setting); iqr the interquartile range where given.
_S22 = "Srivastava et al. (2022, ACP discussion paper acp-2021-705), Table 2"
_R19 = "Renfrew et al. (2019), Table 1"
ICE_ROUGHNESS = {
    "metum_smooth": dict(z0=0.5e-3, source=f"{_R19}: 'Smooth' setting, z0_ice = 0.5e-3 m"),
    "metum_rough": dict(z0=3.0e-3, source=f"{_R19}: 'Rough' setting, z0_ice = 3e-3 m"),
    "eastern_arctic": dict(cdn10=2.1e-3, iqr=(1.3e-3, 2.8e-3), source=f"{_S22}: ACSE ship"),
    "central_arctic": dict(cdn10=2.2e-3, iqr=(1.6e-3, 2.7e-3), source=f"{_S22}: AO2016 ship"),
    "iceland_sea": dict(cdn10=3.4e-3, iqr=(2.5e-3, 4.2e-3),
                        source=f"{_S22}: Elvidge et al. (2021), aircraft (IGP)"),
    "barents_sea_broken_floes": dict(cdn10=2.6e-3, iqr=(2.4e-3, 3.9e-3),
                                     source=f"{_S22}: Elvidge et al. (2016), aircraft"),
    "fram_strait_large_flat_floes": dict(cdn10=0.9e-3, iqr=(0.4e-3, 2.1e-3),
                                         source=f"{_S22}: Elvidge et al. (2016), aircraft"),
}  # fmt: skip
DYNAMIC_ROUGHNESS = {"andreas2010_winter": z0_andreas2010_winter}


def resolve_ice_roughness(spec):
    """Turn a roughness setting into what the solver takes.

    A number or array (m) is returned as it is; a function ``f(u_star, nu)`` too; a name
    from :data:`ICE_ROUGHNESS` gives its z0 (from the drag coefficient where that is what
    is published); "andreas2010_winter" gives :func:`z0_andreas2010_winter`.
    """
    if callable(spec):
        return spec
    if isinstance(spec, str):
        if spec in DYNAMIC_ROUGHNESS:
            return DYNAMIC_ROUGHNESS[spec]
        if spec not in ICE_ROUGHNESS:
            names = sorted([*ICE_ROUGHNESS, *DYNAMIC_ROUGHNESS])
            raise ValueError(f"unknown ice roughness {spec!r}; use a value in m or one of {names}")
        entry = ICE_ROUGHNESS[spec]
        return float(entry["z0"] if "z0" in entry else z0_from_cdn10(entry["cdn10"]))
    z0 = np.asarray(spec, dtype=float)
    if np.any(z0 <= 0):
        raise ValueError("roughness length must be positive")
    return z0
