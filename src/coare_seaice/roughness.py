"""Roughness lengths over sea ice.

Scalar (heat and moisture) roughness after the Andreas (1987) surface-renewal model:

    ln(z0s / z0) = b0 + b1 ln R* + b2 (ln R*)^2,    R* = z0 u* / nu,

with coefficients from Andreas (1987) Table 1 for three regimes of the roughness
Reynolds number R* (smooth, transitional, rough).
"""

from __future__ import annotations

import numpy as np

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
