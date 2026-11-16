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


def roughness_reynolds(z0, u_star, nu):
    """Roughness Reynolds number R* = z0 u* / nu (no lower bound)."""
    return z0 * u_star / nu


def a87_scalar_roughness(z0_i, u_star_i, nu=NU_AIR_V2):
    """Andreas (1987) scalar roughness lengths over ice.

    Parameters are arrays of the same shape: momentum roughness z0_i (m), friction
    velocity u_star_i (m/s) and kinematic viscosity of air nu (m2 s-1, array or number;
    use :func:`coare_seaice.thermo.air_viscosity`; the default is the v2 constant).
    Returns (z0T, z0q, R*), R* bounded below by 1e-6.
    """
    R_star = np.maximum(z0_i * u_star_i / nu, 1e-6)
    lnR = np.log(R_star)
    bt = A87_ZOT[0, 0] + 0 * lnR
    bq = A87_ZOQ[0, 0] + 0 * lnR
    k = (R_star > R_SMOOTH) & (R_star < R_ROUGH)
    bt[k] = A87_ZOT[1, 0] + A87_ZOT[1, 1] * lnR[k]
    bq[k] = A87_ZOQ[1, 0] + A87_ZOQ[1, 1] * lnR[k]
    # AUDIT F3: the rough-regime polynomial is also used beyond R* = 1000 (its stated
    # range) with no guard or flag (to be changed in task F3).
    k = R_star >= R_ROUGH
    bt[k] = A87_ZOT[2, 0] + A87_ZOT[2, 1] * lnR[k] + A87_ZOT[2, 2] * lnR[k] ** 2
    bq[k] = A87_ZOQ[2, 0] + A87_ZOQ[2, 1] * lnR[k] + A87_ZOQ[2, 2] * lnR[k] ** 2
    return z0_i * np.exp(bt), z0_i * np.exp(bq), R_star
