"""Form drag over the marginal ice zone (momentum only).

Lupkes et al. (2012) form-drag parametrisation with the E2016B parameter set of
Elvidge et al. (2016) Table 1, as used by Renfrew et al. (2019).
"""

from __future__ import annotations

import numpy as np

# E2016B parameters (Elvidge et al. 2016, Table 1)
CE_FORM = 0.10  # effective resistance coefficient (Lupkes et al. 2012: 0.3)
BETA_FD = 0.2  # floe-size exponent (Lupkes et al. 2012: 1.0)
S_SHELT = 0.5  # sheltering constant
H_MAX = 0.534  # freeboard at full ice cover (m)
H_MIN = 0.286  # freeboard at the MIZ edge (m)
D_MIN = 8.0  # minimum floe diameter (m)
D_MAX = 300.0  # maximum floe diameter (m)


def lupkes_form_drag(A_ice, z0w):
    """Form-drag part of the neutral 10 m drag coefficient over the MIZ.

    Lupkes et al. (2012) Eq. 22 / Elvidge et al. (2016) Eq. 11, E2016B parameters:

        C_DN10f = A (h_f / D_i) S_c^2 (c_e / 2) [ln(h_f / z0w) / ln(10 / z0w)]^2

    A_ice is the ice concentration (0-1) and z0w the open-water roughness length (m).
    Returns 0 where A < 0.01 and NaN where A is not finite. A is clipped to
    [1e-6, 1 - 1e-6], so at A = 1 the result is very small but not exactly zero.
    """
    # AUDIT F9: citation to be corrected to Lupkes et al. (2012) Eq. 51 (task F9).
    A_in = np.asarray(A_ice, dtype=float)
    A = np.clip(A_in, 1e-6, 1.0 - 1e-6)
    A_star = 1.0 / (1.0 - (D_MIN / D_MAX) ** (1.0 / BETA_FD))
    Di = D_MIN * (A_star / (A_star - A)) ** BETA_FD
    hf = H_MAX * A + H_MIN * (1.0 - A)
    Dw = Di * (1.0 - np.sqrt(A)) / np.sqrt(A)
    Sc = 1.0 - np.exp(-S_SHELT * Dw / hf)
    hf_safe = np.maximum(hf, z0w * 1.01)
    log_ratio = np.clip(np.log(hf_safe / z0w) / np.log(10.0 / z0w), 0.0, 1.0)
    CD10nf = A * (hf / Di) * Sc**2 * (CE_FORM / 2.0) * log_ratio**2
    CD10nf = np.where(A_in < 0.01, 0.0, CD10nf)
    CD10nf = np.where(np.isfinite(A_in), CD10nf, np.nan)
    return CD10nf
