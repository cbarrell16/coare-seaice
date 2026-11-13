"""Radiation constants and a fallback for downwelling longwave."""

from __future__ import annotations

import numpy as np

from .thermo import TK0

SIGMA = 5.67e-8  # Stefan-Boltzmann constant, W m-2 K-4
EPS_SFC = 0.99  # emissivity of the snow / ice surface


def lw_down_konig_langlo(t_air, cloud_fraction):
    """Downwelling longwave (W m-2) after Koenig-Langlo & Augstein (1994), for polar air.

        LWd = (0.765 + 0.22 N^3) sigma T_a^4,   T_a in K, N cloud fraction (clipped to 0-1).

    Only a fallback when no observation or reanalysis value is available.
    """
    n = np.clip(np.asarray(cloud_fraction, dtype=float), 0.0, 1.0)
    ta = np.asarray(t_air, dtype=float) + TK0
    return (0.765 + 0.22 * n**3) * SIGMA * ta**4
