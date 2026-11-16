"""Ice surface temperature assumptions with no energy balance. Temperatures in deg C."""

from __future__ import annotations

import numpy as np

from ..thermo import T_FREEZE_SW


def air_freezing(t_air, t_freeze=T_FREEZE_SW):
    """T_ice = min(T_air, -1.8 C) (the v1 default and the simple reference)."""
    return np.minimum(np.asarray(t_air, dtype=float), t_freeze)


def air_capped(t_air, cap=0.0):
    """T_ice = min(T_air, 0 C): the surface follows the air but cannot exceed melting."""
    return np.minimum(np.asarray(t_air, dtype=float), cap)


def empirical(t_air, offset=2.0, func=None, **kwargs):
    """T_ice = min(T_air + dT, 0 C).

    offset : dT in K (default +2 K, a placeholder until fitted to IGP radiometer data).
    func   : optional callable giving dT = func(t_air, **kwargs) instead of ``offset``.
    """
    t_air = np.asarray(t_air, dtype=float)
    dT = func(t_air, **kwargs) if func is not None else offset
    return np.minimum(t_air + dT, 0.0)
