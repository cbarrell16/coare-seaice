"""Ice surface temperature methods.

- :mod:`.simple`: assumptions with no energy balance (air_freezing, air_capped, empirical).
- :mod:`.seb`: steady surface energy balance with conduction through snow and ice.
- :mod:`.column`: time-dependent snow and ice conduction model (``IceColumn``),
  experimental.
- :mod:`.methods`: the same methods as classes behind one interface (``solve(ctx)``), used
  by :func:`coare_seaice.ice_fluxes`.
"""

from .column import IceColumn
from .methods import (
    SEB,
    AirCapped,
    AirFreezing,
    Column,
    Empirical,
    Prescribed,
    SurfaceTemperature,
    TileContext,
    method_from_name,
)
from .seb import effective_conductance, seb_terms, solve_seb_steady, surface_albedo
from .simple import air_capped, air_freezing, empirical

__all__ = [
    "SEB",
    "AirCapped",
    "AirFreezing",
    "Column",
    "Empirical",
    "IceColumn",
    "Prescribed",
    "SurfaceTemperature",
    "TileContext",
    "method_from_name",
    "air_capped",
    "air_freezing",
    "effective_conductance",
    "empirical",
    "seb_terms",
    "solve_seb_steady",
    "surface_albedo",
]
