"""coare-seaice: COARE bulk fluxes extended to sea ice and the marginal ice zone.

import coare_seaice as cs

res = cs.ice_fluxes(u, t, rh, p, sst, sic, zu=36, zt=22, zq=22, lat=lat,
                    lw_down=lw, sw_down=sw)
res.blended.hsb                     # area-mean sensible heat flux, positive upward
res.tiles["thick_ice"].t_surface    # per-tile results
res.to_xarray()                     # everything, labelled, with units
"""

from importlib.metadata import PackageNotFoundError, version

from .api import ice_fluxes
from .most import ANDREAS2010_LIGHT_WIND, SurfaceLayerOptions
from .results import FluxResult
from .surface_temperature.methods import (
    SEB,
    AirCapped,
    AirFreezing,
    Column,
    Empirical,
    Prescribed,
)
from .tiles import Tile, three_tile, two_tile

try:
    __version__ = version("coare-seaice")
except PackageNotFoundError:  # running from a source tree that is not installed
    __version__ = "0.0.0"

__all__ = [
    "ANDREAS2010_LIGHT_WIND",
    "SEB",
    "AirCapped",
    "AirFreezing",
    "Column",
    "Empirical",
    "FluxResult",
    "Prescribed",
    "SurfaceLayerOptions",
    "Tile",
    "__version__",
    "ice_fluxes",
    "three_tile",
    "two_tile",
]
