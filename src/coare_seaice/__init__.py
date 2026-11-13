"""coare-seaice: COARE bulk fluxes extended to sea ice and the marginal ice zone."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("coare-seaice")
except PackageNotFoundError:  # running from a source tree that is not installed
    __version__ = "0.0.0"

__all__ = ["__version__"]
