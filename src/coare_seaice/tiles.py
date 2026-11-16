"""Surface tiles: what each part of the surface is and how much of it there is.

A tile set is a list of :class:`Tile`. The area fractions come from the total sea-ice
concentration ``sic`` through each tile's ``fraction`` function and must sum to one.
Two sets are provided:

- :func:`two_tile`: open water and ice (the simple reference, as v2);
- :func:`three_tile`: open water, thin ice and thick ice.

**The thin-ice fraction has no routine measurement.** It is an input: for ships it is a
sensitivity parameter. The defaults of the thin-ice tile are provisional (task T10).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, replace

import numpy as np

from .surface_temperature.properties import K_SNOW
from .thermo import T_FREEZE_SW

__all__ = ["Tile", "three_tile", "tile_fractions", "two_tile"]


@dataclass(frozen=True)
class Tile:
    """One surface type.

    name : used in the results (``res.tiles[name]``).
    kind : "water" (fluxes from COARE through pycoare) or "ice".
    fraction : function of the total ice concentration giving this tile's area fraction.
    form_drag : whether the tile counts as ice for the floe-edge form drag.
    Ice tiles only:
    t_surface : a surface temperature method (see
        :mod:`coare_seaice.surface_temperature.methods`); None takes the ``t_ice``
        argument of :func:`coare_seaice.ice_fluxes`.
    z0 : momentum roughness length: a number or one value per sample (m), a name from
        :data:`coare_seaice.roughness.ICE_ROUGHNESS`, "andreas2010_winter", or a
        function of (u*, nu). None takes the ``z0_ice`` argument.
    scalar_scheme : "a87" (Andreas 1987), "ratio" (z0T = z0q = scalar_ratio * z0) or
        "prescribed" (the values ``z0t`` and ``z0q``, m; giving them selects it).
    A water tile with ``z0``, ``z0t`` and ``z0q`` all given is solved with those
    roughness lengths instead of COARE, at the sea surface temperature, with no
    cool-skin correction and no gustiness (unless the options switch it on). This is for
    anchoring the open-water transfer coefficients to observations, as Elvidge et al.
    (2021) do.
    h_ice, h_snow, k_snow, albedo, t_base : thickness (m), snow conductivity, albedo
        (None: 0.85 with more than 2 cm of snow, else 0.65) and ice-base temperature,
        used by energy-balance methods and for the reported conductive flux.
    """

    name: str
    kind: str
    fraction: Callable
    form_drag: bool = False
    t_surface: object | None = None
    z0: float | np.ndarray | str | Callable | None = None
    z0t: float | np.ndarray | None = None
    z0q: float | np.ndarray | None = None
    scalar_scheme: str = "a87"
    scalar_ratio: float = 0.2
    h_ice: float | np.ndarray = 1.0
    h_snow: float | np.ndarray = 0.2
    k_snow: float = K_SNOW
    albedo: float | np.ndarray | None = None
    t_base: float = T_FREEZE_SW

    def __post_init__(self):
        if self.kind not in ("water", "ice"):
            raise ValueError("Tile.kind must be 'water' or 'ice'")
        if (self.z0t is None) != (self.z0q is None):
            raise ValueError("Tile: give both z0t and z0q, or neither")
        if self.z0t is not None:
            object.__setattr__(self, "scalar_scheme", "prescribed")
        if self.scalar_scheme not in ("a87", "ratio", "prescribed"):
            raise ValueError("Tile.scalar_scheme must be 'a87', 'ratio' or 'prescribed'")
        if self.scalar_scheme == "prescribed" and self.z0t is None:
            raise ValueError("Tile: scalar_scheme='prescribed' needs z0t and z0q")
        if self.kind == "water" and (self.z0 is None) != (self.z0t is None):
            raise ValueError("a water tile with prescribed roughness needs z0, z0t and z0q")
        if self.kind == "water" and (isinstance(self.z0, str) or callable(self.z0)):
            raise ValueError("a water tile's z0 must be a number or an array (m)")

    def with_(self, **changes):
        """A copy with some settings changed."""
        return replace(self, **changes)


def _water(sic):
    return 1.0 - sic


def _ice(sic):
    return sic


def two_tile(water=None, **ice):
    """Open water and ice: fractions 1 - sic and sic.

    Keywords set the ice tile; ``water`` is a dict of settings for the water tile.
    """
    return [
        Tile("water", "water", _water, **(water or {})),
        Tile("thick_ice", "ice", _ice, form_drag=True, **ice),
    ]


def three_tile(thin_fraction, of="water", thin=None, water=None, **thick):
    """Open water, thin ice and thick ice.

    thin_fraction : number or one value per sample, between 0 and 1: the share that is
        thin ice, of the part of the surface named by ``of``.
    of : "water" (default): thin ice (nilas, grey ice) covers that share of the area
        that the ice-concentration product sees as open water. Fractions are then
        water (1 - f)(1 - sic), thin ice f (1 - sic), thick ice sic, and the thin ice
        adds no floe-edge form drag. "ice": thin ice is that share of the ice
        concentration: water 1 - sic, thin ice f sic, thick ice (1 - f) sic, and it
        counts for form drag. Caveat: the form drag (Lupkes et al. 2012 via Elvidge et
        al. 2016) then treats thin ice as floes with the freeboard of the MIZ
        parametrisation (about 0.3 to 0.5 m), far above that of 0.1 m ice.
    thin : dict of settings for the thin-ice tile. Provisional defaults: 0.1 m of bare
        ice (h_ice=0.1, h_snow=0.0). Its surface temperature method is the ``t_ice``
        argument of :func:`coare_seaice.ice_fluxes` unless set here; an energy-balance
        method (``SEB()``) is what makes a thin tile differ from the thick one.
    water : dict of settings for the water tile. Other keywords set the thick-ice tile.

    With thin_fraction = 0 the result equals :func:`two_tile`.
    """
    if of not in ("water", "ice"):
        raise ValueError("three_tile: of must be 'water' or 'ice'")
    f = np.asarray(thin_fraction, dtype=float)
    if np.any(~np.isfinite(f)) or np.any((f < 0) | (f > 1)):
        raise ValueError("three_tile: thin_fraction must be between 0 and 1")
    thin_kw = dict(h_ice=0.1, h_snow=0.0)
    thin_kw.update(thin or {})
    if of == "water":
        fractions = (lambda s: (1 - f) * (1 - s), lambda s: f * (1 - s), _ice)
    else:
        fractions = (_water, lambda s: f * s, lambda s: (1 - f) * s)
    return [
        Tile("water", "water", fractions[0], **(water or {})),
        Tile("thin_ice", "ice", fractions[1], form_drag=(of == "ice"), **thin_kw),
        Tile("thick_ice", "ice", fractions[2], form_drag=True, **thick),
    ]


def tile_fractions(tiles, sic):
    """Area fraction of each tile (dict by name) for ice concentration ``sic``.

    Raises if names repeat or the fractions do not sum to one (to 1e-9) where sic is
    finite. NaN ice concentration gives NaN fractions.
    """
    sic = np.asarray(sic, dtype=float)
    names = [t.name for t in tiles]
    if len(set(names)) != len(names):
        raise ValueError("tile names must be unique")
    out = {t.name: np.broadcast_to(np.asarray(t.fraction(sic), dtype=float), sic.shape).copy()
           for t in tiles}  # fmt: skip
    total = sum(out.values())
    ok = np.isfinite(sic)
    if np.any(np.abs(total[ok] - 1.0) > 1e-9):
        raise ValueError("tile fractions do not sum to one")
    if any(np.any((v[ok] < -1e-12) | (v[ok] > 1 + 1e-12)) for v in out.values()):
        raise ValueError("tile fractions must be between 0 and 1")
    return out
