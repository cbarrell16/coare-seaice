"""Surface temperature methods behind one interface.

A method is any object with a ``solve(ctx)`` method that takes a :class:`TileContext`
and returns a :class:`SurfaceTemperature`. The tile code knows nothing else about it, so
a new method needs no change anywhere else: write a class with ``solve`` and pass an
instance as ``t_ice=`` or on a tile. (In software terms this is a "protocol": the
interface is the method name and its arguments, not a shared parent class.)

Methods here: :class:`AirFreezing` (the simple reference), :class:`AirCapped`,
:class:`Empirical`, :class:`Prescribed`, :class:`SEB` and the experimental
:class:`Column`.
"""

from __future__ import annotations

import warnings
from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np

from ..thermo import LAPSE_DRY, T_FREEZE_SW
from . import simple
from .column import IceColumn
from .properties import K_ICE_FRESH
from .seb import effective_conductance, solve_seb_steady

__all__ = [
    "SEB",
    "AirCapped",
    "AirFreezing",
    "Column",
    "Empirical",
    "Prescribed",
    "SurfaceTemperature",
    "TileContext",
    "method_from_name",
]


@dataclass
class TileContext:
    """What a surface temperature method may use. Arrays have one value per sample.

    ``flux(Ts)`` runs the tile's surface-layer solver for surface temperature Ts (deg C)
    and returns its dict (``hsb``, ``hlb``, ``h_s``, ``h_l``, ...). Thickness, albedo and
    conductivity are the tile's values (or the method's own, if it sets them).
    """

    t_air: np.ndarray  # deg C at zt
    zt: np.ndarray  # m
    p: np.ndarray  # mb
    q_air: np.ndarray  # kg/kg
    lw_down: np.ndarray  # W m-2
    sw_down: np.ndarray  # W m-2
    sic: np.ndarray  # total ice concentration, 0-1
    flux: Callable
    h_ice: np.ndarray  # m
    h_snow: np.ndarray  # m
    k_snow: float  # W m-1 K-1
    albedo: np.ndarray
    t_base: float  # deg C
    time: np.ndarray | None = None  # s
    v2_compat: bool = False

    @property
    def n(self):
        return self.t_air.size

    @property
    def conductance(self):
        """Steady conductance of snow over fresh ice (W m-2 K-1)."""
        return effective_conductance(self.h_ice, self.h_snow, K_ICE_FRESH, self.k_snow)


@dataclass
class SurfaceTemperature:
    """Result of a method: ``t_surface`` (deg C) and optional diagnostics.

    ``filled`` marks samples whose value does not come from the method itself.
    ``seb_residual`` (W m-2) and ``f_cond`` (W m-2, positive upward) are given by methods
    that solve an energy balance; others leave them None.
    """

    t_surface: np.ndarray
    filled: np.ndarray | None = None
    seb_residual: np.ndarray | None = None
    f_cond: np.ndarray | None = None


@dataclass(frozen=True)
class AirFreezing:
    """T_s = min(T_air, t_freeze). The simple reference (v1 behaviour)."""

    t_freeze: float = T_FREEZE_SW
    energy_balance = False

    def solve(self, ctx):
        return SurfaceTemperature(simple.air_freezing(ctx.t_air, self.t_freeze))


@dataclass(frozen=True)
class AirCapped:
    """T_s = min(T_air, cap): the surface follows the air up to the melting point."""

    cap: float = 0.0
    energy_balance = False

    def solve(self, ctx):
        return SurfaceTemperature(simple.air_capped(ctx.t_air, self.cap))


@dataclass(frozen=True)
class Empirical:
    """T_s = min(func(T_air), 0), or T_air - offset if no function is given."""

    offset: float = 2.0
    func: Callable | None = None
    energy_balance = False

    def solve(self, ctx):
        return SurfaceTemperature(simple.empirical(ctx.t_air, self.offset, self.func))


@dataclass(frozen=True)
class SEB:
    """Steady surface energy balance with conduction through snow and ice.

    ``h_ice``, ``h_snow`` (m), ``k_snow``, ``albedo`` and ``t_base`` override the tile's
    values when given. Needs downwelling longwave and shortwave.
    """

    h_ice: float | np.ndarray | None = None
    h_snow: float | np.ndarray | None = None
    k_snow: float | None = None
    albedo: float | np.ndarray | None = None
    t_base: float | None = None
    energy_balance = True

    def solve(self, ctx):
        Ts, res = solve_seb_steady(
            ctx.flux, ctx.lw_down, ctx.sw_down, ctx.albedo, ctx.conductance, ctx.t_base, ctx.t_air
        )
        return SurfaceTemperature(Ts, seb_residual=res)


@dataclass(frozen=True)
class Prescribed:
    """A given surface temperature (deg C), one value per sample or a constant.

    Missing values stay NaN unless ``fallback`` (another method) is given; filled samples
    are flagged.
    """

    values: float | np.ndarray = np.nan
    fallback: object | None = None
    energy_balance = False

    def solve(self, ctx):
        T = np.array(self.values, dtype=float)
        if T.size == 1:
            T = np.full(ctx.n, float(T.reshape(-1)[0]))
        elif T.size != ctx.n:
            raise ValueError("Prescribed: surface temperature array has wrong length")
        filled = np.zeros(ctx.n, dtype=bool)
        miss = ~np.isfinite(T)
        if miss.any() and self.fallback is not None:
            T[miss] = self.fallback.solve(ctx).t_surface[miss]
            filled = miss & np.isfinite(T)
        return SurfaceTemperature(T, filled=filled)


@dataclass(frozen=True)
class Column:
    """Time-dependent snow and ice conduction model (:class:`IceColumn`). Experimental.

    One column is one slab carried along the series: thickness and albedo arrays are
    replaced by their medians, with a warning. Needs ``time`` (s). The turbulent
    conductances are those of the steady solution, so the reported terms do not close
    exactly. Outside the ice (SIC below ``sic_reset``) and across gaps the steady value is
    used and flagged. ``column_kwargs`` go to :class:`IceColumn`.
    """

    h_ice: float | np.ndarray | None = None
    h_snow: float | np.ndarray | None = None
    k_snow: float | None = None
    albedo: float | np.ndarray | None = None
    t_base: float | None = None
    sic_reset: float | None = None
    max_gap_s: float | None = None
    column_kwargs: dict = field(default_factory=dict)
    energy_balance = True

    def solve(self, ctx):
        if ctx.time is None:
            raise ValueError("the Column method needs time (s)")
        h_ice, h_snow, albedo = ctx.h_ice, ctx.h_snow, ctx.albedo
        if not ctx.v2_compat:
            varying = [
                name
                for name, x in (("h_ice", h_ice), ("h_snow", h_snow), ("albedo", albedo))
                if np.nanmax(x) - np.nanmin(x) > 1e-9
            ]
            if varying:
                warnings.warn(
                    f"the column (diffusion) method models one slab: {', '.join(varying)} "
                    "vary along the series and are replaced by their medians",
                    stacklevel=2,
                )
            h_ice = np.full(ctx.n, np.nanmedian(h_ice))
            h_snow = np.full(ctx.n, np.nanmedian(h_snow))
            albedo = np.full(ctx.n, np.nanmedian(albedo))
        G = effective_conductance(h_ice, h_snow, K_ICE_FRESH, ctx.k_snow)
        T_ss, _ = solve_seb_steady(
            ctx.flux, ctx.lw_down, ctx.sw_down, albedo, G, ctx.t_base, ctx.t_air
        )
        fl_ss = ctx.flux(T_ss)  # conductances at the steady state
        run_kw = {
            k: v for k, v in (("sic_reset", self.sic_reset), ("max_gap_s", self.max_gap_s))
            if v is not None
        }  # fmt: skip
        col = IceColumn(
            h_ice=float(np.nanmedian(h_ice)),
            h_snow=float(np.nanmedian(h_snow)),
            t_base=ctx.t_base,
            k_snow=ctx.k_snow,
            v2_compat=ctx.v2_compat,
            **self.column_kwargs,
        )
        res = col.run(
            np.asarray(ctx.time, dtype=float), ctx.t_air + LAPSE_DRY * ctx.zt, ctx.q_air, ctx.p,
            ctx.lw_down, ctx.sw_down, fl_ss["h_s"], fl_ss["h_l"], T_ss, sic=ctx.sic,
            albedo=float(np.nanmedian(albedo)), **run_kw,
        )  # fmt: skip
        ran = np.isfinite(res["Ts"])
        Ts = np.where(ran, res["Ts"], T_ss)
        if ctx.v2_compat:
            # v2 back-calculated the flux through an effective conductance, which is wrong
            # (and can be enormous) when the surface is warmer than the base.
            G = np.where(
                np.isfinite(res["F_cond"]),
                res["F_cond"] / np.maximum(ctx.t_base - np.minimum(Ts, 0.0), 1e-6),
                G,
            )
        f_cond = G * (ctx.t_base - np.minimum(Ts, 0.0))  # steady value where not running
        if not ctx.v2_compat:
            f_cond = np.where(np.isfinite(res["F_cond"]), res["F_cond"], f_cond)
        return SurfaceTemperature(
            Ts, filled=~ran & np.isfinite(T_ss), seb_residual=res["melt"], f_cond=f_cond
        )


_NAMES = {
    "air_freezing": AirFreezing,
    "air_capped": AirCapped,
    "empirical": Empirical,
    "seb_steady": SEB,
    "seb": SEB,
    "diffusion": Column,
    "column": Column,
}


def method_from_name(name, **kwargs):
    """A method from its name: 'air_freezing', 'air_capped', 'empirical', 'seb', 'column'
    (the v2 names 'seb_steady' and 'diffusion' also work)."""
    if name not in _NAMES:
        raise ValueError(f"unknown surface temperature method {name!r}; use {sorted(_NAMES)}")
    return _NAMES[name](**kwargs)
