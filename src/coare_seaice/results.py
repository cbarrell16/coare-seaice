"""Named results: per tile, blended, quality flags, and conversion to xarray.

Names follow ``pycoare`` where it has one (``tau``, ``hsb``, ``hlb``, ``usr``). Heat fluxes
are positive upward. Every array has one value per input sample.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

__all__ = ["BlendedResult", "FluxResult", "TileResult"]

# name -> (units, description)
TILE_FIELDS = {
    "fraction": ("1", "area fraction of the tile"),
    "tau": ("N m-2", "wind stress"),
    "hsb": ("W m-2", "sensible heat flux, positive upward"),
    "hlb": ("W m-2", "latent heat flux, positive upward"),
    "usr": ("m s-1", "friction velocity"),
    "t_surface": ("degC", "surface temperature used for the fluxes"),
    "q_surface": ("kg kg-1", "surface specific humidity"),
    "z0": ("m", "roughness length for momentum"),
    "z0t": ("m", "roughness length for temperature"),
    "z0q": ("m", "roughness length for humidity"),
    "cdn10": ("1", "neutral 10 m drag coefficient (skin part)"),
    "chn10": ("1", "neutral 10 m heat transfer coefficient"),
    "cen10": ("1", "neutral 10 m moisture transfer coefficient"),
    "rstar": ("1", "roughness Reynolds number"),
    "L": ("m", "Obukhov length"),
    "zeta": ("1", "stability parameter zu / L"),
    "rib": ("1", "bulk Richardson number"),
    "wind_eff": ("m s-1", "effective wind speed used (floor or light-wind terms)"),
    "lw_net": ("W m-2", "net longwave at the surface, positive downward"),
    "sw_net": ("W m-2", "net shortwave at the surface, positive downward"),
    "f_cond": ("W m-2", "conductive flux to the surface, positive upward"),
    "seb_residual": ("W m-2", "surface energy balance residual (melt flux if melting)"),
}
BLENDED_FIELDS = {
    "tau": ("N m-2", "area-mean wind stress, including form drag"),
    "hsb": ("W m-2", "area-mean sensible heat flux, positive upward"),
    "hlb": ("W m-2", "area-mean latent heat flux, positive upward"),
    "tau_form": ("N m-2", "form-drag part of the stress"),
    "t_surface": ("degC", "area-mean (or the single aggregate) surface temperature"),
    "cdn10": ("1", "blended neutral 10 m drag coefficient, including form drag"),
    "chn10": ("1", "blended neutral 10 m heat transfer coefficient"),
    "cen10": ("1", "blended neutral 10 m moisture transfer coefficient"),
    "cdn10_form": ("1", "form-drag part of the neutral 10 m drag coefficient"),
}
QC_FLAGS = {
    "missing_input": "an input needed for the fluxes is NaN (outputs are NaN)",
    "rib_above_critical": "bulk Richardson number above 0.2 over a tile: outside the "
    "validated range of the stability functions (Grachev et al. 2007)",
    "zeta_limited": "the stability limit acted in the iteration: the result is set by it",
    "height_above_obukhov_length": "stable and zu > L: the measurement is probably above "
    "the surface layer",
    "light_wind": "wind speed below 1 m s-1",
    "rstar_out_of_range": "roughness Reynolds number above 1000: beyond the Andreas (1987) "
    "fit (ratio held at its R* = 1000 value)",
    "not_converged": "the surface-layer iteration was still changing at its last step",
    "t_surface_filled": "a tile's surface temperature does not come from the requested "
    "method (fallback or steady-state value)",
    "melting_surface": "an ice tile is at 0 C",
    "humidity_far_from_ice_saturation": "air humidity far from ice saturation over mostly "
    "ice-covered surface (Andreas et al. 2002): check the humidity measurement",
}

# The 18 ice columns of the legacy array (38-55): (legacy name, source, field, scale)
_LEGACY_ICE = [
    ("tau_ice", "ice", "tau", 1.0), ("hsb_ice", "ice", "hsb", 1.0), ("hlb_ice", "ice", "hlb", 1.0),
    ("Cdn_ice", "ice", "cdn10", 1e3), ("Chn_ice", "ice", "chn10", 1e3),
    ("Cen_ice", "ice", "cen10", 1e3), ("R_star", "ice", "rstar", 1.0),
    ("z0T_ice", "ice", "z0t", 1.0), ("z0q_ice", "ice", "z0q", 1.0),
    ("tau_blend", "blended", "tau", 1.0), ("hsb_blend", "blended", "hsb", 1.0),
    ("hlb_blend", "blended", "hlb", 1.0), ("T_ice", "ice", "t_surface", 1.0),
    ("F_cond", "ice", "f_cond", 1.0), ("LWnet_ice", "ice", "lw_net", 1.0),
    ("SWnet_ice", "ice", "sw_net", 1.0), ("SEB_res", "ice", "seb_residual", 1.0),
    ("tau_form", "blended", "tau_form", 1.0),
]  # fmt: skip


@dataclass
class TileResult:
    """Results for one surface tile. Fields that do not apply to a tile are None."""

    name: str
    kind: str  # "water" or "ice"
    fraction: np.ndarray
    tau: np.ndarray | None = None
    hsb: np.ndarray | None = None
    hlb: np.ndarray | None = None
    usr: np.ndarray | None = None
    t_surface: np.ndarray | None = None
    q_surface: np.ndarray | None = None
    z0: np.ndarray | None = None
    z0t: np.ndarray | None = None
    z0q: np.ndarray | None = None
    cdn10: np.ndarray | None = None
    chn10: np.ndarray | None = None
    cen10: np.ndarray | None = None
    rstar: np.ndarray | None = None
    L: np.ndarray | None = None
    zeta: np.ndarray | None = None
    rib: np.ndarray | None = None
    wind_eff: np.ndarray | None = None
    lw_net: np.ndarray | None = None
    sw_net: np.ndarray | None = None
    f_cond: np.ndarray | None = None
    seb_residual: np.ndarray | None = None
    qc: dict = field(default_factory=dict)  # flag name -> boolean array
    extra: dict = field(default_factory=dict, repr=False)  # not exported (e.g. pycoare object)


@dataclass
class BlendedResult:
    """Area-mean (blended) results."""

    tau: np.ndarray
    hsb: np.ndarray
    hlb: np.ndarray
    tau_form: np.ndarray
    t_surface: np.ndarray | None = None
    cdn10: np.ndarray | None = None
    chn10: np.ndarray | None = None
    cen10: np.ndarray | None = None
    cdn10_form: np.ndarray | None = None


@dataclass
class FluxResult:
    """Everything :func:`coare_seaice.ice_fluxes` returns.

    ``blended`` are the area means, ``tiles`` the per-tile results by name, ``qc`` the
    quality flags (boolean arrays, True = flagged; see ``QC_FLAGS`` for their meaning) and
    ``meta`` the settings used. Flags never change a value.
    """

    blended: BlendedResult
    tiles: dict
    qc: dict = field(default_factory=dict)
    meta: dict = field(default_factory=dict)

    @property
    def n(self):
        return self.blended.tau.size

    @property
    def water(self):
        """The (first) open-water tile."""
        return self._first("water")

    @property
    def ice(self):
        """The first ice tile (the thick ice in the standard tile sets)."""
        return self._first("ice")

    def _first(self, kind):
        for tile in self.tiles.values():
            if tile.kind == kind:
                return tile
        raise KeyError(f"no {kind} tile")

    # -- xarray ---------------------------------------------------------------------------
    def to_xarray(self, time=None):
        """All results as an ``xarray.Dataset`` along the dimension ``sample``.

        Variables are named ``blended_<field>``, ``<tile>_<field>`` and ``qc_<flag>``
        (per-tile flags ``<tile>_qc_<flag>``), with units and descriptions as attributes.
        ``time`` (optional, one value per sample) becomes a coordinate. Needs xarray,
        which is installed with the ``verify`` extra.
        """
        import xarray as xr

        def var(values, units, long_name):
            return ("sample", np.asarray(values), {"units": units, "long_name": long_name})

        data = {}
        for name, (units, desc) in BLENDED_FIELDS.items():
            values = getattr(self.blended, name)
            if values is not None:
                data[f"blended_{name}"] = var(values, units, desc)
        for tname, tile in self.tiles.items():
            for name, (units, desc) in TILE_FIELDS.items():
                values = getattr(tile, name)
                if values is not None:
                    data[f"{tname}_{name}"] = var(values, units, f"{desc} ({tname} tile)")
            for flag, values in tile.qc.items():
                data[f"{tname}_qc_{flag}"] = var(values, "1", QC_FLAGS.get(flag, flag))
        for flag, values in self.qc.items():
            data[f"qc_{flag}"] = var(values, "1", QC_FLAGS.get(flag, flag))
        ds = xr.Dataset(data)
        if time is not None:
            ds = ds.assign_coords(time=("sample", np.asarray(time)))
        ds.attrs.update({k: str(v) for k, v in self.meta.items()})
        ds.attrs["sign_convention"] = "heat fluxes positive upward; radiation positive downward"
        return ds

    # -- legacy array -----------------------------------------------------------------------
    def to_legacy_array(self):
        """The 56-column array of ``legacy.coare35vn_ice`` (two tiles only).

        Needs the 38 open-water columns in ``self.water.extra["legacy_columns"]``; the
        legacy wrapper and :meth:`from_legacy_array` put them there.
        """
        if len(self.tiles) != 2 or "legacy_columns" not in self.water.extra:
            raise ValueError("to_legacy_array needs a two-tile result made by the legacy wrapper")
        src = {"ice": self.ice, "blended": self.blended}
        ice_cols = [getattr(src[s], f) * scale for _, s, f, scale in _LEGACY_ICE]
        return np.column_stack([self.water.extra["legacy_columns"], *ice_cols])

    @classmethod
    def from_legacy_array(cls, array, sic, diagnostics=None):
        """Named results from a 56-column legacy array and the ice concentration used."""
        from .legacy import OUTPUT_COLUMNS

        A = np.asarray(array, dtype=float)
        if A.ndim != 2 or A.shape[1] != len(OUTPUT_COLUMNS):
            raise ValueError("expected an (N, 56) legacy array")
        col = {name: A[:, i] for i, name in enumerate(OUTPUT_COLUMNS)}
        sic = np.broadcast_to(np.asarray(sic, dtype=float), (A.shape[0],)).copy()
        water = TileResult(
            name="water", kind="water", fraction=1.0 - sic, tau=col["tau"], hsb=col["hsb"],
            hlb=col["hlb"], usr=col["usr"], z0t=col["zot"], z0q=col["zoq"], L=col["L"],
            zeta=col["zet"], cdn10=col["Cdn_10"] / 1e3, chn10=col["Chn_10"] / 1e3,
            cen10=col["Cen_10"] / 1e3, q_surface=col["Qs"] / 1e3,
            extra={"legacy_columns": A[:, :38].copy()},
        )  # fmt: skip
        ice = TileResult(name="thick_ice", kind="ice", fraction=sic,
                         qc=dict(diagnostics or {}))  # fmt: skip
        blended = BlendedResult(tau=col["tau_blend"], hsb=col["hsb_blend"],
                                hlb=col["hlb_blend"], tau_form=col["tau_form"])  # fmt: skip
        for name, s, f, scale in _LEGACY_ICE:
            if s == "ice":
                setattr(ice, f, col[name] / scale)
        return cls(blended=blended, tiles={"water": water, "thick_ice": ice},
                   meta={"source": "legacy array"})  # fmt: skip
