"""Task P5: reprocess the ROVER 2025 ship fluxes with coare-seaice and check them against v2.

Runs on the UEA machine. It does two things:

1. **Check.** Runs `coare_seaice.legacy.coare35vn_ice` and the original v2
   `coare35vn_ice.py` on identical inputs and reports, for every one of the 56 output
   columns, the largest difference and whether NaNs fall in the same places.
2. **Reprocess.** Writes the `air_freezing` fluxes (the v2 bug fixes only, no new
   ice-surface-temperature physics) to a new netCDF file with the run settings recorded
   in its attributes.

Inputs are built the same way as in `run_tice_comparison.py` (v2): ship wind, air
temperature, relative humidity, pressure, SST, latitude and AMSR2 sea-ice concentration;
shortwave from the ship where available, else ERA5; longwave from ERA5 along-track.
Missing values are left as NaN (no fills), so a window with any missing input gives NaN.

Paths are read from the config block below. By default everything is found under the
directory named by the environment variable COARE_SEAICE_DATA.

Log: <OUT_DIR>/logs/p5_reprocess_rover_legacy.log
"""

from __future__ import annotations

import contextlib
import io
import logging
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
DATA = Path(os.environ.get("COARE_SEAICE_DATA", "."))
SHIP_NC = DATA / "ROVER2025_KronprinsHaakon_surface_met_fluxes_10min.nc"
ERA5_TRACK_NC = DATA / "era5_along_track_10min.nc"  # optional; set to None to skip
V2_DIR = DATA / "coare_ice_v2"  # folder holding coare35vn_ice.py v2.0.0; None to skip the check
OUT_DIR = Path(os.environ.get("COARE_SEAICE_OUT", "outputs"))

ZU, ZT, ZQ = 36.0, 22.0, 22.0  # wind, temperature, humidity heights (m)
Z0_ICE = 5e-4  # ice roughness length (m), v2 default
JCOOL = 1  # COARE cool skin on
REPROCESS = True  # False: keep an existing output file
# ---------------------------------------------------------------------------

warnings.filterwarnings("ignore", category=RuntimeWarning)
NAME = "p5_reprocess_rover_legacy"
(OUT_DIR / "logs").mkdir(parents=True, exist_ok=True)
(OUT_DIR / "data").mkdir(parents=True, exist_ok=True)
log = logging.getLogger(NAME)
log.setLevel(logging.INFO)
log.handlers.clear()
for h in (
    logging.FileHandler(OUT_DIR / "logs" / f"{NAME}.log", mode="w"),
    logging.StreamHandler(sys.stdout),
):
    h.setFormatter(logging.Formatter("%(message)s"))
    log.addHandler(h)

import coare_seaice  # noqa: E402
from coare_seaice import legacy  # noqa: E402

OUT_NC = OUT_DIR / "data" / "ROVER2025_coare_seaice_air_freezing_10min.nc"
log.info(f"{NAME} | coare-seaice {coare_seaice.__version__}")
log.info(f"Ship file:  {SHIP_NC}\nERA5 file:  {ERA5_TRACK_NC}\nv2 folder:  {V2_DIR}")
log.info(f"Output:     {OUT_NC}\n")

# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
df = xr.open_dataset(SHIP_NC).to_dataframe()
u = df.wind_speed.values
t = df.air_temperature.values
rh = df.relative_humidity.values
P = df.air_pressure.values
ts = df.sea_surface_temperature.values
lat = df.latitude.values
sic = df.sea_ice_concentration.values

Rl = np.full(len(df), np.nan)
Rs_era5 = np.full(len(df), np.nan)
lw_source = "none (open-water cool skin and ice radiation columns will be NaN)"
if ERA5_TRACK_NC is not None and Path(ERA5_TRACK_NC).exists():
    e5 = xr.open_dataset(ERA5_TRACK_NC).to_dataframe().reindex(df.index)
    for lw_name, sw_name in (("era5_lw_down", "era5_sw_down"), ("strd", "ssrd")):
        if lw_name in e5:
            Rl, Rs_era5 = e5[lw_name].values, e5[sw_name].values
            lw_source = f"ERA5 along-track ({lw_name})"
            break
Rs = np.where(np.isfinite(df.shortwave_down.values), df.shortwave_down.values, Rs_era5)

valid = np.isfinite(u) & np.isfinite(t) & np.isfinite(rh) & np.isfinite(P)
valid &= np.isfinite(ts) & np.isfinite(lat)
vi = np.where(valid)[0]
log.info(f"10-min windows: {len(df)}; with valid met inputs: {valid.sum()}")
log.info(f"Longwave: {lw_source}; missing in {np.isnan(Rl[vi]).sum()} valid windows")
log.info(f"Sea-ice concentration missing in {np.isnan(sic[vi]).sum()} valid windows (left NaN)")
log.info("Note: ERA5 downwelling longwave is probably too high here (ERA5 overestimates cloud).\n")


def sub(x):
    return np.asarray(x)[vi]


kw = dict(
    P=sub(P), Rs=sub(Rs), Rl=sub(Rl), zu=ZU, zt=ZT, zq=ZQ, lat=sub(lat), jcool=JCOOL,
    ice_conc=sub(sic), z0_ice=Z0_ICE, t_ice_method="air_freezing", verbose=False,
)  # fmt: skip
# t_ice_method is set explicitly: the v2 default is 'seb_steady', not 'air_freezing'.


def run(fn):
    with contextlib.redirect_stdout(io.StringIO()):
        a = fn(sub(u), sub(t), sub(rh).copy(), sub(ts), **kw)
    full = np.full((len(df), a.shape[1]), np.nan)
    full[vi] = a
    return full


new = run(legacy.coare35vn_ice)
names = list(legacy.OUTPUT_COLUMNS)

# ---------------------------------------------------------------------------
# 1. Check against v2
# ---------------------------------------------------------------------------
v2 = None
if V2_DIR is not None and (Path(V2_DIR) / "coare35vn_ice.py").exists():
    import inspect

    sys.path.insert(0, str(V2_DIR))
    import coare35vn_ice as v2  # noqa: E402

    v2_version = getattr(v2, "__version__", "unknown")
    if "t_ice_method" not in inspect.signature(v2.coare35vn_ice).parameters:
        log.info(
            f"The coare35vn_ice.py in {V2_DIR} (version {v2_version}) has no 't_ice_method' "
            "argument, so it is v1, not v2.0.0. Check against v2 skipped.\n"
            "To run the check, put the v2.0.0 files (coare35vn_ice.py, "
            "ice_surface_temperature.py, meteo.py, util.py) in that folder.\n"
        )
        v2 = None
elif V2_DIR is not None:
    log.info(f"No coare35vn_ice.py in {V2_DIR}: check against v2 skipped.\n")

if v2 is not None:
    old = run(v2.coare35vn_ice)
    log.info(f"Check against v2 {getattr(v2, '__version__', '?')} (air_freezing):")
    log.info(f"{'col':>3s} {'name':12s} {'max |diff|':>11s} {'max rel':>9s} {'NaN mismatch':>12s}")
    ice_only = set(range(38, 47)) | set(range(50, 55))
    worst_ice = 0.0
    for j, nm in enumerate(names):
        d = np.abs(new[:, j] - old[:, j])
        rel = d / np.maximum(np.abs(old[:, j]), 1e-12)
        nan_mis = int((np.isnan(new[:, j]) != np.isnan(old[:, j])).sum())
        dmax = np.nanmax(d) if np.isfinite(d).any() else 0.0
        rmax = np.nanmax(rel) if np.isfinite(rel).any() else 0.0
        if j in ice_only:
            worst_ice = max(worst_ice, dmax)
        log.info(f"{j:3d} {nm:12s} {dmax:11.3g} {rmax:9.2g} {nan_mis:12d}")
    log.info(
        "\nExpected: ice-only columns (38-46, 50-54) identical to v2; open-water columns "
        "differ by up to about 0.1 W m-2 (heat) and 1e-4 N m-2 (stress) because pycoare "
        "and the NOAA code differ slightly in gravity and a few constants; no NaN mismatches."
    )
    log.info(f"Largest difference in an ice-only column: {worst_ice:.3g}\n")
elif V2_DIR is None:
    log.info("V2_DIR is None: check against v2 skipped.\n")

# ---------------------------------------------------------------------------
# 2. Write the reprocessed fluxes
# ---------------------------------------------------------------------------
if OUT_NC.exists() and not REPROCESS:
    log.info(f"{OUT_NC} exists and REPROCESS is False: not overwritten.")
else:
    keep = {
        "tau_blend": ("N m-2", "area-mean momentum flux including form drag"),
        "hsb_blend": ("W m-2", "area-mean sensible heat flux, positive upward"),
        "hlb_blend": ("W m-2", "area-mean latent heat flux, positive upward"),
        "tau_form": ("N m-2", "form drag"),
        "tau": ("N m-2", "open-water momentum flux"),
        "hsb": ("W m-2", "open-water sensible heat flux, positive upward"),
        "hlb": ("W m-2", "open-water latent heat flux, positive upward"),
        "tau_ice": ("N m-2", "ice-tile momentum flux (skin drag)"),
        "hsb_ice": ("W m-2", "ice-tile sensible heat flux, positive upward"),
        "hlb_ice": ("W m-2", "ice-tile latent heat flux, positive upward"),
        "T_ice": ("degC", "ice surface temperature, min(T_air, -1.8 C)"),
        "U10N": ("m s-1", "10 m neutral wind over open water"),
    }
    out = xr.Dataset(coords={"time": df.index.values})
    for nm, (units, desc) in keep.items():
        out[nm] = ("time", new[:, names.index(nm)], {"units": units, "long_name": desc})
    out.attrs = {
        "title": "ROVER 2025 Kronprins Haakon turbulent fluxes, reprocessed with coare-seaice",
        "source": f"coare-seaice {coare_seaice.__version__}, coare_seaice.legacy.coare35vn_ice",
        "t_ice_method": "air_freezing (T_ice = min(T_air, -1.8 C))",
        "z0_ice_m": Z0_ICE,
        "heights_m": f"zu={ZU}, zt={ZT}, zq={ZQ}",
        "ice_scheme": "blended_a87 (Andreas 1987 scalar roughness over ice)",
        "form_drag": "Lupkes et al. (2012) with Elvidge et al. (2016) E2016B parameters",
        "longwave_source": lw_source,
        "sea_ice_concentration": "AMSR2, from the ship file; missing values left NaN",
        "caveats": (
            "ERA5 downwelling longwave is probably too high (ERA5 overestimates cloud). "
            "The non-ice fraction is treated as open water at the SST; across the MIZ this "
            "probably overestimates heat fluxes. Ice roughness is regional and uncertain."
        ),
        "restricted": "ROVER working group only",
    }
    out.to_netcdf(OUT_NC)
    log.info(f"Wrote {OUT_NC} ({len(keep)} variables).")
