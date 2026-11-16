"""ROVER 2025 (00b): recompute the ship fluxes with coare-seaice and compare with the old ones.

Runs on the UEA machine after ROVERview notebook 00. It reads the 00 output file (ship
met, AMSR2 sea-ice concentration and the old fluxes), recomputes the fluxes with
`coare_seaice.ice_fluxes`, and writes them to a **separate** file, so that the 00 output
is never overwritten. It then compares new with old.

Configurations:
- `ref` (always): the simple reference, matching the old processing: two tiles (open
  water and ice), ice surface at min(T_air, -1.8 C), z0_ice = 5e-4 m, COARE 3.5 stability
  functions, blended-A87 scalar roughness over ice, Lupkes/E2016B form drag. This is the
  like-for-like comparison with the old fluxes.
- `seb` (optional, RUN_SEB): ice surface temperature from the steady surface energy
  balance (1 m ice, 0.2 m snow), using ERA5 downwelling longwave. This is new physics, not
  a check. ERA5 longwave is probably too high here (ERA5 overestimates cloud), so these
  are the warm end of a range.

Inputs are built as in P5: ship wind, air temperature, relative humidity, pressure, SST,
latitude and AMSR2 SIC from the 00 file; shortwave from the ship where available, else
ERA5; longwave from ERA5 along-track (optional file). Nothing is filled: a window with a
missing input gives NaN, and is flagged.

Differences to expect between new and old (new minus old):
- Open water: pycoare against the NOAA Python COARE 3.5: up to about 0.1 W m-2 and
  1e-4 N m-2.
- Ice tile heat fluxes: about -2 % (up to -4 %), from the temperature-dependent air
  viscosity (release 0.3.0, F2).
- If the old fluxes come from v1 (`coare35vn_ice.py` before the November fixes), also:
  blended stress up to about 44 % higher in the MIZ (v1 form-drag bug), and blended
  latent heat 1.6 to 2.8 W m-2 higher (v1 humidity units). These were the reason for P5.
- Windows with missing SIC: NaN blended fluxes (the old code may have treated them as
  open water).
Anything else is unexpected: please send the log.

Log: <OUT_DIR>/logs/rover_00b_fluxes_coare_seaice.log
"""

from __future__ import annotations

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
SHIP_NC = DATA / "ROVER2025_KronprinsHaakon_surface_met_fluxes_10min.nc"  # the 00 output
ERA5_TRACK_NC = DATA / "era5_along_track_10min.nc"  # optional; None to skip
# Old fluxes: by default taken from SHIP_NC. To compare with another file (for example
# the P5 output, outputs/data/ROVER2025_coare_seaice_air_freezing_10min.nc), set OLD_NC.
OLD_NC = None
# Names of the old flux variables. None: found automatically from names and long_names
# (the choice is logged). Keys: {tau,hsb,hlb}_{blend,ow,ice}.
OLD_NAMES: dict | None = None
OUT_DIR = Path(os.environ.get("COARE_SEAICE_OUT", "outputs"))

ZU, ZT, ZQ = 36.0, 22.0, 22.0  # wind, temperature, humidity heights (m)
Z0_ICE = 5e-4  # ice roughness length (m)
JCOOL = 1  # COARE cool skin on
RUN_SEB = True  # also run the energy-balance ice surface temperature
SIC_CLASSES = [(0.0, 0.0), (0.0, 0.15), (0.15, 0.5), (0.5, 0.8), (0.8, 0.95), (0.95, 1.0)]
REPROCESS = True
# ---------------------------------------------------------------------------

warnings.filterwarnings("ignore", category=RuntimeWarning)
NAME = "rover_00b_fluxes_coare_seaice"
for sub in ("logs", "data", "figures"):
    (OUT_DIR / sub).mkdir(parents=True, exist_ok=True)
log = logging.getLogger(NAME)
log.setLevel(logging.INFO)
log.handlers.clear()
for h in (logging.FileHandler(OUT_DIR / "logs" / f"{NAME}.log", mode="w"),
          logging.StreamHandler(sys.stdout)):  # fmt: skip
    h.setFormatter(logging.Formatter("%(message)s"))
    log.addHandler(h)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import coare_seaice as cs  # noqa: E402

plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white", "font.size": 9})
SIC_CMAP = "Blues_r"  # dark blue at 0, white at 1
OUT_NC = OUT_DIR / "data" / f"ROVER2025_fluxes_coare_seaice_{cs.__version__}_10min.nc"
log.info(f"{NAME} | coare-seaice {cs.__version__}")
log.info(f"00 file:   {SHIP_NC}\nERA5 file: {ERA5_TRACK_NC}\nOld fluxes from: {OLD_NC or SHIP_NC}")
log.info(f"Output:    {OUT_NC}\n")

# ---------------------------------------------------------------------------
# Inputs (as in P5)
# ---------------------------------------------------------------------------
ds = xr.open_dataset(SHIP_NC)
df = ds.to_dataframe()
n = len(df)
u = df.wind_speed.to_numpy()
t = df.air_temperature.to_numpy()
rh = df.relative_humidity.to_numpy()
P = df.air_pressure.to_numpy()
sst = df.sea_surface_temperature.to_numpy()
lat = df.latitude.to_numpy()
sic = df.sea_ice_concentration.to_numpy()
if np.nanmax(sic) > 1.5:
    raise ValueError("sea_ice_concentration looks like per cent; expected a fraction 0-1")

Rl = np.full(n, np.nan)
Rs_era5 = np.full(n, np.nan)
lw_source = "none (cool skin and energy-balance outputs will be NaN)"
if ERA5_TRACK_NC is not None and Path(ERA5_TRACK_NC).exists():
    e5 = xr.open_dataset(ERA5_TRACK_NC).to_dataframe().reindex(df.index)
    for lw_name, sw_name in (("era5_lw_down", "era5_sw_down"), ("strd", "ssrd")):
        if lw_name in e5:
            Rl, Rs_era5 = e5[lw_name].to_numpy(), e5[sw_name].to_numpy()
            lw_source = f"ERA5 along-track ({lw_name})"
            break
sw_obs = df.shortwave_down.to_numpy() if "shortwave_down" in df else np.full(n, np.nan)
Rs = np.where(np.isfinite(sw_obs), sw_obs, Rs_era5)

met_ok = np.all([np.isfinite(x) for x in (u, t, rh, P, sst, lat)], axis=0)
log.info(f"10-min windows: {n}; with all met inputs: {met_ok.sum()}; of those with SIC: "
         f"{(met_ok & np.isfinite(sic)).sum()}")  # fmt: skip
log.info(f"Longwave: {lw_source}; missing in {(met_ok & ~np.isfinite(Rl)).sum()} met windows")
log.info("Caveat: ERA5 downwelling longwave is probably too high here (ERA5 overestimates "
         "cloud). It affects the cool skin slightly and the 'seb' run strongly.\n")  # fmt: skip

common = dict(zu=ZU, zt=ZT, zq=ZQ, lat=lat, lw_down=Rl, sw_down=Rs, jcool=JCOOL, z0_ice=Z0_ICE)
runs = {"ref": cs.ice_fluxes(u, t, rh, P, sst, sic, **common)}
if RUN_SEB:
    runs["seb"] = cs.ice_fluxes(u, t, rh, P, sst, sic, t_ice=cs.SEB(), **common)

# ---------------------------------------------------------------------------
# Old fluxes
# ---------------------------------------------------------------------------
old_ds = xr.open_dataset(OLD_NC) if OLD_NC else ds
KIND = {"tau": ("momentum", "stress", "tau"), "hsb": ("sensible", "shf", "hsb"),
        "hlb": ("latent", "lhf", "hlb")}  # fmt: skip
SURF = {"blend": ("blend",), "ow": ("open water", "open-water", "_ow", "openwater"),
        "ice": ("ice fraction", "ice tile", "_ice", "over ice")}  # fmt: skip


def find_old(kind, surf):
    """The old variable for e.g. sensible heat, blended: by name and long_name keywords."""
    hits = []
    for name, var in old_ds.data_vars.items():
        text = f"{name} {var.attrs.get('long_name', '')}".lower()
        if any(x in text for x in ("turbulent", "thf", "coefficient", "scaling")) or (
            "form" in name.lower()
        ):
            continue
        if any(k in text for k in KIND[kind]) and any(s in text for s in SURF[surf]):
            if surf == "blend" or "blend" not in text:
                hits.append(name)
    if surf == "ow":  # an unlabelled flux is open-water COARE (e.g. 'tau', 'hsb')
        hits += [nm for nm in (kind,) if nm in old_ds and nm not in hits]
    return hits


old = {}
for kind in KIND:
    for surf in SURF:
        key = f"{kind}_{surf}"
        if OLD_NAMES and key in OLD_NAMES:
            name = OLD_NAMES[key]
        else:
            hits = find_old(kind, surf)
            if len(hits) != 1:
                log.info(f"  old {key}: {'no match' if not hits else f'ambiguous {hits}'} "
                         "(set OLD_NAMES to choose)")  # fmt: skip
                continue
            name = hits[0]
        old[key] = old_ds[name].to_numpy().astype(float)
        log.info(f"  old {key:10s} <- {name}  ({old_ds[name].attrs.get('long_name', '')})")
if not old:
    raise SystemExit("No old flux variables found: set OLD_NAMES.")

ref = runs["ref"]
new = {
    "tau_blend": ref.blended.tau, "hsb_blend": ref.blended.hsb, "hlb_blend": ref.blended.hlb,
    "tau_ow": ref.water.tau, "hsb_ow": ref.water.hsb, "hlb_ow": ref.water.hlb,
    "tau_ice": ref.ice.tau, "hsb_ice": ref.ice.hsb, "hlb_ice": ref.ice.hlb,
}  # fmt: skip


def stats(d):
    d = d[np.isfinite(d)]
    if not d.size:
        return "n=0"
    return (f"n={d.size:5d}  mean {d.mean():+9.4g}  median {np.median(d):+9.4g}  "
            f"5-95% [{np.percentile(d, 5):+.3g}, {np.percentile(d, 95):+.3g}]  "
            f"max|d| {np.abs(d).max():.3g}")  # fmt: skip


log.info("\n" + "=" * 78 + "\nNew (ref) minus old\n")
for key in new:
    if key not in old:
        continue
    d = new[key] - old[key]
    nan_new_only = int((np.isnan(new[key]) & np.isfinite(old[key])).sum())
    nan_old_only = int((np.isfinite(new[key]) & np.isnan(old[key])).sum())
    log.info(f"{key:10s} {stats(d)}   NaN only in new: {nan_new_only}, only in old: "
             f"{nan_old_only}")  # fmt: skip
    if key.endswith("blend") or key.endswith("ice"):
        for lo, hi in SIC_CLASSES:
            m = (sic == 0) if hi == 0 else (sic > lo) & (sic <= hi)
            rel = d[m] / np.where(np.abs(old[key][m]) > 1e-9, old[key][m], np.nan)
            rel = rel[np.isfinite(rel)]
            rtxt = f"  median rel {np.median(rel):+.2%}" if rel.size else ""
            label = "SIC = 0" if hi == 0 else f"SIC {lo:.2f}-{hi:.2f}"
            log.info(f"    {label:16s} {stats(d[m])}{rtxt}")
log.info(
    "\nRead against the expected differences in this script's docstring (pycoare vs NOAA "
    "COARE ~0.1 W m-2; ice heat fluxes ~-2 % from the viscosity fix; and, if the old "
    "fluxes are from v1, the form-drag and humidity-unit fixes).\n"
)
if "seb" in runs:
    s = runs["seb"]
    log.info("seb run (ERA5 longwave; warm end of a range), median by SIC class:")
    for lo, hi in SIC_CLASSES[1:]:
        m = (sic > lo) & (sic <= hi)
        dt_ice = np.nanmedian(s.ice.t_surface[m] - t[m])
        log.info(f"    SIC {lo:.2f}-{hi:.2f}: T_ice - T_air {dt_ice:+.2f} K,"
                 f" blended H {np.nanmedian(s.blended.hsb[m]):+.1f} (ref "
                 f"{np.nanmedian(ref.blended.hsb[m]):+.1f}) W m-2")  # fmt: skip

# ---------------------------------------------------------------------------
# Output (separate file)
# ---------------------------------------------------------------------------
if OUT_NC.exists() and not REPROCESS:
    log.info(f"\n{OUT_NC} exists and REPROCESS is False: not overwritten.")
else:
    VARS = {
        "tau_blend": ("N m-2", "area-mean momentum flux including form drag"),
        "hsb_blend": ("W m-2", "area-mean sensible heat flux, positive upward"),
        "hlb_blend": ("W m-2", "area-mean latent heat flux, positive upward"),
        "tau_ow": ("N m-2", "open-water momentum flux (COARE 3.5)"),
        "hsb_ow": ("W m-2", "open-water sensible heat flux, positive upward"),
        "hlb_ow": ("W m-2", "open-water latent heat flux, positive upward"),
        "tau_ice": ("N m-2", "ice-tile momentum flux (skin drag)"),
        "hsb_ice": ("W m-2", "ice-tile sensible heat flux, positive upward"),
        "hlb_ice": ("W m-2", "ice-tile latent heat flux, positive upward"),
    }
    out = xr.Dataset(coords={"time": df.index.to_numpy()})
    for nm, (units, desc) in VARS.items():
        out[nm] = ("time", new[nm], {"units": units, "long_name": desc})
    out["tau_form"] = ("time", ref.blended.tau_form, {"units": "N m-2", "long_name": "form drag"})
    out["t_ice"] = (
        "time",
        ref.ice.t_surface,
        {"units": "degC", "long_name": "ice surface temperature, min(T_air, -1.8 C)"},
    )
    for flag, values in ref.qc.items():
        out[f"qc_{flag}"] = ("time", values.astype("i1"),
                             {"long_name": cs.results.QC_FLAGS.get(flag, flag),
                              "flag_values": "0, 1", "flag_meanings": "not_set set"})  # fmt: skip
    if "seb" in runs:
        s = runs["seb"]
        for nm, values in (
            ("hsb_blend_seb", s.blended.hsb),
            ("hlb_blend_seb", s.blended.hlb),
            ("tau_blend_seb", s.blended.tau),
            ("t_ice_seb", s.ice.t_surface),
        ):
            out[nm] = ("time", values, {"long_name": f"{nm[:-4]} with the energy-balance ice "
                       "surface temperature (ERA5 longwave; warm end of a range)"})  # fmt: skip
    out.attrs = {
        "title": "ROVER 2025 Kronprins Haakon turbulent fluxes recomputed with coare-seaice",
        "source": f"coare-seaice {cs.__version__}, coare_seaice.ice_fluxes (simple reference)",
        "input_file": str(SHIP_NC),
        "heights_m": f"zu={ZU}, zt={ZT}, zq={ZQ}",
        "z0_ice_m": Z0_ICE,
        "longwave_source": lw_source,
        "caveats": (
            "ERA5 downwelling longwave is probably too high (ERA5 overestimates cloud). The "
            "non-ice fraction is treated as open water at the SST; across the MIZ this "
            "probably overestimates heat fluxes (IGP aircraft evaluation). Ice roughness is "
            "regional and uncertain."
        ),
        "restricted": "ROVER working group only",
    }
    out.to_netcdf(OUT_NC)
    log.info(f"\nWrote {OUT_NC}")

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
pairs = [k for k in ("tau_blend", "hsb_blend", "hlb_blend") if k in old]
if pairs:
    fig, axes = plt.subplots(len(pairs) + 1, 1, figsize=(10, 2.2 * (len(pairs) + 1)),
                             sharex=True)  # fmt: skip
    for ax, k in zip(axes, pairs, strict=False):
        ax.plot(df.index, old[k], color="grey", lw=0.8, label=f"old {k}")
        ax.plot(df.index, new[k], color="tab:red", lw=0.6, label=f"new {k}")
        ax.set_ylabel(k)
        ax.legend(loc="upper left", frameon=False)
    axes[-1].plot(df.index, sic, color="k", lw=0.4)
    axes[-1].scatter(df.index, sic, c=sic, cmap=SIC_CMAP, vmin=0, vmax=1, s=2, label="AMSR2 SIC")
    axes[-1].set_ylabel("SIC")
    axes[-1].legend(loc="upper left", frameon=False)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "figures" / "rover_00b_old_vs_new_timeseries.png", dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(1, len(pairs), figsize=(3.8 * len(pairs), 3.4))
    for ax, k in zip(np.atleast_1d(axes), pairs, strict=False):
        sc = ax.scatter(sic, new[k] - old[k], c=sic, cmap=SIC_CMAP, vmin=0, vmax=1, s=4,
                        edgecolors="k", linewidths=0.1, label=f"new - old {k}")  # fmt: skip
        ax.axhline(0, color="r", lw=0.8)
        ax.set(xlabel="AMSR2 SIC", ylabel=f"new - old {k}")
        ax.legend(loc="upper left", frameon=False)
    fig.colorbar(sc, ax=axes, label="SIC", shrink=0.8)
    fig.savefig(OUT_DIR / "figures" / "rover_00b_difference_vs_sic.png", dpi=150,
                bbox_inches="tight")  # fmt: skip
    plt.close(fig)
    log.info(f"Figures in {OUT_DIR / 'figures'}")
