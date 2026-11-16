"""IGP 2018: verify coare-seaice against the Alliance ship fluxes and the aircraft eddy covariance.

Three parts, each a separate check:

A. **Ship, recomputed.** The Alliance 10-min ship record already holds bulk fluxes from
   COARE 3.0a (written in 2020). They are recomputed here with coare-seaice from the same
   measured inputs and compared. This checks the open-water code path and the input
   handling (units, heights), not the physics: both are bulk estimates.
B. **Ship bulk against aircraft eddy covariance near the ship.** Low-level aircraft runs
   within a distance of the ship (lenient: tens of km, same time) are paired with the
   ship 10-min window nearest in time. The ship bulk fluxes are computed twice: as open
   water (what the ship measured over), and blended with the ice fraction the aircraft
   saw on that run (a "loosely similar area" estimate). The aircraft was not over the
   ship, so this is a loose test; differences by distance are reported.
C. **Aircraft-driven bulk against aircraft eddy covariance, by ice fraction.** Every
   quality-controlled low-level run over the sea, wherever the ship was. The bulk fluxes
   are computed from the aircraft's own mean measurements, so the inputs and the measured
   fluxes describe the same air and surface. This is the cleaner test of how the bulk
   method performs across ice fractions. Two configurations: the simple reference (two
   tiles, ice at min(T_air, -1.8 C)) and one aggregate surface temperature from the
   aircraft radiometer (blend="single_temperature").

All three use `coare_seaice.ice_fluxes`. Nothing is tuned to these data.

Methods details (written to the log and to OUT_DIR/igp_verification_methods.md):
- Ship inputs: measured wind (LIDAR, 40 m), air temperature and relative humidity
  (WeatherPack or calibrated HATPRO; height set by ZT/ZQ below), pressure, SST (bow
  thermometer, CTD-filled). Windows in harbour, or with any input missing, are skipped.
  The stored COARE 3.0a fluxes exist for some windows with a missing measured input
  (filled in the 2020 processing); those windows are not compared.
- Ship position at an aircraft run time: linear interpolation of the 10-min ship
  latitude and longitude in time. Distance: great circle.
- Aircraft runs: the IGP long-run (9 km) database. Over the sea only (orog < 10 m), low
  level (radar altitude <= MAX_ALT for part C; <= MAX_ALT_B for part B, which is lenient
  because few runs near the ship were flown low), wind present. Terrain height (orog)
  is missing for some flights; those runs are kept (the database lists them over sea).
  Quality: z0_flag <= 1 for stress,
  z0T_flag <= 1 for sensible heat, z0q_flag <= 1 for latent heat, as the database readme
  recommends. Ice fraction: icefractionalb (albedo thresholds 0.15 and 0.85).
- Aircraft-driven bulk inputs: wsp, T, rhwater_BUCK at the run's radar altitude; mean
  sea-level pressure; OSTIA SST for the open-water tile; radiometer Tsfc for the
  single-temperature configuration. Radiation is not needed (no cool skin, jcool = 0).
- Aircraft eddy-covariance fluxes are at flight level, not at the surface; no
  flux-divergence correction is applied.
- Statistics: count, bias (bulk minus observed) as mean and median, RMSE, correlation.

Paths: environment variables IGP_DATA (folder with the two input files) and
COARE_SEAICE_OUT (output folder). Log: <OUT_DIR>/logs/igp_verify_ship_aircraft.log
"""

from __future__ import annotations

import logging
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
DATA = Path(os.environ.get("IGP_DATA", "."))
SHIP_NC = DATA / "IGP_Alliance_Surface_Met_QC3_10min.nc"
AIRCRAFT_TXT = DATA / "IGP_flights_database_obs_lr.txt"
OUT_DIR = Path(os.environ.get("COARE_SEAICE_OUT", "outputs"))

ZU = 40.0  # ship wind height (m): "Wind Speed at 40m altitude from WindCube LIDAR"
# Ship temperature and humidity height (m). NOT in the file metadata. 20 m is the height
# that best reproduces the stored COARE 3.0a fluxes (scan of 5-40 m, see the log); please
# replace with the documented height.
ZT = ZQ = 20.0
COARE = "3.5"
JCOOL_SHIP = 1  # cool skin on for the ship, with default radiation (see RS_SHIP, RL_SHIP)
RS_SHIP, RL_SHIP = 150.0, 370.0  # no radiation in the ship file: pycoare-style defaults
MAX_DIST_KM = 50.0  # part B: aircraft run within this distance of the ship
MAX_ALT = 60.0  # part C: aircraft runs at or below this radar altitude (m)
# Part B is lenient: the runs near the ship were mostly flown higher (lowest run of flight
# 295 was 77 m), so allow higher runs there and report results by altitude.
MAX_ALT_B = 150.0
QC_MAX = 1  # database flags: 0 good, 1 ok
ICE_BINS = [-0.001, 0.0001, 0.25, 0.5, 0.75, 0.9999, 1.0]
ICE_LABELS = ["0", "0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1", "1"]
REPROCESS = True
# ---------------------------------------------------------------------------

warnings.filterwarnings("ignore", category=RuntimeWarning)
NAME = "igp_verify_ship_aircraft"
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
FLUXES = [("tau", "windstress", "N m-2", "z0_flag"), ("hsb", "sh", "W m-2", "z0T_flag"),
          ("hlb", "lh", "W m-2", "z0q_flag")]  # fmt: skip
log.info(f"{NAME} | coare-seaice {cs.__version__} | COARE {COARE}")
log.info(f"Ship: {SHIP_NC}\nAircraft: {AIRCRAFT_TXT}\nOutput: {OUT_DIR}\n")


def stats(bulk, obs):
    ok = np.isfinite(bulk) & np.isfinite(obs)
    b, o = np.asarray(bulk)[ok], np.asarray(obs)[ok]
    if ok.sum() < 2:
        return dict(n=int(ok.sum()), bias=np.nan, med=np.nan, rmse=np.nan, r=np.nan)
    d = b - o
    return dict(n=int(ok.sum()), bias=d.mean(), med=np.median(d), rmse=np.sqrt((d**2).mean()),
                r=np.corrcoef(b, o)[0, 1])  # fmt: skip


def fmt(s, unit):
    p = 4 if unit == "N m-2" else 1
    return (f"n={s['n']:4d}  bias {s['bias']:+.{p}f}  median {s['med']:+.{p}f}  "
            f"RMSE {s['rmse']:.{p}f}  r {s['r']:.2f}")  # fmt: skip


def great_circle_km(lat1, lon1, lat2, lon2):
    la1, lo1, la2, lo2 = map(np.radians, (lat1, lon1, lat2, lon2))
    a = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return 2 * 6371.0 * np.arcsin(np.sqrt(a))


# ---------------------------------------------------------------------------
# A. Ship: recompute and compare with the stored COARE 3.0a fluxes
# ---------------------------------------------------------------------------
log.info("=" * 78 + "\nA. Ship fluxes recomputed with coare-seaice vs stored COARE 3.0a\n")
ship = xr.open_dataset(SHIP_NC).to_dataframe()
# Units: the file labels SST 'K' and RH '1', but the values are deg C and per cent.
s_u = ship.measured_wind_speed.to_numpy()
s_t = ship.measured_air_temperature.to_numpy() - 273.15
s_rh = ship.measured_relative_humidity.to_numpy()
s_p = ship.measured_air_pressure.to_numpy()
s_sst = ship.sea_surface_temperature.to_numpy()
s_lat = ship.latitude_nmea.to_numpy()
assert np.nanmax(s_sst) < 40 and np.nanmax(s_rh) > 2, "unexpected SST or RH units"
use = (ship.in_harbour_flag.to_numpy() == 0) & np.all(
    [np.isfinite(x) for x in (s_u, s_t, s_rh, s_p, s_sst, s_lat)], axis=0
)
stored = {"tau": ship.windstress.to_numpy(), "hsb": ship.sensible_heat_flux.to_numpy(),
          "hlb": ship.latent_heat_flux.to_numpy()}  # fmt: skip
n_stored = int(np.isfinite(stored["hsb"]).sum())
log.info(f"10-min windows: {len(ship)}; with stored fluxes: {n_stored}; usable for recompute "
         f"(at sea, all inputs measured): {use.sum()}; both: "
         f"{(use & np.isfinite(stored['hsb'])).sum()}")  # fmt: skip


def ship_run(sic, zt=ZT, idx=None):
    i = np.where(use)[0] if idx is None else idx
    return cs.ice_fluxes(
        s_u[i], s_t[i], s_rh[i], s_p[i], s_sst[i], sic, zu=ZU, zt=zt, zq=zt, lat=s_lat[i],
        lw_down=RL_SHIP, sw_down=RS_SHIP, jcool=JCOOL_SHIP, coare=COARE,
    )  # fmt: skip


log.info(f"\nHeight check (zt = zq scanned; zu = {ZU} m; open water; median new - stored):")
idx = np.where(use)[0]
for zt in (10.0, 15.0, 20.0, 25.0, 30.0):
    r = ship_run(0.0, zt)
    d = r.water.hsb - stored["hsb"][idx]
    log.info(f"  zt={zt:4.0f} m: sensible {np.nanmedian(d):+.2f} W m-2, "
             f"latent {np.nanmedian(r.water.hlb - stored['hlb'][idx]):+.2f} W m-2")  # fmt: skip
log.info(f"Using zt = zq = {ZT} m (not documented in the file: confirm).\n")

res_ship = ship_run(0.0)
new_ship = {"tau": np.full(len(ship), np.nan), "hsb": np.full(len(ship), np.nan),
            "hlb": np.full(len(ship), np.nan)}  # fmt: skip
for k in new_ship:
    new_ship[k][idx] = getattr(res_ship.water, k)
for k, _, unit, _ in FLUXES:
    log.info(f"  {k:4s} ({unit}) new - stored: {fmt(stats(new_ship[k], stored[k]), unit)}")
log.info(
    "Expected: close but not identical. COARE 3.5 changed the wind-speed dependence of the "
    "Charnock parameter and the roughness Reynolds number relation from 3.0a (Edson et al. "
    "2013), so stress (and through it the heat fluxes) differ most at high wind.\n"
)
ow = res_ship.water.extra["open_water"]
out_a = pd.DataFrame(
    {"lat": s_lat, "lon": ship.longitude_nmea.to_numpy(), **{f"stored_{k}": v for k, v in
     stored.items()}, **{f"new_{k}": v for k, v in new_ship.items()}}, index=ship.index,
)  # fmt: skip
out_a.loc[out_a.index[idx], "u10n_new"] = ow.u_n_rf
out_a["u10n_stored"] = ship.u10N_coare.to_numpy()
out_a.to_csv(OUT_DIR / "data" / "igp_ship_recomputed_fluxes.csv")

fig, axes = plt.subplots(2, 3, figsize=(11, 6.5))
for j, (k, _, unit, _) in enumerate(FLUXES):
    ax = axes[0, j]
    ax.scatter(stored[k], new_ship[k], s=3, c="k", alpha=0.4, label="10-min windows")
    lim = np.nanpercentile(np.r_[stored[k], new_ship[k]], [0.5, 99.5])
    ax.plot(lim, lim, "r-", lw=1, label="1:1")
    ax.set(xlabel=f"stored COARE 3.0a {k} ({unit})", ylabel=f"coare-seaice {k} ({unit})")
    ax.legend(loc="upper left", frameon=False)
    ax = axes[1, j]
    ax.plot(ship.index, new_ship[k] - stored[k], "k.", ms=1.5, label="new - stored")
    ax.axhline(0, color="r", lw=0.8)
    ax.set(ylabel=f"difference ({unit})")
    ax.legend(loc="upper left", frameon=False)
    ax.tick_params(axis="x", labelrotation=30)
fig.tight_layout()
fig.savefig(OUT_DIR / "figures" / "igp_A_ship_recomputed_vs_stored.png", dpi=150)
plt.close(fig)

# ---------------------------------------------------------------------------
# Aircraft runs
# ---------------------------------------------------------------------------
air = pd.read_csv(AIRCRAFT_TXT, sep=r"\s+").copy()
air["time"] = pd.to_datetime(air.yyyymmdd.astype(int).astype(str)) + pd.to_timedelta(
    air.meantime, unit="s"
)
sea = ~(air.orog >= 10) & air.wsp.notna() & air.sstOSTIA.notna()  # orog NaN counts as sea
low = sea & (air.altradgps <= MAX_ALT)
log.info("=" * 78)
log.info(f"Aircraft runs: {len(air)}; over the sea, radar altitude <= {MAX_ALT:.0f} m, with "
         f"wind and OSTIA SST: {low.sum()}")  # fmt: skip
for _, obs, _, flag in FLUXES:
    log.info(f"  {obs:10s}: passes {flag} <= {QC_MAX} in {((air[flag] <= QC_MAX) & low).sum()}")

sv = ship[["latitude_nmea", "longitude_nmea"]].dropna()
tsec = sv.index.to_numpy().astype("datetime64[s]").astype(float)
at = air.time.to_numpy().astype("datetime64[s]").astype(float)
slat = np.interp(at, tsec, sv.latitude_nmea.to_numpy(), left=np.nan, right=np.nan)
slon = np.interp(at, tsec, sv.longitude_nmea.to_numpy(), left=np.nan, right=np.nan)
air["dist_ship_km"] = great_circle_km(air.lat, air.lon, slat, slon)

# ---------------------------------------------------------------------------
# B. Ship bulk vs aircraft eddy covariance near the ship
# ---------------------------------------------------------------------------
log.info("\n" + "=" * 78)
log.info(f"B. Ship bulk vs aircraft eddy covariance, runs within {MAX_DIST_KM:.0f} km of the ship"
         f" and at or below {MAX_ALT_B:.0f} m\n")  # fmt: skip
near = air[sea & (air.altradgps <= MAX_ALT_B) & (air.dist_ship_km <= MAX_DIST_KM)].copy()
# nearest ship window in time (10 min windows; keep matches within 10 min)
pos = np.searchsorted(ship.index.to_numpy(), near.time.to_numpy())
pos = np.clip(pos, 1, len(ship) - 1)
prev = ship.index.to_numpy()[pos - 1]
nxt = ship.index.to_numpy()[pos]
pick = np.where(np.abs(near.time.to_numpy() - prev) <= np.abs(nxt - near.time.to_numpy()),
                pos - 1, pos)  # fmt: skip
dt_min = np.abs(ship.index.to_numpy()[pick] - near.time.to_numpy()) / np.timedelta64(1, "m")
near["ship_i"] = pick
ok_b = (dt_min <= 10) & use[pick]
near = near[ok_b]
log.info(f"Runs within {MAX_DIST_KM:.0f} km: {int(ok_b.size)}; with a usable ship window within "
         f"10 min: {len(near)} (flights {sorted(near.flightno.unique())})")  # fmt: skip
log.info(f"Distance to ship: median {near.dist_ship_km.median():.1f} km, range "
         f"{near.dist_ship_km.min():.1f}-{near.dist_ship_km.max():.1f} km; aircraft ice fraction "
         f"median {near.icefractionalb.median():.2f}, "
         f"max {near.icefractionalb.max():.2f}")  # fmt: skip
ii = near.ship_i.to_numpy()
b_ow = ship_run(0.0, idx=ii)
b_ice = ship_run(near.icefractionalb.to_numpy(), idx=ii)
for k, obs, unit, flag in FLUXES:
    q = (near[flag] <= QC_MAX).to_numpy()
    o = near[obs].to_numpy()
    log.info(f"\n  {k} ({unit}) vs aircraft {obs} ({flag} <= {QC_MAX}):")
    log.info(f"    ship open water        {fmt(stats(getattr(b_ow.blended, k)[q], o[q]), unit)}")
    log.info(f"    ship + run ice frac.   {fmt(stats(getattr(b_ice.blended, k)[q], o[q]), unit)}")
    for lo, hi in ((0, 10), (10, 25), (25, MAX_DIST_KM)):
        m = q & (near.dist_ship_km >= lo).to_numpy() & (near.dist_ship_km < hi).to_numpy()
        log.info(f"    {lo:>2.0f}-{hi:<2.0f} km (ice frac.)  "
                 f"{fmt(stats(getattr(b_ice.blended, k)[m], o[m]), unit)}")  # fmt: skip
    for lo, hi in ((0, 60), (60, MAX_ALT_B)):
        m = q & (near.altradgps >= lo).to_numpy() & (near.altradgps < hi).to_numpy()
        log.info(f"    {lo:>3.0f}-{hi:<3.0f} m alt. (ice f.) "
                 f"{fmt(stats(getattr(b_ice.blended, k)[m], o[m]), unit)}")  # fmt: skip
    for f in sorted(near.flightno.unique()):
        m = q & (near.flightno == f).to_numpy()
        log.info(f"    flight {f} (ice frac.)  "
                 f"{fmt(stats(getattr(b_ice.blended, k)[m], o[m]), unit)}")  # fmt: skip
log.info(
    "\nCaveats: the aircraft was tens of km from the ship and sampled a different surface; "
    "aircraft fluxes are at flight level; the ship temperature/humidity height is inferred."
)
out_b = near[["time", "flightno", "legno", "runno", "lat", "lon", "dist_ship_km", "altradgps",
              "icefractionalb", "windstress", "sh", "lh", "z0_flag", "z0T_flag",
              "z0q_flag"]].copy()  # fmt: skip
for k in ("tau", "hsb", "hlb"):
    out_b[f"ship_ow_{k}"] = getattr(b_ow.blended, k)
    out_b[f"ship_ice_{k}"] = getattr(b_ice.blended, k)
out_b.to_csv(OUT_DIR / "data" / "igp_B_aircraft_near_ship.csv", index=False)

fig, axes = plt.subplots(1, 3, figsize=(11, 3.6))
for j, (k, obs, unit, flag) in enumerate(FLUXES):
    ax = axes[j]
    q = (near[flag] <= QC_MAX).to_numpy()
    sc = ax.scatter(near[obs][q], getattr(b_ice.blended, k)[q], c=near.icefractionalb[q],
                    cmap=SIC_CMAP, vmin=0, vmax=1, s=18, edgecolors="k", linewidths=0.4,
                    label="ship bulk, run ice fraction")  # fmt: skip
    ax.scatter(near[obs][q], getattr(b_ow.blended, k)[q], marker="x", c="grey", s=12,
               label="ship bulk, open water")  # fmt: skip
    lim = np.nanpercentile(np.r_[near[obs][q], getattr(b_ow.blended, k)[q]], [0, 100])
    ax.plot(lim, lim, "r-", lw=1, label="1:1")
    ax.set(xlabel=f"aircraft {obs} ({unit})", ylabel=f"ship bulk {k} ({unit})")
    ax.legend(loc="lower right", frameon=False, fontsize=7)
fig.colorbar(sc, ax=axes, label="aircraft ice fraction", shrink=0.8)
fig.savefig(OUT_DIR / "figures" / "igp_B_ship_bulk_vs_aircraft_near_ship.png", dpi=150,
            bbox_inches="tight")  # fmt: skip
plt.close(fig)

# ---------------------------------------------------------------------------
# C. Aircraft-driven bulk vs aircraft eddy covariance, by ice fraction
# ---------------------------------------------------------------------------
log.info("\n" + "=" * 78 + "\nC. Aircraft-driven bulk vs aircraft eddy covariance, by ice "
         "fraction (all low-level runs over the sea)\n")  # fmt: skip
c = air[low].copy()
c_args = (c.wsp.to_numpy(), c["T"].to_numpy() - 273.15, c.rhwater_BUCK.to_numpy(),
          c.mslp.to_numpy(), c.sstOSTIA.to_numpy() - 273.15,
          c.icefractionalb.to_numpy())  # fmt: skip
z = c.altradgps.to_numpy()
c_kw = dict(zu=z, zt=z, zq=z, lat=c.lat.to_numpy(), jcool=0, coare=COARE)
configs = {
    "reference (2 tiles, ice at min(T,-1.8 C))": cs.ice_fluxes(*c_args, **c_kw),
    "single temperature (radiometer Tsfc)": cs.ice_fluxes(
        *c_args, blend="single_temperature", t_surface=c.Tsfc.to_numpy() - 273.15, **c_kw
    ),
}
c["ice_bin"] = pd.cut(c.icefractionalb, ICE_BINS, labels=ICE_LABELS)
rows = []
for cname, r in configs.items():
    log.info(f"[{cname}]")
    for k, obs, unit, flag in FLUXES:
        q = (c[flag] <= QC_MAX).to_numpy()
        bulk, o = getattr(r.blended, k), c[obs].to_numpy()
        log.info(f"  {k} ({unit}), all bins:      {fmt(stats(bulk[q], o[q]), unit)}")
        for lab in ICE_LABELS:
            m = q & (c.ice_bin == lab).to_numpy()
            s = stats(bulk[m], o[m])
            rows.append(dict(config=cname, flux=k, ice_bin=lab, **s))
            log.info(f"    ice {lab:>9s}:            {fmt(s, unit)}")
    log.info("")
pd.DataFrame(rows).to_csv(OUT_DIR / "data" / "igp_C_stats_by_ice_fraction.csv", index=False)
log.info(
    "Caveats: aircraft fluxes are at flight level (15-60 m), so flux divergence is not "
    "removed; open-sea scalar coefficients in these aircraft data are known to be below "
    "COARE's (Elvidge et al. 2021); one radiometer temperature represents the whole run."
)

fig, axes = plt.subplots(1, 3, figsize=(11, 3.6))
x = np.arange(len(ICE_LABELS))
for j, (k, obs, unit, flag) in enumerate(FLUXES):
    ax = axes[j]
    q = (c[flag] <= QC_MAX).to_numpy()
    obs_med = [np.nanmedian(c[obs].to_numpy()[q & (c.ice_bin == lab).to_numpy()])
               if (q & (c.ice_bin == lab).to_numpy()).any() else np.nan
               for lab in ICE_LABELS]  # fmt: skip
    ax.plot(x, obs_med, "ko-", label="aircraft eddy covariance")
    for (cname, r), col in zip(configs.items(), ("tab:red", "tab:blue"), strict=True):
        b = getattr(r.blended, k)
        med = [np.nanmedian(b[q & (c.ice_bin == lab).to_numpy()])
               if (q & (c.ice_bin == lab).to_numpy()).any() else np.nan
               for lab in ICE_LABELS]  # fmt: skip
        ax.plot(x, med, "s--", color=col, label=cname.split(" (")[0])
    ax.set_xticks(x, ICE_LABELS)
    ax.set(xlabel="ice fraction (icefractionalb)", ylabel=f"median {k} ({unit})")
    ax.legend(loc="best", frameon=False, fontsize=7)
fig.tight_layout()
fig.savefig(OUT_DIR / "figures" / "igp_C_bulk_vs_aircraft_by_ice_fraction.png", dpi=150)
plt.close(fig)

# ---------------------------------------------------------------------------
# Methods file
# ---------------------------------------------------------------------------
(OUT_DIR / "igp_verification_methods.md").write_text(
    "# IGP verification of coare-seaice: methods\n\n"
    + (__doc__ or "").split("Paths:")[0].strip()
    + f"\n\n## Settings used\n\n- coare-seaice {cs.__version__}, COARE {COARE}\n"
    f"- ship heights zu = {ZU} m, zt = zq = {ZT} m (zt inferred, to confirm)\n"
    f"- ship cool skin {'on' if JCOOL_SHIP else 'off'}, radiation {RS_SHIP}/{RL_SHIP} W m-2 "
    "(none in the ship file)\n"
    f"- aircraft runs: radar altitude <= {MAX_ALT} m, over the sea, database flags <= {QC_MAX}\n"
    f"- part B: within {MAX_DIST_KM} km of the ship, ship window within 10 min\n"
)
log.info(f"\nWrote figures, CSV tables and methods to {OUT_DIR}")
