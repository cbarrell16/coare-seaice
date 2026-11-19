# Scripts

Data retrieval and campaign runs. Not part of the installed package. Each script has a
config block with `REPROCESS = True` by default and writes a log to
`outputs/logs/<name>.log`. Data are read from the directory named by the environment
variable `COARE_SEAICE_DATA`.

| Script | Runs on | What it does |
|---|---|---|
| `p5_reprocess_rover_legacy.py` | UEA machine | P5: ROVER fluxes with the legacy wrapper, checked against v2 |
| `rover_00b_fluxes_coare_seaice.py` | UEA machine | ROVER fluxes with `ice_fluxes` (simple reference, plus an energy-balance run), written to a separate file and compared with the old fluxes in the 00 output |
| `igp_verify_ship_aircraft.py` | anywhere (`IGP_DATA`) | IGP: ship fluxes recomputed against the stored COARE 3.0a fluxes; ship bulk against aircraft eddy covariance near the ship; aircraft-driven bulk against aircraft eddy covariance by ice fraction. Results: `docs/verification/igp_2018.md` |
