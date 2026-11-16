# Task backlog

Each task is small, has acceptance criteria, names the model tier to use and lists what it
depends on. One task per branch and pull request (`task/<id>-<short-name>`).

**Tier:** `sonnet` for implementation, tests, docs and data work; `opus` for physics review
and careful reasoning; `fable` only for design decisions.
**Where:** `cloud` runs in a Claude Code web session; `UEA` needs Chris's machine and data.
**Chris:** tasks marked 🔶 need Chris's judgement or action before or during the task.

Every task that touches equations, constants or results ends with a `physics-reviewer`
(`opus`) pass before merging.

## Release 0.2.0 — v2 ported, behaviour-identical

**Status (14 Nov 2026):** P1–P4 done on branch `task/P1-P4-port-v2`, awaiting Chris's review. The stored reference uses 300 synthetic cases per run (15 runs), not 2000, to keep the fixtures small; the port was also checked against v2 on the full ROVER record outside the repository. P5 is open.

| ID | Task | Tier | Depends | Where | Acceptance criteria |
|---|---|---|---|---|---|
| P1 | Port the ice physics from v2 into modules (`thermo`, `stability`, `roughness`, `form_drag`, `most`, `surface_temperature`) with no change in behaviour | opus | – | cloud | Each function gives the same output as its v2 counterpart on stored reference inputs (relative tolerance 1e-12). No NOAA or Blomquist code copied. |
| P2 | `openwater.py`: adapter over the public `pycoare` interface returning the quantities the ice code needs (fluxes, final z0, 10 m neutral wind, air density, humidity) | opus | – | cloud | Open-water fluxes equal plain `pycoare` exactly. Works for COARE 3.5 and 3.6. |
| P3 | `legacy.coare35vn_ice`: the v2 call signature and 56-column array, built on P1 and P2 | opus | P1, P2 | cloud | Columns 38–55 match stored v2 outputs to 1e-3 relative on 2000 synthetic cases for every `t_ice_method`; columns 0–37 match v2 within the documented `pycoare` differences; any column that cannot be reproduced is listed. |
| P4 | Store v2 reference outputs as a test fixture (outputs only, generated from synthetic inputs) | sonnet | – | cloud | Fixture under 300 kB with a README recording how it was made. |
| P5 | 🔶 Reprocess the shared ROVER flux file with `scripts/p5_reprocess_rover_legacy.py` (`air_freezing`, v2 bug fixes only) and check it against v2 | Chris | P3 | UEA | Log shows ice-only columns identical to v2, open-water columns within about 0.1 W m⁻², no NaN mismatches. Chris signs off the port. See `docs/guides/install_and_reprocess.md`. |

## Release 0.3.0 — audit fixes, one per pull request

Each fix states in `CHANGELOG.md` how much results move on the reference fixture.

| ID | Task | Tier | Depends | Where | Acceptance criteria |
|---|---|---|---|---|---|
| F1 | Replace the hard-coded stability clip with a named, documented limit and a QC flag based on bulk Richardson number (0.2, Grachev et al. 2007) | opus | P3 | cloud | Default results unchanged; flag set where the limit acts; test covers both sides. |
| F2 | Use the temperature-dependent air viscosity in the Andreas (1987) roughness Reynolds number | sonnet | P1 | cloud | Matches the COARE viscosity formula; change in H reported (expected 1–2 %). |
| F3 | Guard the Andreas (1987) polynomial beyond R* = 1000, with a flag | sonnet | P1 | cloud | No extrapolation; flagged; test at R* = 5000. |
| F4 | Column model: output its own conductive flux; fix the wrong values when the surface is warmer than the base | sonnet | P1 | cloud | Test with surface at −1 °C and base at −1.8 °C gives a finite, correctly signed flux. |
| F5 | Column model: initialise with its own (salinity-dependent) properties; reset `melt` each pass; decide and document array `h_ice` behaviour | sonnet | P1 | cloud | Steady-state test shows no drift after initialisation. |
| F6 | Add a convergence diagnostic to the ice surface-layer iteration and update u* in the final pass | sonnet | P1 | cloud | Diagnostic output; results change by less than 0.01 W m⁻² on the fixture. |
| F7 | No silent fills: missing longwave or other inputs propagate as NaN with a flag | sonnet | P3 | cloud | NaN-in gives NaN-out test for every input. |
| F8 | Accept array heights and any real number type for heights | sonnet | P1 | cloud | Per-sample heights give the same result as looping over scalars. |
| F9 | Correct constants and citations: fresh-ice heat capacity 2110 J kg⁻¹ K⁻¹ (Bitz & Lipscomb 1999); conductivity cited to Untersteiner (1964); form-drag term cited as Lüpkes et al. (2012) Eq. 51 | sonnet | P1 | cloud | Changelog records the size of the change. |
| F10 | Add gustiness (unstable) and windless (stable) terms over ice after Andreas et al. (2010) Eqs. 2.4–2.6, off by default | opus | P1 | cloud | Equations quoted in the docstring; off reproduces previous results. |

## Release 0.4.0 — tiles, stability options, new interface

| ID | Task | Tier | Depends | Where | Acceptance criteria |
|---|---|---|---|---|---|
| T1 | `results.py`: named outputs per tile and blended, QC flags, `to_xarray()` | sonnet | P3 | cloud | Round-trips to the legacy array for the two-tile case. |
| T2 | Surface-temperature interface with `AirFreezing`, `Prescribed`, `SEB`; `Column` as experimental | opus | P1 | cloud | Same numbers as the legacy methods; a new method can be added without touching other modules. |
| T3 | `tiles.py`: tile definition, `two_tile()`, `three_tile()`, `flux` blending | opus | T1, T2 | cloud | Two tiles equals legacy; three tiles with zero thin-ice fraction equals two tiles. |
| T4 | `single_temperature` blend mode (one aggregate surface temperature, as Elvidge et al. 2021 Eq. 12) | opus | T3 | cloud | Equals `flux` mode when all tiles have the same temperature and stability is neutral. |
| T5 | `api.ice_fluxes`: the one entry point, with defaults that reproduce v2 `air_freezing` | opus | T3 | cloud | Signature as in `docs/design/software_design.md`; legacy wrapper now calls it. |
| T6 | SHEBA stability functions (Grachev et al. 2007 Eqs. 9, 12, 13) as an option | sonnet | P1 | cloud | psi functions satisfy phi = 1 − ζ dψ/dζ numerically; COARE remains default. |
| T7 | Ice roughness options: number, array, Andreas et al. (2010) Eq. 5.1, named settings from Renfrew et al. (2019) Table 1 and Srivastava et al. (2022) Table 2 | sonnet | P1 | cloud | Each named value cites its table; default unchanged. |
| T8 | COARE 3.6 selectable, wave inputs passed through | sonnet | P2 | cloud | Identical to 3.5 results without wave inputs, within `pycoare` differences. |
| T9 | Humidity QC: flag air far from ice saturation over ice (Andreas et al. 2002) | sonnet | T1 | cloud | Flag only; no adjustment. |
| T10 | 🔶 Thin-ice tile defaults: thickness, snow, how the thin fraction is specified | fable | T3, V2 | cloud | Proposal with evidence from the IGP and ACCACIA experiments; Chris decides. |
| U1 | `uncertainty.ensemble`: perturbed-input reruns returning percentiles | sonnet | T5 | cloud | Reproducible with a seed; zero perturbation returns the central result. |

## Release 0.5.0 — verification toolkit

| ID | Task | Tier | Depends | Where | Acceptance criteria |
|---|---|---|---|---|---|
| V1 | `verify/format.py`: the common table format and its validation | sonnet | T1 | cloud | Documented variable list; validation errors name the missing field. |
| V2 | Reader for the IGP and ACCACIA aircraft run databases; small public subset as a fixture | sonnet | V1 | cloud | Reads both files; flags and ice fraction (`icefractionalb`) handled as in the databases readme. |
| V3 | `verify/metrics.py`: bias, median error, RMSE, correlation, slope with bootstrap intervals; stratification by ice fraction, stability, wind, longwave | sonnet | V1 | cloud | Unit tests on synthetic data with known answers. |
| V4 | Experiments E0–E6, E9, E10 on IGP and ACCACIA (see `docs/design/verification_design.md`) | opus | V2, V3, T5 | cloud | Tables and figures regenerate from one command; methods written into the docs. |
| V5 | 🔶 Download MOSAiC and SHEBA tower data; record variables, heights and licences | data-pipeline | – | UEA | Inventory written to the docs; open items updated. |
| V6 | Tower reader (MOSAiC, SHEBA) | sonnet | V1, V5 | cloud + UEA | Written against a synthetic sample; run on real data by Chris. |
| V7 | Experiments E1, E2, E4, E8 on MOSAiC and SHEBA | opus | V6, T6, F10 | UEA | As V4. SHEBA reported as an implementation check, not independent evidence. |
| V8 | Ship reader (ROVER, IGP Alliance) and the ROVER application with uncertainty bands | sonnet | V1, U1 | UEA | ROVER data never leave the UEA machine. |
| V9 | `verify/benchmark`: one command that reruns all available campaigns and reports changes from the previous version | sonnet | V4 | cloud | Output is a table of metric changes; skips campaigns whose data are absent. |
| V10 | Standard figures following the figure rules in `CLAUDE.md` | sonnet | V3 | cloud | White backgrounds, legends inside panels, no explanatory text. |
| V11 | 🔶 Satellite ice surface temperature: retrieval script tested against the live Copernicus catalogue; comparison with IGP radiometer temperature | data-pipeline | V2 | UEA | Product and dates recorded; clear-sky caveat in the docs. |
| V12 | 🔶 Multi-product sea-ice concentration along the ROVER track | data-pipeline | V8 | UEA | Uncertainty band written alongside the central value. |

## Documentation and release

| ID | Task | Tier | Depends | Where | Acceptance criteria |
|---|---|---|---|---|---|
| D1 | Methods documentation, one page per module, equations with paper references | docs-writer | P1 | cloud | `mkdocs build --strict` passes; every equation cites its source. |
| D2 | User guide with a worked example on synthetic data | docs-writer | T5 | cloud | Example runs in CI. |
| D3 | Migration guide from `coare35vn_ice` | docs-writer | T5 | cloud | Chris's notebook calls map one to one. |
| R1a | 🔶 Post the `pycoare` bug report (`docs/drafts/pycoare_issue_rh_mutation.md`) as an issue on pyCOARE/coare | Chris | – | – | Issue posted; link recorded here. |
| R1b | 🔶 Send the courtesy note to the `pycoare` maintainer (`docs/drafts/pycoare_maintainer_note.md`), linking R1a | Chris | R1a | – | Reply recorded; name confirmed before any PyPI registration. |
| R2 | 🔶 Release 1.0.0: tag, Zenodo archive, PyPI | fable + Chris | all | – | Only on Chris's explicit go-ahead. |

## Open questions that block tasks

- **T10:** no source yet for the thin-ice fraction on ROVER.
- **V5:** MOSAiC and SHEBA licences not yet read.
- **V8:** ROVER wind is flagged in many ice windows; decide whether a correction exists.
