# COARE-ice: Phase B software design

*13 Nov 2026. Claude (Fable) for Chris Barrell. Checkpoint B of the project brief. Nothing has been created yet; this is the design for your approval.*

Each software-engineering choice has a short "why", as you asked.

## Design in one paragraph

A small Python package, separate from `pycoare` and depending on it for the open-water fluxes. The surface is a list of **tiles** (open water, thick ice, thin ice). Each tile has a fraction, a roughness, a surface-temperature method and a stability treatment. One function computes fluxes per tile and blends them. The simplest configuration (two tiles, ice at min(T_air, −1.8 °C)) reproduces v2 and stays the reference that every richer option is compared against. A separate optional part of the package holds the verification tools.

## 1. Name and layout

**Name:** `coare-seaice` to install, `coare_seaice` to import. It was free on PyPI when I checked today. It says what it is without looking like an official part of `pycoare`. I will not register it until you have contacted the `pycoare` maintainer.

**Layout** ("src layout": the code sits under `src/`, so tests always run against the installed package, not a stray local copy):

```
coare-seaice/
  pyproject.toml          package definition, dependencies, tool settings
  README.md  LICENSE  CITATION.cff  CHANGELOG.md  CLAUDE.md
  .gitignore  .pre-commit-config.yaml
  .github/workflows/      ci.yml, upstream.yml
  .claude/agents/         subagent definitions
  src/coare_seaice/
    __init__.py           public names and version
    api.py                ice_fluxes(): the one entry point
    openwater.py          thin adapter over pycoare (public API only)
    most.py               surface-layer solver for ice tiles, with QC flags
    stability.py          stability functions: COARE and SHEBA
    roughness.py          z0 options; Andreas (1987) scalar roughness
    form_drag.py          Lüpkes (2012) / Elvidge (2016)
    thermo.py             humidity, saturation, constants
    surface_temperature/  one interface; simple, prescribed, seb, column (experimental)
    tiles.py              Tile definition, standard tile sets, blending
    uncertainty.py        perturbed-input ensembles
    results.py            named outputs, QC flags, conversion to xarray
    legacy.py             coare35vn_ice(): the v2 call and 56-column array
    verify/               optional: data format, readers, colocation, metrics, figures, benchmark
  tests/                  unit/, regression/, property/, data/ (small synthetic files)
  docs/                   methods, caveats, API reference
  scripts/                data retrieval and campaign runs (not installed)
  notebooks/              build_*.py builders only
```

Modules are small and single-purpose so that `roughness`, `form_drag`, `stability` and `surface_temperature` could later move into `pycoare/ice/` unchanged.

## 2. Core concepts

- **Tile.** A surface type with: name, fraction, roughness option, scalar-roughness scheme, surface-temperature method, and thickness and snow where relevant. Open water is a tile whose fluxes come from `pycoare`.
- **Surface-temperature method.** One interface with a single job: given the forcing and a tile, return the surface temperature and its diagnostics. `AirFreezing`, `Prescribed`, `SEB` and `Column` all implement it. Adding a method does not touch the rest of the code.
- **Blend mode.** Two modes, because the IGP first look showed both are needed:
  - `flux`: each tile has its own temperature and stability; fluxes are area-weighted, plus form drag. This is v2.
  - `single_temperature`: all tiles share one aggregate surface temperature, as in Elvidge et al. (2021). For use when a radiometer measures the footprint mean.
- **QC flags.** Every output carries flags, never silent adjustments: stability outside the validated range (bulk Richardson number above 0.2), light wind, measurement height large relative to the Obukhov length, scalar-roughness formula out of range, humidity far from ice saturation, melting surface, missing input.
- **Missing data.** NaN in gives NaN out. No default fills.

## 3. Public API

```python
import coare_seaice as cs

res = cs.ice_fluxes(
    u=u, t=t, rh=rh, p=p, sst=sst, sic=sic,       # arrays or scalars
    zu=36, zt=22, zq=22,                           # scalars or arrays
    lat=lat, lw_down=lw, sw_down=sw,
    coare="3.5",                                   # or "3.6", with wave inputs passed through
    tiles=cs.two_tile(),                           # default; cs.three_tile(thin_fraction=...) for thin ice
    t_ice=cs.AirFreezing(),                        # or cs.SEB(h_ice=1.0, h_snow=0.2), cs.Prescribed(values)
    z0_ice=5e-4,                                   # number, array, or a named option
    stability="coare",                             # or "sheba"
    blend="flux",                                  # or "single_temperature"
)

res.blended.hsb          # area-mean sensible heat flux
res.water.tau            # open-water tile
res.tiles["thick_ice"].t_surface
res.qc                   # flags
res.to_xarray()          # everything as a labelled Dataset
```

- **Named outputs, not column numbers.** Names follow `pycoare` (`tau`, `hsb`, `hlb`, `usr`, …) so the two read alike.
- **Sign convention:** heat fluxes positive upward, as in COARE and v2.
- **Heights may be arrays.** The aircraft data need this; v2 allows scalars only.
- **Defaults reproduce v2 with `air_freezing`.** The default surface-temperature method is decided after Phase C.
- **Uncertainty:** `cs.ensemble(ice_fluxes_kwargs, perturb={"sic": 0.05, "lw_down": 20, ...}, n=200)` reruns with perturbed inputs and returns percentiles. Simple first: Monte Carlo, no analytic propagation.

## 4. Dependencies

| Dependency | Pin | Why |
|---|---|---|
| Python | ≥ 3.10 | Your UEA machine runs 3.10. |
| `numpy` | ≥ 1.23 | Core arrays. |
| `pycoare` | ≥ 0.4.3, < 0.5 | It is pre-1.0, so a minor release may change its interface. The upper bound stops a new release breaking users; the scheduled CI job tells us when to raise it. |
| `scipy` | ≥ 1.9 | Banded solver for the experimental column model only. |
| `xarray`, `pandas`, `netCDF4`, `matplotlib` | optional extra `[verify]` | Only the verification tools need them. |

`pycoare` exposes what we need through its public objects (final z0, 10 m neutral wind, gust factor, u*), so no upstream change is required. Its `util` module (stability functions, humidity) is not in the package's declared public list, so we keep our own copies of those few published formulas and test them against `pycoare`.

**Licence:** MIT, the same as `pycoare`. The new package contains no NOAA or Blomquist code, because the open-water part comes from `pycoare`. That keeps the licence clean.

## 5. Tooling

| Tool | What it does | Why |
|---|---|---|
| `pyproject.toml` with `hatchling` | One file that defines the package and how to build it | The current standard; replaces `setup.py`. |
| `uv` (pip still works) | Creates environments and installs quickly; writes a lock file | Reproducible installs for CI and for Claude sessions. On the UEA machine `pip install -e .` is enough. |
| `ruff` | Checks style and common errors, and formats the code | One fast tool in place of several. |
| `pytest`, with `hypothesis` | Runs the tests; `hypothesis` generates many random inputs for property tests | The standard test runner. |
| `pre-commit` | Runs `ruff` and file checks before each commit | Stops mistakes early. It will also block PDFs, netCDF files and large files, which protects the restricted data and the papers. |
| GitHub Actions | Runs the tests on every push, on Python 3.10 to 3.13 | Nothing merges with failing tests. |
| Scheduled CI job (weekly) | Runs the tests against the newest `pycoare` | We learn of upstream changes before users do. |
| MkDocs with `mkdocstrings` | Builds the documentation from Markdown and docstrings | Simpler than Sphinx and Markdown-native, which suits your methods docs. Built in CI from the start; published at the first public release. |
| `hatch-vcs` | Takes the version number from the git tag | One source of truth for the version. |

**Workflow:** `main` is protected. Each backlog task is done on its own branch and merged through a pull request, so every change has a description, a review and a test run.

## 6. Test strategy

| Kind | What it checks |
|---|---|
| Open-water regression | Our open-water tile equals plain `pycoare` exactly. |
| Legacy regression | `legacy.coare35vn_ice` equals stored v2 outputs on synthetic inputs, before any audit fix is applied. |
| Published values | A87 ratios at set R*; form drag against Elvidge (2016) Eq. 11 and Table 1; Buck ice saturation; C_DN10i from z0; SHEBA functions against Grachev (2007) Eq. 9. |
| Limits | SIC 0 gives open water; SIC 1 gives ice with zero form drag; neutral case; three tiles with zero thin-ice fraction equals two tiles. |
| Energy balance | SEB closes; one root; monotonic in surface temperature. |
| Column model | Analytic slab solution; energy conservation; conduction output when the surface is warmer than the base. |
| Properties (random inputs) | Finite outputs for valid inputs; NaN propagation; scalar and array inputs agree; height consistency. |
| Benchmark against observations | Metrics per campaign, compared with the previous version. Runs on your machine or where data are available, not in public CI. |

The Phase A audit scripts become the first of these tests.

## 7. Versions and releases

- **Semantic versioning** (major.minor.patch). 0.x while developing; 1.0.0 is the version described in the paper.
- **Planned steps:** 0.2.0 is v2 ported, behaviour-identical. 0.3.0 applies the audit fixes, one per pull request. 0.4.0 adds tiles and stability options. 0.5.0 adds the verification toolkit.
- **`CHANGELOG.md`** records every change that alters results, with its size.
- **Zenodo DOI per release.** GMD needs the exact archived version. The automatic GitHub–Zenodo link works for public repositories, so either the repository becomes public at submission, or we upload the release archive to Zenodo by hand. This is your decision at that time.

## 8. Attribution

- **`CITATION.cff`:** authors are you and, if he agrees, Ian. Claude is named in a contributors note, not as an author.
- **AI-assistance statement** in the README and docs: which parts were drafted by Claude models, and that all physics and results were reviewed by the authors.
- **Commit trailers:** `Co-Authored-By` on commits made by Claude.

## 9. Data handling

| Where | What |
|---|---|
| Git repository | Code, docs, small synthetic test files. No campaign data, no PDFs. |
| UEA machine | All raw and restricted data: ROVER, ERA5, satellite, campaign files. The package finds them through one environment variable, `COARE_SEAICE_DATA`. |
| CEDA and other archives | Campaign originals stay where they are published; readers document how to fetch them. |
| Zenodo | Release archives; later, a verification bundle of derived public data if licences allow. |

**Claude Code on the web:** the cloud environment needs only PyPI and GitHub, which are allowed by default. Its setup script is `pip install -e ".[dev,verify]"`. Data downloads and campaign runs happen on the UEA machine, not in the cloud environment.

Your conventions carry over and go into `CLAUDE.md`: logs to `outputs/logs/<name>.log`; notebooks from `build_*.py`; `REPROCESS = True` in script configs; methods written into the docs; figure rules; flag, don't adjust.

## 10. Migration from v2

1. **Your notebooks keep working.** `from coare_seaice.legacy import coare35vn_ice` has the v2 arguments and returns the same 56 columns.
2. **Behaviour-identical first.** The port is proven equal to v2 before anything changes.
3. **Audit fixes arrive one at a time**, each with a changelog line saying how much results move.
4. **New work uses `ice_fluxes`.** The legacy wrapper stays until your notebooks have moved.

## Decisions I need from you

1. **Name:** `coare-seaice` / `coare_seaice`?
2. **Licence:** MIT?
3. **Authors in `CITATION.cff`:** you and Ian?
4. **Docs tool:** MkDocs, as proposed?
5. **GitHub:** the account name and an empty private repository for Phase E.
