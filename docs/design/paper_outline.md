# COARE-ice: Phase D paper outline

*13 Nov 2026. Claude (Fable) for Chris Barrell, to share with Ian Renfrew. Checkpoint D of the project brief.*

The outline is journal-agnostic. It follows the structure Geoscientific Model Development (GMD) asks for, because that is the most demanding of the likely venues; a paper that meets it can be cut down for others.

**One constraint to know first.** The Copernicus AI policy (read today on a sister journal's site) says generative AI must not be used to produce a manuscript's text or interpretations, and that generative use, such as in writing code, must be declared in the Methods or Acknowledgements. So for GMD or any Copernicus journal the paper text has to be written by you and Ian. This outline gives structure and content only.

## Working title

**coare-seaice v1.0: a tile-based extension of the COARE bulk flux algorithm to sea ice and the marginal ice zone**

GMD requires the model name and version number in the title.

## What the paper claims

1. A documented, tested, open tool that gives area-mean turbulent fluxes over mixed ice and water from standard surface-layer measurements.
2. Momentum: form drag with a regional ice roughness reproduces observed stress across the MIZ.
3. Heat and moisture: the surface temperature of the non-thick-ice fraction controls the MIZ fluxes. Treating it as open water at the freezing point overestimates them.
4. An energy-balance ice surface temperature improves on the air-temperature assumption, within stated limits.
5. Error bars: the uncertainty from ice concentration, longwave, thickness and roughness is quantified and delivered with the fluxes.

Claims 3 and 4 are provisional. They rest on the feasibility checks so far and must be confirmed by the Phase C experiments.

## Sections

### Abstract
Purpose; what the tool adds to COARE; the evaluation datasets; the main results in numbers; availability.

### 1. Introduction
- Why turbulent fluxes over the MIZ and pack matter, and why ship and aircraft bulk estimates need an ice-aware algorithm.
- COARE as the open-ocean standard (Fairall et al. 2003; Edson et al. 2013).
- The physics available: form drag (Lüpkes et al. 2012; Elvidge et al. 2016; Renfrew et al. 2019); scalar exchange (Andreas 1987; Elvidge et al. 2021, 2023; Srivastava et al. 2024); stable boundary layers over ice (Grachev et al. 2007; Andreas et al. 2010).
- The gap: these exist in models and papers, not as a tool for observational use with uncertainty.
- Aims and structure.

### 2. Model description
- **2.1 Overview and tile framework.** Tiles, fractions, the two blend modes. The simple two-tile reference.
- **2.2 Open water.** COARE 3.5 and 3.6 through `pycoare`; what is taken from it.
- **2.3 Surface-layer solution over ice.** Equations, iteration, stability functions (COARE and SHEBA), light-wind terms, validity limits.
- **2.4 Roughness.** Momentum roughness options; Andreas (1987) scalar roughness.
- **2.5 Form drag.** Lüpkes et al. (2012) with the E2016B parameters.
- **2.6 Ice surface temperature.** Air-temperature assumption; prescribed; steady surface energy balance with conduction through ice and snow. The time-dependent column as an experimental option.
- **2.7 Uncertainty estimation.** Perturbed-input ensembles.
- **2.8 Quality flags.** What is flagged and why.

Enough detail for someone to re-implement it, as GMD asks.

### 3. Implementation and verification
- Software design, interface, dependencies, licence, hardware and software requirements.
- **Verification** in GMD's sense (the equations are solved correctly): tests against published values, limiting cases, the analytic conduction solution, energy closure, and equality with `pycoare` over open water.
- Migration from the earlier `coare35vn_ice` and the errors corrected.

### 4. Evaluation data and methods
- Table of campaigns: platform, period, region, ice regime, variables, measurement heights, access.
- The common data format, quality control, matching and averaging.
- Metrics and stratification.
- Which parameters were fixed in advance and which were fitted, and on what data.

### 5. Evaluation
- **5.1 Momentum across the MIZ.** IGP and ACCACIA; roughness options.
- **5.2 Heat and moisture across the MIZ.** Two tiles, three tiles, one aggregate temperature; anchoring of the open-water coefficients.
- **5.3 Ice surface temperature.** Methods against radiometer and tower measurements; thickness and snow sensitivity.
- **5.4 Stable conditions over pack ice.** MOSAiC as the independent test; SHEBA as the check on the implementation of SHEBA-derived formulas.
- **5.5 Sensitivity and uncertainty.** Ice concentration, longwave source, thickness, snow, roughness.

### 6. Application: ROVER 2025
Fluxes along the Kronprins Haakon track with uncertainty bands; by ice-fraction class; comparison with ERA5. If the ROVER data remain restricted, this section can be replaced by the IGP ship record, which is public.

### 7. Discussion
- Scope: cold season only; no melt ponds or wet snow.
- What limits accuracy: thin-ice fraction, regional roughness, longwave input, measurement height in stable conditions.
- The ERA5 cloud and longwave caveat.
- Roadmap.

### 8. Conclusions

### Code and data availability, author contributions, competing interests, acknowledgements, references

## Figures and tables

Figure numbers match the Phase C design.

| Figure | Content | Section |
|---|---|---|
| 1 | Schematic of the tile framework and the two blend modes | 2.1 |
| 2 | Map of campaigns and tracks | 4 |
| 3 | Drag coefficient against ice fraction, observations and scheme | 5.1 |
| 4 | Heat and moisture coefficients against ice fraction and roughness Reynolds number | 5.2 |
| 5 | Modelled against observed stress and heat fluxes by campaign, reference and best configuration | 5.1, 5.2 |
| 6 | Flux error against ice fraction for the three surface descriptions | 5.2 |
| 7 | Modelled against measured surface temperature by method | 5.3 |
| 8 | Flux error against bulk Richardson number by stability option | 5.4 |
| 9 | Sensitivity summary | 5.5 |
| 10 | ROVER time series with uncertainty bands | 6 |

| Table | Content |
|---|---|
| 1 | Symbols, parameters and default values, with sources |
| 2 | Campaign datasets |
| 3 | Experiments |
| 4 | Summary metrics by campaign and configuration |
| 5 | Quality flags |

## Draft statements

These are drafts for you to rewrite in your own words.

**Code availability.** The exact version of coare-seaice described here (v1.0.0) is archived on Zenodo (DOI to be added) under the MIT licence. Development takes place at https://github.com/cbarrell16/coare-seaice. The package depends on pycoare (version to be added), which provides the COARE 3.5 and 3.6 open-water algorithms. The scripts that produce the figures are included in the archive.

**Data availability.** The IGP and ACCACIA aircraft and ship observations are available from CEDA (identifiers to be added). MOSAiC data are available from the Arctic Data Center and SHEBA data from the NCAR Earth Observing Laboratory (identifiers to be added). ERA5 is available from the Copernicus Climate Change Service. The ROVER 2025 ship observations are currently restricted to the ROVER working group and will be made available at (repository) on (date); until then they are available from the authors on request.

**AI use** (for the Methods or Acknowledgements). Parts of the coare-seaice software, its tests and its documentation were drafted with the assistance of Claude large language models (Anthropic). An independent review of the code against the source literature was also carried out with these tools. All equations, parameter values and results were checked against the cited papers and verified by the authors, who take full responsibility for the content. The models are credited in the software repository and are not authors of this paper.

## GMD requirements to meet

| Requirement | Status |
|---|---|
| Model name and version in the title | In the working title |
| Description detailed enough to re-implement | Section 2 |
| User manual | The package documentation |
| Source code on a persistent public archive with an identifier, exact version | Zenodo at release 1.0.0; planned |
| Code availability section | Drafted above |
| Licence, web address, hardware and software requirements in the text | Section 3 |
| Verification kept distinct from evaluation | Sections 3 and 5 |
| Scope and limitations discussed | Section 7 |
| Generative AI use declared; not used for text or interpretation | Statement drafted; text to be written by the authors |

The manuscript-type requirements are from the GMD site, read today. I could not open GMD's separate code and data policy page from this session; from memory it requires a persistent archive with a DOI at submission and does not accept a GitHub link alone. Please check it before submission.

**Manuscript type.** "Model description paper" fits best. "Development and technical paper" is the alternative if the editors see it as a parameterisation added to COARE.

## If you choose a lighter venue

Not checked today; these are suggestions from general knowledge for you and Ian to weigh.

- **A short software paper** (for example the Journal of Open Source Software): a one- to two-page summary, with review focused on the repository, tests and documentation. It gets the tool citable quickly and leaves the science for a later paper.
- **A measurement or boundary-layer journal** (for example the Journal of Atmospheric and Oceanic Technology or Boundary-Layer Meteorology, where Srivastava et al. 2024 and Andreas 1987 appeared): the method and its evaluation, with lighter software requirements.

A reduced paper would keep Sections 2, 4 and 5.1–5.3 with the IGP and ACCACIA aircraft data and one tower, and move the uncertainty work and ROVER application to a follow-up.

## Questions for you and Ian

1. Is the central message right: a tool paper whose main scientific finding is that MIZ heat fluxes depend on the temperature of the non-thick-ice surface?
2. Full GMD-style paper, or a short software paper first with the science to follow?
3. Should ROVER be the application, given the data restriction, or the public IGP ship record?
