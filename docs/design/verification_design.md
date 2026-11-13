# COARE-ice: Phase C verification study design

*13 Nov 2026. Claude (Fable) for Chris Barrell. Checkpoint C of the project brief. A design for your approval; the numbers quoted from the IGP database are a feasibility check, not results.*

## Aims

The study has to answer six questions, in this order of importance for the paper:

1. **Momentum:** does the form-drag scheme with a chosen ice roughness reproduce observed stress across the MIZ and pack?
2. **Heat and moisture across the MIZ:** which surface description works best: two tiles, three tiles with thin ice, or one aggregate surface temperature?
3. **Ice surface temperature:** how well do the methods (`air_freezing`, `seb`, satellite or reanalysis prescribed) reproduce measured surface temperature, and does a better temperature give better fluxes?
4. **Stable conditions over pack ice:** which stability functions and light-wind terms work, and where should results be flagged as unreliable?
5. **Uncertainty:** how large are the flux error bars from ice concentration, longwave, thickness, snow and roughness?
6. **Ship geometry:** does the method hold at ship measurement heights (20–40 m)?

## Common method

- **One data format.** Every campaign is read into the same table: one row per sample (aircraft run, tower hour, ship 10-min window), with fixed variable names, measurement heights per row, and the campaign's own quality flags kept.
- **Quality control.** Use each campaign's published flags. For IGP that is `z0_flag`, `z0T_flag`, `z0q_flag` ≤ 1 and `icefractionalb`, as in the readme. Flag, never adjust.
- **Reference configuration.** Two tiles, ice at min(T_air, −1.8 °C), z0_ice = 5 × 10⁻⁴ m, COARE stability functions. Every experiment is reported as a change from this.
- **Metrics.** Bias, median error, RMSE, correlation and regression slope, each with a bootstrap confidence interval, for stress, sensible heat, latent heat and surface temperature. Means and percentiles are always reported with medians (audit finding 2).
- **Stratification.** Ice fraction (the Elvidge bins: 0, 0–0.25, 0.25–0.5, 0.5–0.75, 0.75–1, 1); stability (bulk Richardson number); wind speed; downwelling longwave or cloud class.
- **No tuning on the test data.** Parameters are fixed from the literature before each comparison. Where a parameter is fitted (for example a regional roughness), it is fitted on one campaign and tested on another, and the paper says so.
- **Methods are written into the docs** for every campaign: matching, averaging, binning, height handling.

## Experiments

Each is run on every dataset where the inputs exist.

| # | Experiment | Options | Main question |
|---|---|---|---|
| E0 | Reference | as above | baseline |
| E1 | Ice surface temperature | `air_freezing`; `seb`; prescribed measured; prescribed satellite; prescribed ERA5 ice tile | 3 |
| E2 | Thickness and snow | ice 0.3, 0.5, 1, 2 m; snow 0, 0.05, 0.2 m; snow conductivity 0.14 and 0.31 W m⁻¹ K⁻¹ | 3, 5 |
| E3 | Ice roughness | 5 × 10⁻⁴ m; Renfrew et al. (2019) "smooth" and "rough" settings; campaign-median observed; Andreas et al. (2010) Eq. 5.1 | 1 |
| E4 | Stability and light wind | COARE; SHEBA (Grachev et al. 2007); with and without Andreas et al. (2010) gustiness and windless terms | 4 |
| E5 | Blend mode | separate tile temperatures; one aggregate temperature | 2 |
| E6 | Thin-ice tile | thin fraction 0 to all of the non-thick fraction; thin-ice thickness 0.05–0.3 m | 2 |
| E7 | Ice concentration | campaign-measured; AMSR2; OSI SAF; ± 0.05 | 5 |
| E8 | Longwave source | measured; ERA5; König-Langlo | 3, 5 |
| E9 | Scalar roughness | Andreas (1987); fixed ratio 0.2; MetUM constants from Renfrew et al. (2019) | 2 |
| E10 | Open-water base | COARE 3.5; 3.6; water coefficients anchored to the campaign's open-sea median | 2 |

**E10 matters more than I expected.** In the IGP database the observed open-sea scalar coefficients are well below COARE's (median C_HN 0.68 × 10⁻³, C_EN 0.36 × 10⁻³). Elvidge et al. (2021) anchor the water coefficients to the observed sea-bin medians for this reason. Without that step, an open-water mismatch is carried into every MIZ bin and looks like an ice error.

## Datasets

### IGP aircraft, Feb–Mar 2018 (in hand) — essential

- **Tests:** questions 1, 2, 3, 5. Rough ice, cold-air outbreaks, MIZ to pack.
- **Inputs available:** wind, temperature, humidity, pressure and height per run; radiometer surface temperature; measured longwave and shortwave; albedo; ice fraction; OSTIA SST; eddy-covariance stress, sensible and latent heat; quality flags. 339 long runs; 190, 239 and 216 pass the momentum, heat and moisture flags.
- **Matching:** none needed; one row per 9 km run. Heights are per run (median 54 m, so the low-level subset below 100 m is used for bulk comparison).
- **Missing:** ice thickness and snow depth (use E2 as a range); short-run surface variables (to come later).
- **Feasibility check done today:** the bin medians and the stress and heat-flux errors in the A2 doc. The mean bias in total turbulent heat flux across the MIZ is +19 W m⁻² with one aggregate temperature, +22 with the reference, +34 with two tiles and `seb`. Elvidge et al. (2021) report 13 W m⁻² for Blended A87 on these data, with anchored water coefficients, which I have not yet applied.

### ACCACIA aircraft, March 2013 (in hand since 13 Nov) — essential, and the independent test

- **Tests:** questions 1 and 2 over smoother ice (Elvidge et al. 2021 give a median z0_ice of 4.3 × 10⁻⁴ m, against 10⁻² m for IGP).
- **Inputs available:** the same format and variables as the IGP database. 201 long runs at 45–100 m; 170, 178 and 174 pass the momentum, heat and moisture flags; 172 have some ice.
- **Use:** any regional parameter fitted on IGP is tested here, and the reverse.
- **Feasibility check done:** it confirms the IGP result with a larger MIZ sample. Mean bias in total turbulent heat flux across the MIZ (123 runs): +1 W m⁻² with one aggregate surface temperature, +30 with the two-tile reference, +39 with two tiles and `seb`. The drag coefficient peaks at 75–100 % ice (2.1 × 10⁻³, against 1.3 × 10⁻³ over full ice), which the anchored E2016B curve underestimates there and overestimates at lower ice fractions.
- **Unexplained:** the open-sea drag coefficient is very low in the flagged runs (median 0.48 × 10⁻³). To check against Elvidge et al. (2016) before use.

### MOSAiC tower, 2019–20 — essential

- **Tests:** questions 3 and 4 over a full winter on pack ice, with standard heights (2–10 m).
- **Inputs expected** (from memory of Cox et al. 2023; to confirm on download): multi-level wind, temperature and humidity; eddy-covariance fluxes; skin temperature; all four radiation components; snow depth. Ice thickness from other MOSAiC datasets.
- **Matching:** the tower's own averaging period; use one level as input and compare with the measured fluxes; repeat with another level as a height-consistency test.
- **Why essential:** independent of the formulas being tested.

### SHEBA tower, 1997–98 — essential (your request)

- **Tests:** questions 3 and 4. It is the foundation of the Grachev et al. (2007) functions and the Andreas et al. (2010) algorithm, so agreement here is a check on our implementation, not independent evidence. The paper must say that.
- **Inputs expected** (to confirm): five-level hourly tower data, radiation, surface temperature, snow depth.
- **Target to reproduce:** the seasonal fluxes in Persson et al. (2002), including the small downward winter sensible heat flux.

### Ship campaigns in ice (AO2016, AO2018; ACSE when obtained) — desirable

- **Tests:** question 6, the ship-in-ice geometry closest to ROVER, mostly summer and autumn.
- **Target:** the Blended A87 results in Srivastava et al. (2024) and the drag values in Srivastava et al. (2022) Table 2.
- **Access:** AO2016 and AO2018 on the Bolin Centre database; ACSE to request.

### IGP ship (NRV Alliance), 2018 (in hand) — desirable

- **Tests:** question 6 and the open-water path. The ship entered the MIZ on several occasions (Renfrew et al. 2021). I earlier assumed it had not; the A2 doc is corrected.
- **Limit:** bulk measurements only as far as I can see from the inventory, so it tests consistency with the aircraft and ERA5, not against ship eddy covariance.

### N-ICE2015 — desirable

- **Tests:** questions 3 and 4 over young, thin winter ice north of Svalbard, the regime closest to ROVER.

### ROVER 2025 — application, not verification

- **Use:** the worked example, with uncertainty bands from E2, E3, E7 and E8.
- **Limits:** restricted data; no surface temperature; many ice windows carry a wind-quality flag; drone and lidar data will come later and can be added as verification then.

### Satellite ice surface temperature — optional

- **Product:** Copernicus Marine Level 4 reprocessed, daily, 0.05°, 1982 to June 2026. The Level 3 product (clear-sky pixels only) is better for strict comparison if it covers the dates; to check.
- **Use:** compare with radiometer temperature on IGP flight days first, to learn the product's error; only then use it as a reference for ROVER.
- **Caveat:** infrared retrievals are clear-sky only, so the sample is biased to cold, clear conditions, and Level 4 fills cloudy areas by interpolation.

## Comparison with published results

These are the paper's anchor points. Each becomes a figure or a table row.

| Paper | Result to compare against |
|---|---|
| Elvidge et al. (2016) | Drag coefficient against ice fraction; the E2016B curve. |
| Elvidge et al. (2021) | Scalar coefficients against ice fraction and roughness Reynolds number; heat-flux bias across the MIZ of 13 W m⁻² (Blended A87) against 48 (MetUM) and 80 (IFS). |
| Renfrew et al. (2019) | The MetUM "smooth", "rough" and form-drag settings as named roughness options. |
| Renfrew et al. (2021) | ERA5 against IGP observations over open water and the MIZ; larger errors over the MIZ from an over-smooth ice field. |
| Srivastava et al. (2022) | Table 2 of ice drag coefficients by campaign. |
| Srivastava et al. (2024) | Blended A87 against ship eddy covariance. |
| Andreas et al. (2010); Grachev et al. (2007); Persson et al. (2002) | SHEBA bulk algorithm, stability functions and seasonal fluxes. |

**First check on the drag curve (IGP, flagged runs):** observed median C_DN10 rises from 1.06 × 10⁻³ over open sea to 3.35 × 10⁻³ over full ice, with no peak in the MIZ. With the ice end anchored to the observed value the E2016B curve fits at 25–50 % and 75–100 % ice but is high at 0–25 % (1.5 against 0.9 × 10⁻³, 34 runs) and at 50–75 % (3.2 against 2.0 × 10⁻³, 14 runs). With the v2 default roughness it is low at high ice fractions (1.9 against 3.1 × 10⁻³).

## Figures

1. Map of all campaigns and tracks.
2. Drag coefficient against ice fraction: observations by campaign, with the scheme for each roughness option.
3. Heat and moisture coefficients against ice fraction and roughness Reynolds number.
4. Modelled against observed stress, sensible and latent heat, by campaign: reference and best configuration side by side.
5. Flux error against ice fraction for the three surface descriptions (E5, E6).
6. Surface temperature: modelled against measured, by method (E1), for IGP, MOSAiC and SHEBA.
7. Stable regime: flux error against bulk Richardson number for each stability option (E4), MOSAiC and SHEBA.
8. Sensitivity summary: change in flux for each input perturbation (E2, E3, E7, E8).
9. ROVER application: time series with uncertainty bands, and by ice-fraction class.
10. Optional: satellite against radiometer surface temperature.

Your figure rules apply: white backgrounds, legends inside panels, ice-concentration colours dark blue at 0 to white at 1, no explanatory text on figures.

## What makes it publishable

**Needed:** IGP and ACCACIA aircraft (E0–E6, E9, E10); one independent tower (MOSAiC: E1, E2, E4, E8); SHEBA as the implementation check; the ROVER application with uncertainty; the comparisons with Elvidge (2016, 2021) and Srivastava (2022).

**Strengthens it:** the ship campaigns; N-ICE2015.

**Optional:** satellite surface temperature; ROVER drone data.

## Caveats that stay in the paper

- **ERA5 longwave** is probably too high because ERA5 overestimates cloud. Results with ERA5 longwave are the warm end of a range. Campaigns with measured longwave (IGP, MOSAiC, SHEBA, N-ICE2015) are the primary tests for that reason.
- **Aircraft fluxes** are measured at about 30–50 m with a sampling error on each 9 km run, and the open-sea scalar coefficients are below COARE's.
- **SHEBA is not independent** of the SHEBA-derived formulas.
- **Ice roughness is regional** and rarely measured; it is an input with an uncertainty range, not a constant.
- **Thin-ice fraction** has no routine measurement; for ships it is a sensitivity parameter.

## Order of work, and where it runs

| Step | Work | Where |
|---|---|---|
| 1 | IGP aircraft reader and experiments E0–E6, E9, E10 | Cloud session; data file in hand (kept out of the repository) |
| 2 | Download MOSAiC and SHEBA; confirm variables and licences | UEA machine (you) |
| 3 | Tower readers; E1, E2, E4, E8 | UEA machine for runs; code written in cloud sessions against small synthetic samples |
| 4 | ROVER application with uncertainty | UEA machine |
| 5 | ACCACIA, ship campaigns, N-ICE2015 | As data arrive |
| 6 | Satellite surface temperature | UEA machine |

## Decisions (13 Nov)

- Design approved by Chris, including SHEBA as an implementation check and MOSAiC as the independent test.
- ACCACIA run database supplied.
- IGP data may go in the public repository; it is already public on CEDA. ACCACIA: to confirm the same before committing any of it.
