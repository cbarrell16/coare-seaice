# Changelog

All changes that alter results are recorded here, with their size.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the
project uses [semantic versioning](https://semver.org/).

## [0.3.0] - unreleased (branch `release/0.3.0`)

Audit fixes F1-F10. "Reference fixture" below means the 300 synthetic polar cases in
`tests/fixtures/v2_reference` (heights 36/22/22 m); changes are new minus old. Sizes on
the IGP and ACCACIA aircraft data will be added when the reader exists (task V2).

### Added
- **F1.** The limit on the stability parameter in the ice surface-layer iteration is now
  a named, documented option (`most.SurfaceLayerOptions.zeta_max`, default 10 as in v2)
  in place of a hard-coded number. New diagnostics for the ice fraction: `zeta`, the
  bulk Richardson number `Rib` (Grachev et al. 2007, Eq. 8) and two flags:
  `zeta_limited` (the limit acted, so the result is set by it) and `qc_rib`
  (`Rib` > 0.2, outside the validated range of the stability functions).
  `legacy.coare35vn_ice(..., return_diagnostics=True)` returns them.
  **No change in results.** On the reference fixture the flags are set in 3.0 %
  (`zeta_limited`) and 5.3 % (`qc_rib`) of cases with `air_freezing`, and 14.0 % and
  12.3 % with `seb_steady`.
- **F6.** Convergence diagnostics for the ice surface-layer iteration: `iter_dH` and
  `iter_dLE` (change in the ice heat fluxes in the last of the 10 iterations, W m-2)
  and the flag `not_converged` (either above `SurfaceLayerOptions.conv_tol`, default
  0.1 W m-2). **No change in results.** On the reference fixture every case converges;
  the largest last-iteration change is 0.013 W m-2. The audit also asked for u* to be
  updated in the final pass. On inspection no update is needed: the u* of the last
  iteration is already computed from the final Obukhov length. A test now checks this.

### Changed (results change)
- **F2.** The roughness Reynolds number R* = z0 u* / nu in the Andreas (1987) scalar
  roughness now uses the temperature-dependent air viscosity of COARE
  (`thermo.air_viscosity`, 1.24e-5 m2 s-1 at -10 C) in place of the constant 1.5e-5.
  Colder air gives a larger R* and a smaller z0T and z0q, so the magnitudes of the ice
  heat fluxes fall. Reference fixture, `air_freezing`: ice sensible and latent heat
  change by -2.4 % and -2.2 % (medians where the flux exceeds 5 W m-2; range -0.1 to
  -4 %); the largest changes are 1.8 W m-2 in the ice fluxes and 1.2 to 1.5 W m-2 in
  the blended fluxes; ice stress changes by at most 0.001 N m-2 (through stability).
  `seb_steady`: ice surface temperature changes by at most 0.27 K. This is larger than
  the 1-2 % expected in the audit. The `simple` scheme (fixed ratio) is unchanged
  except for the reported R*. `most.V2_OPTIONS` restores the old value.
- **F3.** The Andreas (1987) rough-regime polynomial is no longer extrapolated beyond
  its fitted range (R* <= 1000): above that the ratios z0T/z0 and z0q/z0 are held at
  their R* = 1000 values, and the new flag `qc_rstar` is set. **No change on the
  reference fixture** (z0_ice = 5e-4 m gives R* of order 10 to 100). It matters only for
  rough ice and strong wind, for example z0_ice = 0.05 m with u* above 0.3 m s-1.

### Fixed (results change)
- **F4.** `t_ice_method="diffusion"`: the `F_cond` output (column 51) is now the
  conductive flux of the column model itself. v2 back-calculated it through an
  effective conductance, dividing by max(t_base - T_ice, 1e-6), which gave wrong and
  sometimes enormous values when the surface was warmer than the ice base (-1.8 C).
  Reference fixture, `diffusion` run: 7 of 300 cases change, all with T_ice above
  -1.8 C; the largest |F_cond| falls from 1.1e8 to 176 W m-2. Fluxes and T_ice are
  unchanged, and so are all other methods.
- **F5.** Column model (`IceColumn`, `t_ice_method="diffusion"`), three changes.
  (1) The initial "steady" profile is now the steady state of the column's own layers
  and salinity-dependent conductivity; v2 used fresh-ice conductivity, a 2.6 %
  conductance mismatch that made the profile drift after every (re)initialisation. The
  conductive flux reported at (re)initialisation is the same steady flux, with the
  column's own snow conductivity. (2) The melt flux is reset at each Picard pass (no
  case where this mattered was found in 20 000 random steps). (3) **Decision on
  thickness arrays:** one column is one slab, so `h_ice`, `h_snow` and `albedo` arrays
  are replaced by their medians everywhere in the `diffusion` method, with a warning;
  v2 used the arrays for the initial state and the medians for the column.
  Reference fixture, `diffusion` run (constant thickness): T_ice changes by at most
  0.04 K (mean -0.008 K), ice sensible heat by at most 0.35 W m-2 (mean -0.10),
  `F_cond` by at most 0.51 W m-2 (mean -0.21). Other methods are unchanged.
- **F7.** No silent fills; NaN in gives NaN out for every output that depends on the
  input (tested for every input). Three fills were found and removed or flagged.
  (1) A NaN stability parameter (missing humidity, pressure, latitude or surface
  temperature) was treated as neutral, so ice stress, R* and the scalar roughness
  stayed finite; they are now NaN. (2) A missing snow depth gave the bare-ice albedo
  0.65; it now gives NaN. (3) In `diffusion`, samples where the column is not running
  (SIC below `sic_reset`) take the steady-state T_ice as before, but are now flagged.
  New diagnostics `t_ice_filled` (T_ice does not come from the requested method,
  including the `t_ice_fallback` of `prescribed`) and `missing_input`. The default of
  `t_ice_fallback` is unchanged (`seb_steady`); pass `None` for no fill.
  **No change on the reference fixture** (all ten runs identical).

## [0.2.0] - 16 Nov 2026 (not tagged)

### Added
- Repository layout, packaging, continuous integration, conventions (`CLAUDE.md`),
  task backlog (`TASKS.md`) and design documents.
- The v2 ice physics ported into modules with no change in behaviour, and
  `coare_seaice.legacy.coare35vn_ice` with the v2 call and 56-column output (tasks P1-P4).
  Ice-only columns are identical to v2. Open-water columns now come from `pycoare` and
  differ from v2 by up to about 0.1 W m-2 in heat fluxes and 1e-4 N m-2 in stress,
  because `pycoare` and the NOAA code differ slightly in gravity and in a few constants.
  The defects found in the audit are deliberately still present and are marked
  `# AUDIT F<n>` in the code; they are fixed one at a time in tasks F1-F10.
- `scripts/p5_reprocess_rover_legacy.py`: reprocesses ship fluxes with `air_freezing` and
  checks the package against v2 (task P5).
- Guides for reviewing and merging pull requests, and for installing and reprocessing.
- Draft bug report for `pycoare` (in-place change of the caller's `rh` array).

### Changed
- Copyright holder in `LICENSE`: Chris Barrell.

## History before this repository

- **v2 (10 Nov 2026), `coare35vn_ice.py` 2.0.0.** Four errors in v1 corrected: humidity
  over ice 1000 times too small; measured relative humidity treated as relative to ice;
  latent heat of sublimation set to the vaporisation value; form drag ten times too small
  (c_e squared instead of c_e). Surface energy balance and conduction-column ice surface
  temperature methods added.
- **v1 (May 2026), `coare35vn_ice.py`.** First sea-ice extension of the NOAA Python
  COARE 3.5 code.
