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
