# Changelog

All changes that alter results are recorded here, with their size.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the
project uses [semantic versioning](https://semver.org/).

## [Unreleased]

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

## History before this repository

- **v2 (10 Nov 2026), `coare35vn_ice.py` 2.0.0.** Four errors in v1 corrected: humidity
  over ice 1000 times too small; measured relative humidity treated as relative to ice;
  latent heat of sublimation set to the vaporisation value; form drag ten times too small
  (c_e squared instead of c_e). Surface energy balance and conduction-column ice surface
  temperature methods added.
- **v1 (May 2026), `coare35vn_ice.py`.** First sea-ice extension of the NOAA Python
  COARE 3.5 code.
