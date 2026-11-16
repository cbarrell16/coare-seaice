# Test fixtures

Small files only (under 300 kB each). Synthetic inputs, stored reference outputs, and
small subsets of **public** campaign data (IGP and ACCACIA aircraft run databases, both
available from CEDA). Never restricted data.

## `v2_reference/`: stored outputs of the v2 code

Reference outputs of `coare35vn_ice.py` version 2.0.0 (Nov 2026) and its helper modules,
used to prove that the port in `src/coare_seaice/` behaves identically. Outputs only,
from synthetic inputs; no v2 code is stored.

- **How they were made:** `python tests/fixtures/make_v2_reference.py /path/to/v2`, run
  once by hand (not by CI) on 13 Nov 2026 with Python 3.13, numpy 2.5.3, scipy 1.18.1.
  The script documents every input and call. Re-running it gives identical files.
- **Inputs:** 300 synthetic polar cases, seed 20261113: wind 0.5-25 m/s, air -30 to
  +5 C, RH 50-100 %, SST -1.9 to 6 C, pressure 960-1040 hPa, latitude 60-82 N,
  SW 0-200 and LW 150-320 W m-2, SIC 0-1 with exact 0, exact 1, NaN and values below
  0.01; prescribed T_ice with NaNs; a 10-min time axis with two gaps for `diffusion`.
- **Function-level files:** `thermo_stability`, `roughness_form_drag`, `most`, `seb`,
  `column` (inputs and outputs of each v2 function on grids and the cases above).
- **Whole-function files:** `legacy_openwater_<heights>.npz` hold columns 0-37 for
  heights 36/22/22 m and 10/10/10 m; `legacy_<run>.npz` hold columns 38-55 and the
  call's keyword arguments as JSON (`"@name"` refers to an array in `inputs.npz`).
  Runs cover every `t_ice_method`, `ice_scheme='simple'`, `rh_reference='ice'`,
  thickness arrays, and a run with rain and `jcool=0` (which stores all 56 columns).
- **Size:** each file under 120 kB; about 1 MB in total.
