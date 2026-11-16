# Installing coare-seaice and reprocessing ship fluxes

## Install

Python 3.10 or later. The package pulls in `numpy`, `scipy` and `pycoare` from PyPI.
Activate the environment you use for analysis first (for a conda environment:
`conda activate <name>`).

**From GitHub (after the port pull request is merged):**

```bash
pip install "coare-seaice @ git+https://github.com/cbarrell16/coare-seaice@main"
```

**From a wheel file** (a `.whl` file you were sent):

```bash
pip install coare_seaice-<version>-py3-none-any.whl
```

**For development** (you want to edit the code):

```bash
git clone https://github.com/cbarrell16/coare-seaice
cd coare-seaice
pip install -e ".[dev,verify]"
pytest
```

Check it worked:

```bash
python -c "import coare_seaice, pycoare; print(coare_seaice.__version__)"
```

To upgrade later, run the same install command again (add `--upgrade` for a wheel).

## Switching existing code from `coare35vn_ice.py`

Replace

```python
from coare35vn_ice import coare35vn_ice
```

with

```python
from coare_seaice.legacy import coare35vn_ice
```

The call and the 56 output columns are the same as v2.

**Three things to check in existing notebooks:**

1. **Set `t_ice_method` explicitly.** The v2 default is `'seb_steady'`, not the v1
   behaviour. For the bug fixes only, pass `t_ice_method='air_freezing'`.
2. **Missing sea-ice concentration.** A NaN concentration gives NaN blended fluxes. If
   older code filled missing values with 0 (open water), decide which you want and
   record it in the methods.
3. **Results differ from v1** because of the four v1 errors corrected in v2 (see
   `CHANGELOG.md`). The largest changes are a higher blended momentum flux in the MIZ
   (the form-drag correction) and a higher latent heat flux over ice (the humidity
   correction).

## Reprocessing the ROVER 2025 ship fluxes (task P5)

`scripts/p5_reprocess_rover_legacy.py` reprocesses the ship fluxes with the v2 bug fixes
only (`air_freezing`) and checks the new package against the original v2 code.

1. Put these in one folder:
   - `ROVER2025_KronprinsHaakon_surface_met_fluxes_10min.nc`
   - `era5_along_track_10min.nc` (for longwave and missing shortwave)
   - a folder `coare_ice_v2/` holding the v2 code (`coare35vn_ice.py`,
     `ice_surface_temperature.py`, `meteo.py`, `util.py`), for the check.
2. Run, from wherever you want the `outputs/` folder:

   ```bash
   export COARE_SEAICE_DATA=/path/to/that/folder
   python /path/to/p5_reprocess_rover_legacy.py
   ```

   Paths and heights can also be changed in the config block at the top of the script.
3. Read `outputs/logs/p5_reprocess_rover_legacy.log`:
   - **Ice-only columns (38–46, 50–54):** should be identical to v2 (difference 0).
   - **Open-water columns:** differ by up to about 0.1 W m⁻² (heat) and 1e-4 N m⁻²
     (stress), because `pycoare` and the NOAA code differ slightly in gravity and a few
     constants.
   - **NaN mismatches:** should be 0 in every column.
4. The reprocessed fluxes are in
   `outputs/data/ROVER2025_coare_seaice_air_freezing_10min.nc`, with the settings in its
   attributes.

**Method, as written into the output attributes:** ship wind, temperature, humidity,
pressure, SST and AMSR2 sea-ice concentration; heights 36 m (wind) and 22 m (temperature
and humidity); shortwave from the ship where available, else ERA5; longwave from ERA5
along-track; no gap filling; z0_ice = 5 × 10⁻⁴ m; T_ice = min(T_air, −1.8 °C). See
[caveats](../caveats.md), in particular the ERA5 longwave caveat.

## A `pycoare` pitfall

`pycoare` 0.4.3 divides a float `rh` array by 100 **in place**. A second call with the
same array then uses humidities 100 times too small, with no warning. `coare-seaice`
protects itself, but if you call `pycoare` directly, pass a copy: `rh=rh.copy()`.
