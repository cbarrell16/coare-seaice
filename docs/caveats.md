# Caveats

These apply to every result from this package and are repeated wherever results are shown.

## Inputs dominate the uncertainty

- **Sea-ice concentration.** Area-mean fluxes scale with the open-water fraction.
  Satellite concentration cannot resolve ship-scale leads and may be biased high in the
  pack.
- **What the non-ice fraction is.** In cold conditions the dark fraction of the surface is
  often thin ice at close to the air temperature, not open water at the freezing point.
  Treating it as open water overestimates the heat fluxes across the marginal ice zone.
- **Ice roughness** is regional and rarely measured. Published medians of the neutral
  drag coefficient over ice span roughly 1 to 4 x 10^-3 (Srivastava et al. 2022, Table 2).
- **Thickness and snow depth** control the energy-balance surface temperature and are
  usually unknown along a ship or aircraft track.

## ERA5 downwelling longwave

ERA5 is thought to overestimate cloud in this region, so its downwelling longwave is
probably too high. This is not corrected. Where ERA5 longwave drives the surface energy
balance, the resulting ice surface temperature is the warm end of a range. Campaigns with
measured longwave are the primary tests of the method.

## Stable conditions

- The stability functions are validated up to a bulk Richardson number of about 0.2
  (Grachev et al. 2007). Results beyond that are flagged.
- Ship measurement heights of 20 to 40 m are often above the surface layer in stable
  conditions.

## Scope

- **Cold season only.** No melt ponds or wet snow.
- **SHEBA is not independent evidence** for the formulas that were fitted to SHEBA data
  (Grachev et al. 2007; Andreas et al. 2010).
