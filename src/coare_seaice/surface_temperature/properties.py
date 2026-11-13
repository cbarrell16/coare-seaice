"""Thermal properties of sea ice and snow (v2 ice_surface_temperature.py)."""

# AUDIT F9: citation to be corrected to Untersteiner (1964) (task F9).
K_ICE_FRESH = 2.03  # W m-1 K-1, conductivity of fresh ice (Untersteiner 1961)
BETA_K = 0.117  # W m-1 ppt-1, salinity term: k = 2.03 + 0.117 S / T (Untersteiner 1961)
K_SNOW = 0.31  # W m-1 K-1, snow conductivity (Maykut & Untersteiner 1971; CICE default 0.30)
RHO_ICE = 917.0  # kg m-3
RHO_SNOW = 330.0  # kg m-3
# AUDIT F9: fresh-ice heat capacity to be changed to 2110 J kg-1 K-1 (Bitz & Lipscomb 1999)
# in task F9.
C_ICE0 = 2106.0  # J kg-1 K-1, fresh ice heat capacity
L0_MU = 3.34e5 * 0.054  # J K kg-1 ppt-1, brine term: c = c0 + L0 mu S / T^2 (Bitz & Lipscomb 1999)
C_SNOW = 2106.0  # J kg-1 K-1
S_ICE_DEFAULT = 5.0  # ppt, bulk salinity of first-year ice

ALBEDO_SNOW = 0.85  # dry snow, winter
ALBEDO_ICE = 0.65  # bare cold ice
