"""Thermal properties of sea ice and snow.

Conductivity and heat capacity of sea ice follow Bitz & Lipscomb (1999):

    k(S, T) = k0 + beta S / T        which they cite to Untersteiner (1964)
    c(S, T) = c0 + L0 mu S / T^2     their Eq. 1 with gamma = L0 mu (their Eq. 2)

with S in ppt and T in deg C. Untersteiner (1961) gives the heat-capacity form with a
constant 5 % smaller (17.2 against 18.0 kJ K kg-1 ppt-1); the Bitz & Lipscomb value is
used here.
"""

K_ICE_FRESH = 2.03  # W m-1 K-1, conductivity of fresh ice, k0 (Untersteiner 1964, via
#                     Bitz & Lipscomb 1999, who give 2.034; the v2 value 2.03 is kept)
BETA_K = 0.117  # W m-1 ppt-1, salinity term beta (Untersteiner 1964, via Bitz & Lipscomb 1999)
K_SNOW = 0.31  # W m-1 K-1, snow conductivity (Bitz & Lipscomb 1999; CICE default 0.30)
RHO_ICE = 917.0  # kg m-3
RHO_SNOW = 330.0  # kg m-3 (Bitz & Lipscomb 1999)
C_ICE0 = 2110.0  # J kg-1 K-1, heat capacity of fresh ice, c0 (Bitz & Lipscomb 1999, Eq. 1)
C_ICE0_V2 = 2106.0  # the v2 value, kept only for the regression tests against v2
L0_MU = 3.34e5 * 0.054  # J K kg-1 ppt-1, brine term L0 mu (Bitz & Lipscomb 1999, Eq. 2)
C_SNOW = 2106.0  # J kg-1 K-1, heat capacity of snow (v2 value; no source given)
S_ICE_DEFAULT = 5.0  # ppt, bulk salinity of first-year ice

ALBEDO_SNOW = 0.85  # dry snow, winter
ALBEDO_ICE = 0.65  # bare cold ice
