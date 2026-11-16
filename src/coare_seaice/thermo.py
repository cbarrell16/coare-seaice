"""Physical constants, humidity and air properties.

Units: temperature in deg C, pressure in mb (hPa), specific humidity in kg/kg unless a
name says otherwise.

The water-reference formulas below are the ones used by COARE 3.5 (and by v2
``coare35vn_ice``). They are written here from the published formulas, not copied from
the NOAA code, and the tests show they give the same numbers as v2.
"""

from __future__ import annotations

import numpy as np

# --- COARE constants (as used in v2 and in pycoare) ------------------------------------
VON = 0.4  # von Karman constant
FDG = 1.00  # turbulent Prandtl number factor
TDK = 273.16  # deg C -> K offset used by COARE
RGAS = 287.1  # gas constant for dry air, J kg-1 K-1
CPA = 1004.67  # specific heat of air, J kg-1 K-1
LAPSE_DRY = 0.0098  # dry adiabatic lapse rate used for potential temperature, K m-1

# --- Ice-surface constants (v2 ice_surface_temperature.py) -----------------------------
TK0 = 273.15  # deg C -> K offset used by the ice code
LS_ICE = 2.834e6  # latent heat of sublimation, J kg-1 (~0 C; varies < 0.2 % over -30..0 C)
T_FREEZE_SW = -1.8  # deg C, freezing point of sea water (ice base)


# --- Saturation over liquid water and air humidity -------------------------------------
def esat_water(t, p):
    """Saturation vapour pressure over liquid water (mb), Buck (1981).

    es = (1.0007 + 3.46e-6 p) * 6.1121 exp(17.502 t / (240.97 + t)),
    with t in deg C and p in mb. The first factor is Buck's enhancement factor.
    """
    t = np.asarray(t, dtype=float)
    p = np.asarray(p, dtype=float)
    return 6.1121 * np.exp(17.502 * t / (240.97 + t)) * (1.0007 + 3.46e-6 * p)


def qsat_sea(ts, p, salinity=None):
    """Saturation specific humidity at the sea surface (kg/kg).

    The saturation vapour pressure over pure water is reduced for salt: by 2 %
    (factor 0.98, as COARE 3.5 and v2) when ``salinity`` is None, otherwise by
    (1 - 0.02 S / 35), the salinity form that pycoare's COARE 3.6 uses.
    q = 0.622 e / (p - 0.378 e).
    """
    p = np.asarray(p, dtype=float)
    factor = 0.98 if salinity is None else 1.0 - 0.02 * np.asarray(salinity, dtype=float) / 35.0
    e = factor * esat_water(ts, p)
    return 622.0 * e / (p - 0.378 * e) / 1000.0


def vapour_pressure(t, p, rh):
    """Vapour pressure of the air (mb) from relative humidity (%) over liquid water."""
    rh = np.asarray(rh, dtype=float)
    return (rh / 100.0) * esat_water(t, p)


def specific_humidity(t, p, rh):
    """Air specific humidity (kg/kg) from t (deg C), p (mb) and RH (%) over liquid water.

    q = 0.62197 e / (p - 0.378 e) with e = (RH/100) es(t, p), es from :func:`esat_water`.
    """
    p = np.asarray(p, dtype=float)
    e = vapour_pressure(t, p, rh)
    return 621.97 * e / (p - 0.378 * e) / 1000.0


# --- Saturation over ice (project's own code, from v2 ice_surface_temperature.py) ------
def esat_ice(t, p):
    """Saturation vapour pressure over ice (mb), Buck (1981) 'ei2' with enhancement factor.

    es_i = (1.0003 + 4.18e-6 p) * 6.1115 exp(22.452 t / (272.55 + t)), t in deg C, p in mb.
    Same paper and form as :func:`esat_water`, so the two references are consistent.
    """
    t = np.asarray(t, dtype=float)
    p = np.asarray(p, dtype=float)
    return (1.0003 + 4.18e-6 * p) * 6.1115 * np.exp(22.452 * t / (272.55 + t))


def qsat_ice(t, p):
    """Saturation specific humidity over ice (kg/kg), t in deg C, p in mb."""
    e = esat_ice(t, p)
    return 0.622 * e / (p - 0.378 * e)


def dqsat_ice_dt(t, p, dt=0.01):
    """d(qsat_ice)/dT by centred difference (kg kg-1 K-1)."""
    return (qsat_ice(t + dt, p) - qsat_ice(t - dt, p)) / (2 * dt)


# --- Gravity, air density, latent heat --------------------------------------------------
_G_EQUATOR = 9.7803267715
_G_SERIES = (0.0052790414, 0.0000232718, 0.0000001262, 0.0000000007)


def gravity(lat):
    """Normal gravity (m s-2) at latitude ``lat`` (degrees).

    g = g_e (1 + c1 sin^2 phi + c2 sin^4 phi + c3 sin^6 phi + c4 sin^8 phi),
    g_e = 9.7803267715 m s-2, with the coefficients used by COARE 3.5 and v2
    (``meteo.grv``). v2 gives no reference; the numbers are those of the GRS 1980
    normal-gravity series (citation from memory, unverified). pycoare uses a closed
    WGS84 formula instead, which differs by up to about 1.5e-7 relative.
    """
    s2 = np.sin(np.deg2rad(np.asarray(lat, dtype=float))) ** 2
    c1, c2, c3, c4 = _G_SERIES
    return _G_EQUATOR * (1 + c1 * s2 + c2 * s2**2 + c3 * s2**3 + c4 * s2**4)


def air_density(t, p, q):
    """Moist air density (kg m-3) as COARE computes it.

    rho = 100 p / (R_d (t + 273.16) (1 + 0.61 q)), p in mb, t in deg C, q in kg/kg.
    """
    t = np.asarray(t, dtype=float)
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    return p * 100.0 / (RGAS * (t + TDK) * (1 + 0.61 * q))


def air_viscosity(t):
    """Kinematic viscosity of air (m2 s-1) at temperature t (deg C), as COARE computes it.

        nu = 1.326e-5 (1 + 6.542e-3 t + 8.301e-6 t^2 - 4.84e-9 t^3)

    The same polynomial as in pycoare's COARE 3.5 and 3.6 (tested against it). COARE
    takes it from Andreas (1989, CRREL Report 89-11); that citation is from memory and
    unverified. Values: 1.326e-5 at 0 C, 1.240e-5 at -10 C, 1.156e-5 at -20 C.
    """
    t = np.asarray(t, dtype=float)
    return 1.326e-5 * (1 + 6.542e-3 * t + 8.301e-6 * t**2 - 4.84e-9 * t**3)


def latent_heat_vaporisation(ts):
    """Latent heat of vaporisation (J kg-1) at water temperature ts (deg C), as COARE 3.5.

    Le = (2.501 - 0.00237 ts) 1e6.
    """
    return (2.501 - 0.00237 * np.asarray(ts, dtype=float)) * 1e6
