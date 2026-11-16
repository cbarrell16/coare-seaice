"""COARE 3.5 stability (profile) functions psi(zeta), zeta = z / L.

Stable side (zeta >= 0): the Beljaars-Holtslag type form used in COARE 3.5.
Unstable side (zeta < 0): a blend of a Kansas-type form and a free-convection form,
weighted by f = zeta^2 / (1 + zeta^2).
COARE references as given in v2: Fairall et al. (2003), Edson et al. (2013). pycoare's
source comments attribute the stable form to Beljaars & Holtslag (1991) and the
convective form to Grachev et al. (2000).

The numerical constants are those of the COARE 3.5 code used by v2 (e.g. exponent
0.3333, constants 14.28 and 8.525 in psit_26). pycoare uses exact fractions (1/3, 2/3)
in places, so its psi values differ from these: by less than 1e-3 for psi_u; for psi_t
by about 1e-4 near neutral, up to about 5e-3 for |zeta| <= 1, and more in very stable
air (about 1.3 at zeta = 1000).
NaN in gives NaN out.
"""

from __future__ import annotations

import numpy as np

_DZETA_MAX = 50.0  # cap on 0.35 zeta inside exp(-0.35 zeta)


def _free_convective(x):
    """Free-convection psi for x = (1 - c zeta)^0.3333."""
    return (
        1.5 * np.log((1.0 + x + x**2) / 3.0)
        - np.sqrt(3.0) * np.arctan((1.0 + 2.0 * x) / np.sqrt(3.0))
        + np.pi / np.sqrt(3.0)
    )


def _blend(zeta, psi_kansas, psi_convective):
    f = zeta**2 / (1.0 + zeta**2)
    return (1 - f) * psi_kansas + f * psi_convective


def psiu_26(zeta):
    """Velocity stability function psi_u(zeta), COARE 3.5 ('psiu_26').

    Stable:   psi = -(a zeta + b (zeta - c/d) exp(-min(d zeta, 50)) + b c/d),
              a = 0.7, b = 3/4, c = 5, d = 0.35.
    Unstable: Kansas part with x = (1 - 15 zeta)^(1/4),
              psi_K = 2 ln((1+x)/2) + ln((1+x^2)/2) - 2 atan(x) + pi/2;
              convective part with x = (1 - 10.15 zeta)^0.3333,
              psi_C = 1.5 ln((1+x+x^2)/3) - sqrt(3) atan((1+2x)/sqrt(3)) + pi/sqrt(3).
    """
    zeta = np.asarray(zeta, dtype=float)
    z = np.array(zeta, ndmin=1).ravel()
    psi = np.empty_like(z)

    unstable = z < 0
    zs = z[~unstable]  # stable, neutral and NaN
    a, b, c, d = 0.7, 3.0 / 4.0, 5.0, 0.35
    dz = np.minimum(d * zs, _DZETA_MAX)
    psi[~unstable] = -(a * zs + b * (zs - c / d) * np.exp(-dz) + b * c / d)

    zu = z[unstable]
    x = (1 - 15 * zu) ** 0.25
    psi_k = 2.0 * np.log((1.0 + x) / 2.0) + np.log((1.0 + x * x) / 2.0) - 2.0 * np.arctan(x)
    psi_k = psi_k + np.pi / 2.0
    psi_c = _free_convective((1 - 10.15 * zu) ** 0.3333)
    psi[unstable] = _blend(zu, psi_k, psi_c)
    return psi.reshape(zeta.shape)


def psit_26(zeta):
    """Temperature (and humidity) stability function psi_t(zeta), COARE 3.5 ('psit_26').

    Stable:   psi = -((1 + 0.6667 zeta)^1.5 + 0.6667 (zeta - 14.28) exp(-min(0.35 zeta, 50))
              + 8.525).
    Unstable: Kansas part with x = (1 - 15 zeta)^(1/2), psi_K = 2 ln((1+x)/2);
              convective part as in :func:`psiu_26` with x = (1 - 34.15 zeta)^0.3333.
    """
    zeta = np.asarray(zeta, dtype=float)
    z = np.array(zeta, ndmin=1).ravel()
    psi = np.empty_like(z)

    unstable = z < 0
    zs = z[~unstable]
    dz = np.minimum(0.35 * zs, _DZETA_MAX)
    psi[~unstable] = -((1 + 0.6667 * zs) ** 1.5 + 0.6667 * (zs - 14.28) * np.exp(-dz) + 8.525)

    zu = z[unstable]
    psi_k = 2 * np.log((1 + (1 - 15 * zu) ** 0.5) / 2.0)
    psi_c = _free_convective((1 - 34.15 * zu) ** 0.3333)
    psi[unstable] = _blend(zu, psi_k, psi_c)
    return psi.reshape(zeta.shape)


# ---- SHEBA functions for stable conditions (Grachev et al. 2007) -----------------------
_AM, _BM = 5.0, 5.0 / 6.5  # a_m, b_m = a_m / 6.5
_AH, _BH, _CH = 5.0, 5.0, 3.0  # a_h, b_h, c_h


def phim_sheba(zeta):
    """Non-dimensional wind gradient, stable side (Grachev et al. 2007, Eq. 9a).

    phi_m = 1 + a_m zeta (1 + zeta)^(1/3) / (1 + b_m zeta), a_m = 5, b_m = a_m / 6.5.
    """
    zeta = np.asarray(zeta, dtype=float)
    return 1.0 + _AM * zeta * np.cbrt(1.0 + zeta) / (1.0 + _BM * zeta)


def phih_sheba(zeta):
    """Non-dimensional temperature gradient, stable side (Grachev et al. 2007, Eq. 9b).

    phi_h = 1 + (a_h zeta + b_h zeta^2) / (1 + c_h zeta + zeta^2), a_h = b_h = 5, c_h = 3.
    """
    zeta = np.asarray(zeta, dtype=float)
    return 1.0 + (_AH * zeta + _BH * zeta**2) / (1.0 + _CH * zeta + zeta**2)


def _psim_sheba_stable(zeta):
    """Grachev et al. (2007) Eq. 12, for zeta >= 0."""
    x = np.cbrt(1.0 + zeta)
    B = np.cbrt((1.0 - _BM) / _BM)
    s3 = np.sqrt(3.0)
    return -3.0 * _AM / _BM * (x - 1.0) + _AM * B / (2.0 * _BM) * (
        2.0 * np.log((x + B) / (1.0 + B))
        - np.log((x**2 - x * B + B**2) / (1.0 - B + B**2))
        + 2.0 * s3 * (np.arctan((2.0 * x - B) / (s3 * B)) - np.arctan((2.0 - B) / (s3 * B)))
    )


def _psih_sheba_stable(zeta):
    """Grachev et al. (2007) Eq. 13, for zeta >= 0."""
    B = np.sqrt(_CH**2 - 4.0)  # sqrt(5)
    return -_BH / 2.0 * np.log(1.0 + _CH * zeta + zeta**2) + (-_AH / B + _BH * _CH / (2.0 * B)) * (
        np.log((2.0 * zeta + _CH - B) / (2.0 * zeta + _CH + B)) - np.log((_CH - B) / (_CH + B))
    )


def _stable_or_coare(zeta, stable_fn, coare_fn):
    zeta = np.asarray(zeta, dtype=float)
    stable = zeta >= 0
    with np.errstate(invalid="ignore"):
        out = np.where(stable, stable_fn(np.where(stable, zeta, 0.0)), coare_fn(zeta))
    return np.where(np.isnan(zeta), np.nan, out)


def psim_sheba(zeta):
    """psi_m: SHEBA form for zeta >= 0 (Grachev et al. 2007, Eq. 12), COARE for zeta < 0.

    Grachev et al. give stable functions only. They are fitted to SHEBA data and stated
    to be valid for a bulk Richardson number below about 0.2. Compared with COARE's
    functions: psi_h is less negative at all zeta > 0 (-10 against -30 at zeta = 10), so
    scalar transfer is stronger; psi_m is slightly less negative below zeta of about 3,
    more negative from there to about 60 (-22 against -18 at zeta = 10), and less
    negative beyond.
    """
    return _stable_or_coare(zeta, _psim_sheba_stable, psiu_26)


def psih_sheba(zeta):
    """psi_h: SHEBA form for zeta >= 0 (Grachev et al. 2007, Eq. 13), COARE for zeta < 0."""
    return _stable_or_coare(zeta, _psih_sheba_stable, psit_26)


# name -> (psi for momentum, psi for scalars)
STABILITY_FUNCTIONS = {"coare": (psiu_26, psit_26), "sheba": (psim_sheba, psih_sheba)}
