"""Task T6: SHEBA stability functions (Grachev et al. 2007)."""

import numpy as np
import pytest

import coare_seaice as cs
from coare_seaice.stability import (
    phih_sheba,
    phim_sheba,
    psih_sheba,
    psim_sheba,
    psit_26,
    psiu_26,
)

ZETA = np.array([1e-3, 0.01, 0.1, 0.5, 1.0, 3.0, 10.0, 30.0, 100.0])


@pytest.mark.parametrize("psi,phi", [(psim_sheba, phim_sheba), (psih_sheba, phih_sheba)])
def test_psi_is_the_integral_of_phi(psi, phi):
    # phi = 1 - zeta dpsi/dzeta (from Eq. 6): check Eqs. 12 and 13 against Eqs. 9a and 9b.
    h = ZETA * 1e-5
    dpsi = (psi(ZETA + h) - psi(ZETA - h)) / (2 * h)
    np.testing.assert_allclose(1.0 - ZETA * dpsi, phi(ZETA), rtol=1e-6)
    assert psi(0.0) == 0.0


def test_published_limits():
    assert np.isclose(phim_sheba(1e6) / 1e6 ** (1 / 3), 6.5, rtol=5e-3)  # phi_m -> 6.5 zeta^1/3
    assert np.isclose(phih_sheba(1e6), 6.0, rtol=1e-4)  # phi_h -> 1 + b_h = 6
    np.testing.assert_allclose(phim_sheba(1e-4), 1 + 5e-4, rtol=1e-6)  # 1 + 5 zeta near neutral
    np.testing.assert_allclose(phih_sheba(1e-4), 1 + 5e-4, rtol=1e-6)


def test_unstable_side_is_coare_and_nan_propagates():
    z = np.array([-5.0, -0.5, -0.01])
    np.testing.assert_array_equal(psim_sheba(z), psiu_26(z))
    np.testing.assert_array_equal(psih_sheba(z), psit_26(z))
    assert np.isnan(psim_sheba(np.nan)) and np.isnan(psih_sheba(np.array([np.nan, 1.0]))[0])


def test_comparison_with_coare_as_documented():
    z = np.array([0.1, 1.0, 10.0, 100.0])
    assert np.all(psih_sheba(z) > psit_26(z))  # scalar transfer stronger at all zeta
    assert psim_sheba(1.0) > psiu_26(1.0) and psim_sheba(10.0) < psiu_26(10.0)
    assert psim_sheba(100.0) > psiu_26(100.0)
    np.testing.assert_allclose([psim_sheba(10.0), psih_sheba(10.0)], [-21.8, -10.3], atol=0.05)


def test_option_in_ice_fluxes(inputs):
    sic = np.where(np.isfinite(inputs["sic"]), np.clip(inputs["sic"], 0, 1), np.nan)
    args = (inputs["u"], inputs["t"], inputs["rh"], inputs["P"], inputs["ts"], sic)
    kw = dict(zu=36, zt=22, zq=22, lat=inputs["lat"], lw_down=inputs["Rl"], sw_down=inputs["Rs"])
    coare = cs.ice_fluxes(*args, **kw)
    same = cs.ice_fluxes(*args, stability="coare", **kw)
    sheba = cs.ice_fluxes(*args, stability="sheba", **kw)
    np.testing.assert_array_equal(coare.ice.hsb, same.ice.hsb)
    np.testing.assert_array_equal(coare.water.hsb, sheba.water.hsb)  # open water unchanged
    stable, unstable = coare.ice.zeta > 0.05, coare.ice.zeta < -0.05
    np.testing.assert_array_equal(coare.ice.hsb[unstable], sheba.ice.hsb[unstable])
    assert np.nanmax(np.abs(coare.ice.hsb[stable] - sheba.ice.hsb[stable])) > 0.1
    with pytest.raises(ValueError, match="stability must be"):
        cs.ice_fluxes(*args, stability="kansas", **kw)
