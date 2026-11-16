"""Task F10: gustiness (unstable) and windless (stable) terms, Andreas et al. (2010)."""

import numpy as np
import pytest

from coare_seaice import thermo
from coare_seaice.most import ANDREAS2010_LIGHT_WIND as BOTH
from coare_seaice.most import SurfaceLayerOptions, ice_surface_layer


def _solve(u, t, t_ice, options=None, nits=10):
    u, t = np.atleast_1d(u).astype(float), np.atleast_1d(t).astype(float)
    n = u.size
    P = np.full(n, 1000.0)
    q = thermo.specific_humidity(t, P, np.full(n, 80.0))
    kw = dict(
        u=u, t=t, Q_air=q, t_ice=np.full(n, float(t_ice)), z0_ice=5e-4, P=P,
        grav=np.full(n, 9.82), zu=np.full(n, 10.0), zt=np.full(n, 10.0), zq=np.full(n, 10.0),
        rhoa=thermo.air_density(t, P, q), cpa=thermo.CPA, von=thermo.VON, fdg=thermo.FDG,
        ta=t + thermo.TDK, ice_scheme="blended_a87", scalar_ratio=0.2, nits=nits,
    )  # fmt: skip
    if options is not None:
        kw["options"] = options
    return ice_surface_layer(**kw)


def test_off_by_default_and_off_reproduces_previous_results():
    u = [0.2, 1.0, 5.0, 12.0]
    a = _solve(u, [-10.0] * 4, -20.0)
    b = _solve(u, [-10.0] * 4, -20.0, SurfaceLayerOptions(gustiness=False, windless=False))
    for k in ("tau", "hsb", "hlb", "L"):
        np.testing.assert_array_equal(a[k], b[k])
    np.testing.assert_array_equal(a["S"], np.maximum(u, 0.5))


def test_windless_term_eq_2_6():
    u = np.array([0.0, 0.5, 2.2, 10.0])
    fl = _solve(u, [-10.0] * 4, -20.0, BOTH)  # stable: surface colder than the air
    np.testing.assert_allclose(fl["S"], u + 0.5 / np.cosh(u), rtol=1e-14)
    assert np.isclose(fl["S"][0], 0.5)  # sech(0) = 1
    # "when U_r reaches 2.2 m/s, S_r is already within 5% of U_r" (text after Eq. 2.6)
    assert 0.045 < fl["S"][2] / 2.2 - 1.0 < 0.05
    assert np.isfinite(fl["hsb"]).all() and (fl["hsb"] < 0).all()


def test_gustiness_eq_2_4_and_2_5():
    u = np.array([0.0, 1.0, 4.0, 12.0])
    fl = _solve(u, [-25.0] * 4, -5.0, BOTH, nits=60)  # unstable: surface warmer
    wstar = fl["usr"] * np.cbrt(-600.0 / (0.4 * fl["L"]))  # Eq. 2.5
    np.testing.assert_allclose(fl["S"], np.sqrt(u**2 + (1.25 * wstar) ** 2), rtol=1e-6)
    # calm: free convection, a finite upward flux carried by the gusts alone
    assert 1.0 < fl["S"][0] < 3.0 and fl["hsb"][0] > 100.0 and not fl["zeta_limited"][0]
    assert fl["S"][3] / 12.0 - 1.0 < 0.03  # small in strong wind, even 20 K unstable
    rho = thermo.air_density(-25.0, 1000.0, thermo.specific_humidity(-25.0, 1000.0, 80.0))
    np.testing.assert_allclose(fl["tau"], rho * fl["usr"] ** 2, rtol=1e-12)  # Eq. 2.1a


def test_each_term_acts_only_in_its_own_regime():
    only_gust = SurfaceLayerOptions(gustiness=True, zeta_min=-1000.0)
    only_calm = SurfaceLayerOptions(windless=True)
    stable = dict(u=[1.0], t=[-10.0], t_ice=-20.0)
    unstable = dict(u=[1.0], t=[-25.0], t_ice=-5.0)
    assert _solve(**stable, options=only_gust)["hsb"][0] == _solve(**stable)["hsb"][0]
    assert _solve(**unstable, options=only_calm)["hsb"][0] == _solve(**unstable)["hsb"][0]
    assert _solve(**stable, options=only_calm)["hsb"][0] != _solve(**stable)["hsb"][0]
    assert _solve(**unstable, options=only_gust)["hsb"][0] != _solve(**unstable)["hsb"][0]


@pytest.mark.parametrize("t,t_ice", [(-10.0, -20.0), (-25.0, -5.0)])
def test_ten_iterations_are_enough_with_the_terms_on(t, t_ice):
    u = np.array([0.3, 1.0, 3.0, 8.0])
    a = _solve(u, [t] * 4, t_ice, BOTH)
    b = _solve(u, [t] * 4, t_ice, BOTH, nits=200)
    ok = ~a["not_converged"]
    # the flag bounds the last step, not the remaining error: allow 0.5 W m-2 (0.2 %)
    np.testing.assert_allclose(a["hsb"][ok], b["hsb"][ok], atol=0.5)
    assert ok.sum() >= 2


def test_nan_propagates_with_the_terms_on():
    fl = _solve([np.nan, 3.0], [-10.0, np.nan], -20.0, BOTH)
    assert np.isnan(fl["hsb"]).all() and np.isnan(fl["tau"]).all()


def test_default_unstable_limit_throttles_the_gust_term_and_is_flagged():
    kw = dict(u=[0.3], t=[-25.0], t_ice=-5.0, nits=100)
    limited = _solve(options=SurfaceLayerOptions(gustiness=True), **kw)
    free = _solve(options=BOTH, **kw)
    assert limited["zeta_limited"][0] and not free["zeta_limited"][0]
    assert limited["hsb"][0] < 0.5 * free["hsb"][0]
