"""Task T4: blending with one aggregate surface temperature."""

import numpy as np
import pytest

from coare_seaice import thermo
from coare_seaice.engine import Forcing, compute
from coare_seaice.surface_temperature import SEB, Prescribed
from coare_seaice.tiles import three_tile, two_tile


def _forcing(n=40, seed=2, z=10.0, **over):
    rng = np.random.default_rng(seed)
    d = dict(
        u=rng.uniform(3, 18, n), t=rng.uniform(-25, -3, n), rh=rng.uniform(60, 95, n),
        p=rng.uniform(980, 1030, n), sst=np.full(n, -1.5), sic=rng.uniform(0.05, 0.95, n),
        zu=np.full(n, z), zt=np.full(n, z), zq=np.full(n, z), lat=np.full(n, 75.0),
        lw_down=rng.uniform(170, 290, n), sw_down=np.zeros(n), zi=np.full(n, 600.0), jcool=0,
    )  # fmt: skip
    d.update(over)
    return Forcing(**d)


def test_equals_flux_mode_when_all_tiles_share_one_temperature():
    F = _forcing()
    T = np.linspace(-20.0, -2.0, F.n)
    single = compute(F, three_tile(0.3), blend="single_temperature", t_surface=T)
    F_same = _forcing(sst=T)
    flux = compute(F_same, three_tile(0.3), t_ice=Prescribed(T))
    for name in ("tau", "hsb", "hlb", "tau_form"):
        np.testing.assert_allclose(getattr(single.blended, name), getattr(flux.blended, name),
                                   rtol=1e-12, err_msg=name)  # fmt: skip
    np.testing.assert_array_equal(single.blended.t_surface, T)
    assert all(np.allclose(r.t_surface, T) for r in single.tiles.values())


def test_neutral_limit_is_coefficient_blending_elvidge_2021_eq_12():
    # Strong wind and a small temperature difference: near neutral, at 10 m.
    n = 30
    t = np.linspace(-12.0, -6.0, n)
    F = _forcing(n=n, u=np.full(n, 20.0), t=t, rh=np.full(n, 85.0), p=np.full(n, 1000.0))
    T = t + 0.2  # 0.1 K warmer than the 10 m potential temperature
    res = compute(F, two_tile(), blend="single_temperature", t_surface=T)
    assert np.max(np.abs(res.ice.zeta)) < 0.01
    q = thermo.specific_humidity(t, F.p, F.rh)
    rho = thermo.air_density(t, F.p, q)
    u10n = res.water.extra["open_water"].u_n_rf
    chn = (1 - F.sic) * res.water.chn10 + F.sic * res.ice.chn10  # Eq. 12
    np.testing.assert_allclose(res.blended.chn10, chn, rtol=1e-12)
    expected = rho * thermo.CPA * chn * u10n * (T - (t + thermo.LAPSE_DRY * 10.0))
    np.testing.assert_allclose(res.blended.hsb, expected, rtol=0.03)


def test_default_aggregate_is_the_area_mean_of_the_tile_temperatures():
    F = _forcing()
    flux = compute(F, two_tile(), t_ice=SEB())
    single = compute(F, two_tile(), t_ice=SEB(), blend="single_temperature")
    np.testing.assert_allclose(single.blended.t_surface, flux.blended.t_surface, rtol=1e-12)
    # one cold aggregate surface loses less heat than warm water plus cold ice
    assert np.mean(single.blended.hsb) < np.mean(flux.blended.hsb)


def test_arguments():
    F = _forcing(n=5)
    with pytest.raises(ValueError, match="single_temperature"):
        compute(F, two_tile(), t_surface=-5.0)
    with pytest.raises(ValueError, match="blend must be"):
        compute(F, two_tile(), blend="coefficients")
    nan_t = np.array([-5.0, np.nan, -5.0, -5.0, -5.0])
    res = compute(F, two_tile(), blend="single_temperature", t_surface=nan_t)
    assert np.isnan(res.blended.hsb[1]) and res.qc["missing_input"][1]
    assert np.isfinite(res.blended.hsb[[0, 2, 3, 4]]).all()


def test_turbulent_fluxes_without_cool_skin_do_not_need_radiation():
    F = _forcing()
    F_nan = _forcing(lw_down=np.full(F.n, np.nan), sw_down=np.full(F.n, np.nan))
    a, b = compute(F, two_tile()), compute(F_nan, two_tile())
    for name in ("tau", "hsb", "hlb"):
        np.testing.assert_array_equal(getattr(a.blended, name), getattr(b.blended, name))
    assert np.isnan(b.ice.lw_net).all() and not b.qc["missing_input"].any()
    seb = compute(F_nan, two_tile(), t_ice=SEB())  # the energy balance does need it
    assert np.isnan(seb.blended.hsb).all() and seb.qc["missing_input"].all()
