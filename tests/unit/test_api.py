"""Task T5: the public ice_fluxes function."""

import numpy as np
import pytest

import coare_seaice as cs
from coare_seaice.legacy import OUTPUT_COLUMNS as COLS
from coare_seaice.legacy import coare35vn_ice


def _args(inputs):
    return (inputs["u"], inputs["t"], inputs["rh"], inputs["P"], inputs["ts"],
            np.where(np.isfinite(inputs["sic"]), np.clip(inputs["sic"], 0, 1), np.nan))  # fmt: skip


def _kw(inputs):
    return dict(zu=36, zt=22, zq=22, lat=inputs["lat"], lw_down=inputs["Rl"],
                sw_down=inputs["Rs"])  # fmt: skip


def test_defaults_are_the_simple_reference(inputs):
    res = cs.ice_fluxes(*_args(inputs), **_kw(inputs))
    A = coare35vn_ice(inputs["u"], inputs["t"], inputs["rh"], inputs["ts"], P=inputs["P"],
                      lat=inputs["lat"], Rs=inputs["Rs"], Rl=inputs["Rl"], zu=36, zt=22, zq=22,
                      ice_conc=inputs["sic"], t_ice_method="air_freezing",
                      verbose=False)  # fmt: skip
    for col, values in (("tau_blend", res.blended.tau), ("hsb_blend", res.blended.hsb),
                        ("hlb_blend", res.blended.hlb), ("T_ice", res.ice.t_surface)):  # fmt: skip
        np.testing.assert_allclose(values, A[:, COLS.index(col)], rtol=1e-12, atol=1e-12)
    assert list(res.tiles) == ["water", "thick_ice"] and res.meta["blend"] == "flux"
    np.testing.assert_array_equal(res.ice.z0, 5e-4)
    np.testing.assert_array_equal(res.ice.t_surface, np.minimum(inputs["t"], -1.8))


def test_surface_temperature_can_be_given_three_ways(inputs):
    a = cs.ice_fluxes(*_args(inputs), t_ice=cs.SEB(), **_kw(inputs))
    b = cs.ice_fluxes(*_args(inputs), t_ice="seb", **_kw(inputs))
    np.testing.assert_array_equal(a.blended.hsb, b.blended.hsb)
    T = np.full(inputs["u"].size, -9.0)
    c = cs.ice_fluxes(*_args(inputs), t_ice=T, **_kw(inputs))
    d = cs.ice_fluxes(*_args(inputs), t_ice=cs.Prescribed(-9.0), **_kw(inputs))
    np.testing.assert_array_equal(c.blended.hsb, d.blended.hsb)
    assert np.all(c.ice.t_surface == -9.0)


def test_scalars_and_lists():
    res = cs.ice_fluxes(8.0, -12.0, 80.0, 1000.0, -1.5, 0.5, zu=10, zt=10, zq=10, lat=75,
                        lw_down=230.0, sw_down=0.0)  # fmt: skip
    assert res.n == 1 and np.isfinite(res.blended.hsb[0])
    res = cs.ice_fluxes([8.0, 9.0], [-12.0, -5.0], [80, 85], 1000, -1.5, [0.5, 0.9], zu=36,
                        zt=22, zq=22, lat=75, lw_down=230, sw_down=0)  # fmt: skip
    assert res.blended.hsb.shape == (2,)


def test_input_checks():
    ok = dict(zu=10, zt=10, zq=10, lat=75)
    with pytest.raises(ValueError, match="one value per sample"):
        cs.ice_fluxes([8.0, 9.0], [-12.0, -5.0, 1.0], 80, 1000, -1.5, 0.5, **ok)
    with pytest.raises(ValueError, match="fraction, not per cent"):
        cs.ice_fluxes(8.0, -12.0, 80, 1000, -1.5, 50.0, **ok)
    with pytest.raises(ValueError, match="zt must be positive"):
        cs.ice_fluxes(8.0, -12.0, 80, 1000, -1.5, 0.5, zu=10, zt=-2, zq=10, lat=75)
    with pytest.raises(ValueError, match="unknown surface temperature method"):
        cs.ice_fluxes(8.0, -12.0, 80, 1000, -1.5, 0.5, t_ice="warm", **ok)
    with pytest.raises(ValueError, match="needs time"):
        cs.ice_fluxes(8.0, -12.0, 80, 1000, -1.5, 0.5, t_ice=cs.Column(), lw_down=230,
                      sw_down=0, **ok)  # fmt: skip


def test_no_default_radiation():
    # Cool skin on and no radiation given: the water fluxes are NaN, not computed with a
    # made-up value. With the cool skin off they do not need it.
    kw = dict(zu=10, zt=10, zq=10, lat=75)
    on = cs.ice_fluxes(8.0, -12.0, 80, 1000, -1.5, 0.5, **kw)
    off = cs.ice_fluxes(8.0, -12.0, 80, 1000, -1.5, 0.5, jcool=0, **kw)
    assert np.isnan(on.blended.hsb[0]) and np.isfinite(off.blended.hsb[0])
    assert np.isfinite(on.ice.hsb[0])  # the ice tile with AirFreezing does not need it


def test_inputs_are_not_modified(inputs):
    args = [np.array(a, dtype=float) for a in _args(inputs)]
    before = [a.copy() for a in args]
    cs.ice_fluxes(*args, t_ice=cs.SEB(), **_kw(inputs))
    for a, b in zip(args, before, strict=True):
        np.testing.assert_array_equal(a, b)


def test_light_wind_options_pass_through(inputs):
    a = cs.ice_fluxes(*_args(inputs), **_kw(inputs))
    b = cs.ice_fluxes(*_args(inputs), options=cs.ANDREAS2010_LIGHT_WIND, **_kw(inputs))
    assert np.nanmax(np.abs(a.ice.hsb - b.ice.hsb)) > 0.01
    np.testing.assert_array_equal(a.water.hsb, b.water.hsb)
