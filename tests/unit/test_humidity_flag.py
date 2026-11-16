"""Task T9: flag air far from ice saturation over ice. A flag only: nothing is adjusted."""

import numpy as np

import coare_seaice as cs
from coare_seaice.engine import relative_humidity_over_ice

KW = dict(zu=10, zt=10, zq=10, lat=75, lw_down=220.0, sw_down=0.0)


def test_relative_humidity_over_ice():
    rh_i = relative_humidity_over_ice(np.array([-10.0, -20.0, 0.0]), 1000.0, 100.0)
    np.testing.assert_allclose(rh_i, [110.2, 121.5, 100.0], atol=0.6)


def test_flag_is_set_only_over_ice_and_only_outside_the_range():
    #            ice sat.  dry    supersat.  dry but open water
    t = np.full(4, -15.0)
    rh = np.array([86.0, 55.0, 99.0, 55.0])
    sic = np.array([0.95, 0.95, 0.95, 0.3])
    res = cs.ice_fluxes(8.0, t, rh, 1000.0, -1.5, sic, **KW)
    flag = res.qc["humidity_far_from_ice_saturation"]
    assert list(flag) == [False, True, True, False]


def test_flag_does_not_change_any_value():
    t, sic = np.full(3, -15.0), np.full(3, 0.95)
    res = cs.ice_fluxes(8.0, t, np.array([86.0, 55.0, 99.0]), 1000.0, -1.5, sic, **KW)
    one = cs.ice_fluxes(8.0, -15.0, 55.0, 1000.0, -1.5, 0.95, **KW)
    assert res.blended.hlb[1] == one.blended.hlb[0]  # the dry, flagged sample is as computed
    assert np.isfinite(res.blended.hlb).all()
    assert not cs.ice_fluxes(8.0, -15.0, np.nan, 1000.0, -1.5, 0.95, **KW).qc[
        "humidity_far_from_ice_saturation"
    ][0]
