"""Findings of the physics review of release 0.4.0."""

import numpy as np
import pytest

import coare_seaice as cs
from coare_seaice.engine import relative_humidity_over_ice
from coare_seaice.uncertainty import ensemble, normal

KW = dict(zu=10, zt=10, zq=10, lat=75)


def test_anchored_water_tile_ignores_the_ice_stability_options():
    tiles = cs.two_tile(water=dict(z0=1e-4, z0t=1e-5, z0q=1e-5))
    args = (2.0, -2.0, 80.0, 1000.0, -8.0, 0.5)  # stable over water
    a = cs.ice_fluxes(*args, tiles=tiles, jcool=0, **KW)
    b = cs.ice_fluxes(*args, tiles=tiles, jcool=0, stability="sheba", **KW)
    c = cs.ice_fluxes(*args, tiles=tiles, jcool=0, options=cs.ANDREAS2010_LIGHT_WIND, **KW)
    assert a.water.hsb[0] == b.water.hsb[0] == c.water.hsb[0]


def test_nan_blended_flux_is_always_flagged():
    no_rad = cs.ice_fluxes(8.0, -12.0, 80, 1000, -1.5, [0.5, 1.0], **KW)  # cool skin on
    assert np.isnan(no_rad.blended.hsb).all() and no_rad.qc["missing_input"].all()
    sst_nan = cs.ice_fluxes(8.0, -12.0, 80, 1000, np.nan, 1.0, jcool=0, **KW)
    assert np.isnan(sst_nan.blended.hsb[0]) and sst_nan.qc["missing_input"][0]
    with pytest.raises(ValueError):
        cs.three_tile(np.nan)


def test_humidity_flag_respects_rh_reference():
    over_ice = cs.ice_fluxes(8.0, -20.0, 100.0, 1000, -1.5, 0.95, jcool=0,
                             rh_reference="ice", **KW)  # fmt: skip
    assert not over_ice.qc["humidity_far_from_ice_saturation"][0]
    assert np.isnan(relative_humidity_over_ice(2.0, 1000.0, 90.0))


def test_ensemble_reports_valid_members_and_keeps_speed_positive():
    inp = dict(u=1.0, t=-12.0, rh=80.0, p=1000.0, sst=-1.5, sic=0.5)
    ens = ensemble(inp, dict(u=normal(2.0)), n=40, seed=3, jcool=0, **KW)
    assert (ens.n_valid["blended.hsb"] == 40).all()


def test_water_tile_rejects_named_roughness():
    with pytest.raises(ValueError, match="number or an array"):
        cs.two_tile(water=dict(z0="iceland_sea", z0t=1e-5, z0q=1e-5))
