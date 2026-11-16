"""Task T3: tiles, fractions and flux blending."""

import numpy as np
import pytest

from coare_seaice.engine import Forcing, compute
from coare_seaice.legacy import OUTPUT_COLUMNS as COLS
from coare_seaice.legacy import coare35vn_ice
from coare_seaice.surface_temperature import SEB, AirFreezing
from coare_seaice.tiles import Tile, three_tile, tile_fractions, two_tile


@pytest.fixture(scope="module")
def forcing(inputs):
    n = inputs["u"].size
    full = lambda v: np.full(n, float(v))  # noqa: E731
    sic = np.where(np.isfinite(inputs["sic"]), np.clip(inputs["sic"], 0, 1), np.nan)
    return Forcing(
        u=inputs["u"], t=inputs["t"], rh=inputs["rh"], p=inputs["P"], sst=inputs["ts"],
        sic=sic, zu=full(36), zt=full(22), zq=full(22), lat=inputs["lat"],
        lw_down=inputs["Rl"], sw_down=inputs["Rs"], zi=full(600),
    )  # fmt: skip


def _legacy(inputs, method, **kw):
    return coare35vn_ice(
        inputs["u"], inputs["t"], inputs["rh"], inputs["ts"], P=inputs["P"], lat=inputs["lat"],
        Rs=inputs["Rs"], Rl=inputs["Rl"], zu=36, zt=22, zq=22, ice_conc=inputs["sic"],
        t_ice_method=method, verbose=False, **kw,
    )  # fmt: skip


@pytest.mark.parametrize("method,name", [(AirFreezing(), "air_freezing"), (SEB(), "seb_steady")])
def test_two_tiles_equal_legacy(forcing, inputs, method, name):
    res = compute(forcing, two_tile(), t_ice=method)
    A = _legacy(inputs, name)
    pairs = [("tau_blend", res.blended.tau), ("hsb_blend", res.blended.hsb),
             ("hlb_blend", res.blended.hlb), ("tau_form", res.blended.tau_form),
             ("tau_ice", res.ice.tau), ("hsb_ice", res.ice.hsb), ("T_ice", res.ice.t_surface),
             ("F_cond", res.ice.f_cond), ("SEB_res", res.ice.seb_residual),
             ("tau", res.water.tau), ("hsb", res.water.hsb)]  # fmt: skip
    for col, values in pairs:
        np.testing.assert_allclose(values, A[:, COLS.index(col)], rtol=1e-12, atol=1e-12,
                                   err_msg=col)  # fmt: skip


def test_three_tiles_with_no_thin_ice_equal_two_tiles(forcing):
    two = compute(forcing, two_tile(), t_ice=SEB())
    for of in ("water", "ice"):
        three = compute(forcing, three_tile(0.0, of=of), t_ice=SEB())
        for name in ("tau", "hsb", "hlb", "tau_form", "cdn10", "chn10"):
            np.testing.assert_allclose(getattr(three.blended, name), getattr(two.blended, name),
                                       rtol=1e-13, atol=1e-13, err_msg=name)  # fmt: skip
        assert np.nanmax(three.tiles["thin_ice"].fraction) == 0.0


def test_fractions():
    sic = np.array([0.0, 0.4, 1.0, np.nan])
    f = tile_fractions(three_tile(0.25, of="water"), sic)
    np.testing.assert_allclose(f["water"][:3], [0.75, 0.45, 0.0])
    np.testing.assert_allclose(f["thin_ice"][:3], [0.25, 0.15, 0.0])
    np.testing.assert_allclose(f["thick_ice"][:3], [0.0, 0.4, 1.0])
    f = tile_fractions(three_tile(0.25, of="ice"), sic)
    np.testing.assert_allclose(f["thin_ice"][:3], [0.0, 0.1, 0.25])
    np.testing.assert_allclose(f["thick_ice"][:3], [0.0, 0.3, 0.75])
    assert all(np.isnan(v[3]) for v in f.values())
    bad = [Tile("water", "water", lambda s: 1 - s), Tile("ice", "ice", lambda s: 0.5 * s)]
    with pytest.raises(ValueError, match="sum to one"):
        tile_fractions(bad, sic)
    with pytest.raises(ValueError):
        three_tile(1.5)


def test_blend_is_the_area_mean_plus_form_drag(forcing):
    res = compute(forcing, three_tile(0.5, of="water", thin=dict(t_surface=SEB())))
    tiles = res.tiles.values()
    for name in ("hsb", "hlb"):
        np.testing.assert_allclose(getattr(res.blended, name),
                                   sum(t.fraction * getattr(t, name) for t in tiles))  # fmt: skip
    np.testing.assert_allclose(
        res.blended.tau, sum(t.fraction * t.tau for t in tiles) + res.blended.tau_form
    )
    # thin bare ice is warmer than thick snow-covered ice at min(T_air, -1.8)
    cold = forcing.t < -10
    thin, thick = res.tiles["thin_ice"], res.tiles["thick_ice"]
    assert np.nanmedian(thin.t_surface[cold] - thick.t_surface[cold]) > 0


def test_thin_ice_from_the_water_fraction_adds_no_form_drag(forcing):
    two = compute(forcing, two_tile())
    three = compute(forcing, three_tile(0.6, of="water"))
    np.testing.assert_allclose(three.blended.tau_form, two.blended.tau_form, rtol=1e-13)


def test_per_tile_settings_and_flags(forcing):
    tiles = two_tile(z0=1e-2, scalar_scheme="ratio", t_surface=AirFreezing(-5.0))
    res = compute(forcing, tiles)
    np.testing.assert_allclose(res.ice.z0t, 0.2 * 1e-2)
    assert np.nanmax(res.ice.t_surface) == -5.0
    expected = {"missing_input", "rib_above_critical", "zeta_limited", "light_wind",
                "height_above_obukhov_length", "rstar_out_of_range", "not_converged",
                "t_surface_filled", "melting_surface"}  # fmt: skip
    assert set(res.qc) == expected and all(v.dtype == bool for v in res.qc.values())
    np.testing.assert_array_equal(res.qc["light_wind"], forcing.u < 1.0)
    assert res.qc["missing_input"][np.isnan(forcing.sic)].all()
    ds = res.to_xarray()
    assert "thick_ice_hsb" in ds and "qc_rib_above_critical" in ds
