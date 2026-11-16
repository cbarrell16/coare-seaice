"""Limits of the legacy interface: SIC = 0, SIC = 1, no ice, argument checks."""

import numpy as np
import pytest

from coare_seaice.legacy import OUTPUT_COLUMNS, coare35vn_ice

C = {name: i for i, name in enumerate(OUTPUT_COLUMNS)}


def _forcing(n=40):
    rng = np.random.default_rng(3)
    return dict(
        u=rng.uniform(0.5, 25, n), t=rng.uniform(-30, 5, n), rh=rng.uniform(50, 100, n),
        ts=rng.uniform(-1.9, 6, n), P=rng.uniform(960, 1040, n), lat=rng.uniform(60, 82, n),
        Rs=rng.uniform(0, 200, n), Rl=rng.uniform(150, 320, n),
    )  # fmt: skip


def _call(sic, **kw):
    f = _forcing()
    u, t, rh, ts = f.pop("u"), f.pop("t"), f.pop("rh"), f.pop("ts")
    return coare35vn_ice(u, t, rh, ts, zu=36, zt=22, zq=22, ice_conc=sic, verbose=False,
                         **f, **kw)  # fmt: skip


def test_no_ice_gives_38_columns():
    A = _call(None)
    assert A.shape == (40, 38)
    assert len(OUTPUT_COLUMNS) == 56


def test_open_ocean_scheme_gives_38_columns():
    assert _call(0.5, ice_scheme="open_ocean").shape == (40, 38)


@pytest.mark.parametrize("method", ["air_freezing", "seb_steady"])
def test_sic_zero_is_open_water(method):
    A = _call(0.0, t_ice_method=method)
    assert A.shape == (40, 56)
    assert np.all(A[:, C["tau_form"]] == 0.0)
    for w, b in (("tau", "tau_blend"), ("hsb", "hsb_blend"), ("hlb", "hlb_blend")):
        np.testing.assert_array_equal(A[:, C[b]], A[:, C[w]])


@pytest.mark.parametrize("method", ["air_freezing", "seb_steady"])
def test_sic_one_is_ice(method):
    A = _call(1.0, t_ice_method=method)
    # v2 clips the concentration to 1 - 1e-6 inside the form-drag formula, so the form
    # drag at SIC = 1 is tiny but not exactly zero.
    assert np.all(np.abs(A[:, C["tau_form"]]) < 1e-9)
    np.testing.assert_allclose(A[:, C["tau_blend"]], A[:, C["tau_ice"]], rtol=0, atol=1e-9)
    np.testing.assert_array_equal(A[:, C["hsb_blend"]], A[:, C["hsb_ice"]])
    np.testing.assert_array_equal(A[:, C["hlb_blend"]], A[:, C["hlb_ice"]])


def test_sic_nan_gives_nan_blends():
    A = _call(np.nan, t_ice_method="air_freezing")
    for col in ("tau_blend", "hsb_blend", "hlb_blend", "tau_form"):
        assert np.all(np.isnan(A[:, C[col]]))
    assert np.all(np.isfinite(A[:, C["hsb_ice"]]))


def test_argument_checks():
    with pytest.raises(ValueError):
        coare35vn_ice([5.0, 6.0], [-5.0], [80.0, 80.0], [0.0, 0.0])
    with pytest.raises(ValueError):
        _call(0.5, t_ice_method="nope")
    with pytest.raises(ValueError):
        _call(0.5, t_ice_method="prescribed")
    with pytest.raises(ValueError):
        _call(0.5, t_ice_method="diffusion")
    with pytest.raises(ValueError):
        _call(0.5, rh_reference="air")
    with pytest.raises(ValueError):
        coare35vn_ice([5.0], [-5.0], [80.0], [0.0], zu=np.array([10.0]))


def test_inputs_not_modified():
    f = _forcing()
    before = {k: v.copy() for k, v in f.items()}
    coare35vn_ice(f["u"], f["t"], f["rh"], f["ts"], P=f["P"], lat=f["lat"], ice_conc=0.5,
                  t_ice_method="air_freezing", verbose=False)  # fmt: skip
    for k in f:
        np.testing.assert_array_equal(f[k], before[k])
