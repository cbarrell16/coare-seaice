"""Task T1: named results, round trip to the legacy array, conversion to xarray."""

import numpy as np
import pytest

from coare_seaice.legacy import OUTPUT_COLUMNS, coare35vn_ice
from coare_seaice.results import BLENDED_FIELDS, QC_FLAGS, TILE_FIELDS, FluxResult, TileResult


@pytest.fixture(scope="module")
def legacy(inputs):
    kw = dict(P=inputs["P"], lat=inputs["lat"], Rs=inputs["Rs"], Rl=inputs["Rl"], zu=36, zt=22,
              zq=22, ice_conc=inputs["sic"], verbose=False, return_diagnostics=True)  # fmt: skip
    A, diag = coare35vn_ice(inputs["u"], inputs["t"], inputs["rh"], inputs["ts"],
                            t_ice_method="seb_steady", **kw)  # fmt: skip
    return A, diag, inputs["sic"]


def test_round_trip_to_the_legacy_array(legacy):
    A, diag, sic = legacy
    res = FluxResult.from_legacy_array(A, sic, diag)
    np.testing.assert_array_equal(res.to_legacy_array(), A)
    assert res.n == A.shape[0] and set(res.tiles) == {"water", "thick_ice"}


def test_names_and_signs(legacy):
    A, _, sic = legacy
    res = FluxResult.from_legacy_array(A, sic)
    c = OUTPUT_COLUMNS.index
    np.testing.assert_array_equal(res.blended.hsb, A[:, c("hsb_blend")])
    np.testing.assert_array_equal(res.water.tau, A[:, c("tau")])
    np.testing.assert_array_equal(res.ice.t_surface, A[:, c("T_ice")])
    np.testing.assert_allclose(res.ice.cdn10, A[:, c("Cdn_ice")] / 1e3)
    np.testing.assert_array_equal(res.tiles["thick_ice"].fraction, sic)


def test_to_xarray(legacy):
    A, diag, sic = legacy
    ds = FluxResult.from_legacy_array(A, sic, diag).to_xarray(time=np.arange(A.shape[0]))
    assert ds.sizes["sample"] == A.shape[0] and "time" in ds.coords
    assert ds["blended_hsb"].attrs["units"] == "W m-2"
    assert "positive upward" in ds["thick_ice_hsb"].attrs["long_name"]
    assert ds["thick_ice_qc_qc_rib"].dtype == bool
    np.testing.assert_array_equal(ds["water_tau"].values, A[:, 1])


def test_every_field_is_documented():
    fields = set(TileResult.__dataclass_fields__) - {"name", "kind", "qc", "extra"}
    assert fields == set(TILE_FIELDS)
    assert all(len(v) == 2 for v in BLENDED_FIELDS.values()) and "missing_input" in QC_FLAGS


def test_to_legacy_needs_two_tiles(legacy):
    A, _, sic = legacy
    res = FluxResult.from_legacy_array(A, sic)
    res.tiles["thin_ice"] = TileResult(name="thin_ice", kind="ice", fraction=0 * sic)
    with pytest.raises(ValueError):
        res.to_legacy_array()
