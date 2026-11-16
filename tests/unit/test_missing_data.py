"""Task F7: NaN in gives NaN out; nothing is filled silently."""

import numpy as np
import pytest

from coare_seaice.legacy import OUTPUT_COLUMNS as COLS
from coare_seaice.legacy import coare35vn_ice
from coare_seaice.surface_temperature import surface_albedo

N = 5
BASE = dict(u=7.0, t=-12.0, rh=80.0, ts=-1.0, P=1000.0, Rs=50.0, Rl=230.0, lat=75.0,
            ice_conc=0.6, h_ice=1.0, h_snow=0.2)  # fmt: skip
TURB_ICE = ("tau_ice", "hsb_ice", "hlb_ice", "R_star", "z0T_ice", "z0q_ice")
BLEND = ("tau_blend", "hsb_blend", "hlb_blend")
SEB = ("T_ice", "F_cond", "LWnet_ice", "SEB_res")
# input -> outputs that must be NaN, for each surface-temperature method
EXPECT = {
    "air_freezing": {
        "u": TURB_ICE + BLEND, "t": TURB_ICE + BLEND + ("T_ice",), "rh": TURB_ICE + BLEND,
        "P": TURB_ICE + BLEND, "lat": TURB_ICE + BLEND, "ts": BLEND, "Rs": BLEND + ("SWnet_ice",),
        "Rl": BLEND + ("LWnet_ice",), "ice_conc": BLEND + ("tau_form",),
        "h_ice": ("F_cond", "SEB_res"), "h_snow": ("F_cond", "SEB_res", "SWnet_ice"),
    },
    "seb_steady": {
        "u": TURB_ICE + BLEND + SEB, "t": TURB_ICE + BLEND + SEB, "rh": TURB_ICE + BLEND + SEB,
        "P": TURB_ICE + BLEND + SEB, "lat": TURB_ICE + BLEND + SEB, "ts": BLEND,
        "Rs": TURB_ICE + BLEND + SEB + ("SWnet_ice",), "Rl": TURB_ICE + BLEND + SEB,
        "ice_conc": BLEND + ("tau_form",), "h_ice": TURB_ICE + ("hsb_blend", "hlb_blend") + SEB,
        "h_snow": TURB_ICE + ("hsb_blend", "hlb_blend", "SWnet_ice") + SEB,
    },
}  # fmt: skip


def _run(method, bad=None, **extra):
    d = {k: np.full(N, v) for k, v in BASE.items()}
    if bad is not None:
        d[bad][2] = np.nan
    args = [d.pop(k) for k in ("u", "t", "rh", "ts")]
    return coare35vn_ice(*args, zu=10, zt=10, zq=10, t_ice_method=method, verbose=False,
                         return_diagnostics=True, **d, **extra)  # fmt: skip


@pytest.mark.parametrize("method", sorted(EXPECT))
@pytest.mark.parametrize("name", sorted(BASE))
def test_nan_in_nan_out(method, name):
    A, diag = _run(method, bad=name)
    clean, _ = _run(method)
    others = [0, 1, 3, 4]
    np.testing.assert_array_equal(A[others], clean[others])  # other samples untouched
    for col in EXPECT[method][name]:
        assert np.isnan(A[2, COLS.index(col)]), f"{col} should be NaN when {name} is"
    seb_only = ("Rs", "Rl", "h_ice", "h_snow")  # not needed by air_freezing
    if name != "ts" and not (method == "air_freezing" and name in seb_only):
        assert diag["missing_input"][2] and not diag["missing_input"][others].any()


def test_prescribed_t_ice_is_not_filled_unless_asked():
    t_ice = np.array([-10.0, -9.0, np.nan, -8.0, -7.0])
    A, diag = _run("prescribed", t_ice=t_ice, t_ice_fallback=None)
    for col in TURB_ICE + BLEND + ("T_ice",):
        assert np.isnan(A[2, COLS.index(col)])
    assert diag["missing_input"][2] and not diag["t_ice_filled"].any()
    B, diag = _run("prescribed", t_ice=t_ice, t_ice_fallback="air_freezing")
    assert np.isfinite(B[2, COLS.index("hsb_ice")])
    assert list(diag["t_ice_filled"]) == [False, False, True, False, False]


def test_diffusion_flags_samples_where_the_column_is_not_running():
    sic = np.array([0.9, 0.9, 0.05, 0.9, 0.9])  # out of ice at sample 2
    d = {k: np.full(N, v) for k, v in BASE.items()}
    d["ice_conc"] = sic
    args = [d.pop(k) for k in ("u", "t", "rh", "ts")]
    A, diag = coare35vn_ice(*args, zu=10, zt=10, zq=10, t_ice_method="diffusion",
                            time_s=600.0 * np.arange(N), verbose=False,
                            return_diagnostics=True, **d)  # fmt: skip
    assert list(diag["t_ice_filled"]) == [False, False, True, False, False]
    assert np.isfinite(A[:, COLS.index("T_ice")]).all()


def test_albedo_is_nan_for_missing_snow_depth():
    np.testing.assert_array_equal(surface_albedo([0.0, 0.1, np.nan]), [0.65, 0.85, np.nan])


@pytest.mark.parametrize("name", ["h_ice", "h_snow", "albedo"])
def test_diffusion_does_not_fill_a_missing_slab_value_with_the_median(name):
    # Found in the physics review of release 0.3.0.
    d = {k: np.full(N, v) for k, v in BASE.items()}
    d["albedo"] = np.full(N, 0.85)
    d[name][2] = np.nan
    args = [d.pop(k) for k in ("u", "t", "rh", "ts")]
    A, diag = coare35vn_ice(*args, zu=10, zt=10, zq=10, t_ice_method="diffusion",
                            time_s=600.0 * np.arange(N), verbose=False,
                            return_diagnostics=True, **d)  # fmt: skip
    for col in ("T_ice", "hsb_ice", "hsb_blend"):
        assert np.isnan(A[2, COLS.index(col)])
    assert np.isfinite(A[[0, 1, 3, 4], COLS.index("hsb_ice")]).all()
    assert diag["missing_input"][2] and not diag["t_ice_filled"][2]
