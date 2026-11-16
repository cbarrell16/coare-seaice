"""roughness and form_drag equal their v2 counterparts."""

import numpy as np
import pytest

from coare_seaice.form_drag import lupkes_form_drag
from coare_seaice.roughness import R_ROUGH, R_SMOOTH, a87_scalar_roughness
from conftest import assert_port, load_ref


@pytest.fixture(scope="module")
def ref():
    return load_ref("roughness_form_drag")


def test_a87_scalar_roughness(ref):
    z0T, z0q, R = a87_scalar_roughness(ref["a87_z0"], ref["a87_usr"], guard=False)
    assert_port(z0T, ref["a87_z0T"])
    assert_port(z0q, ref["a87_z0q"])
    assert_port(R, ref["a87_R_star"])


def test_reference_covers_all_regimes(ref):
    R = ref["a87_R_star"]
    assert (R <= R_SMOOTH).any()
    assert ((R > R_SMOOTH) & (R < R_ROUGH)).any()
    assert ((R >= R_ROUGH) & (R <= 1000)).any()
    assert (R > 1000).any()  # beyond the stated range (guarded since task F3)


def test_form_drag(ref):
    assert_port(lupkes_form_drag(ref["fd_A"], ref["fd_z0w"]), ref["fd_CD10nf"])


def test_form_drag_limits():
    cd = lupkes_form_drag(np.array([0.0, 0.009, np.nan, 1.0]), 2e-4)
    assert cd[0] == 0.0 and cd[1] == 0.0 and np.isnan(cd[2])
    assert 0.0 < cd[3] < 1e-11  # A is clipped to 1 - 1e-6, so not exactly zero at A = 1


def test_a87_guard_beyond_fitted_range():
    # Task F3: R* = 5000 (z0 = 0.05 m, u* = 1.5 m/s, nu = 1.5e-5). No extrapolation.
    z0, us = np.array([0.05, 0.05]), np.array([1.5, 0.3])
    z0T, z0q, R = a87_scalar_roughness(z0, us)
    assert np.isclose(R[0], 5000.0) and np.isclose(R[1], 1000.0)
    np.testing.assert_allclose(z0T[0], z0T[1], rtol=1e-12)  # held at the R* = 1000 value
    np.testing.assert_allclose(z0q[0], z0q[1], rtol=1e-12)
    lnr = np.log(1000.0)
    ratio = np.exp(0.317 - 0.565 * lnr - 0.183 * lnr**2)
    np.testing.assert_allclose(z0T[1] / z0[1], ratio, rtol=1e-12)
    z0T_x, _, _ = a87_scalar_roughness(z0, us, guard=False)
    assert z0T_x[0] < z0T[0] and z0T_x[1] == z0T[1]


def test_published_a87_values():
    # Andreas (1987) Table 1: smooth regime ln(z0T/z0) = 1.250, ln(z0q/z0) = 1.610.
    z0T, z0q, R = a87_scalar_roughness(np.array([1e-5]), np.array([0.1]))
    assert R[0] < R_SMOOTH
    np.testing.assert_allclose([z0T[0] / 1e-5, z0q[0] / 1e-5], np.exp([1.25, 1.61]), rtol=1e-12)
