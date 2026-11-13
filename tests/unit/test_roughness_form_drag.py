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
    z0T, z0q, R = a87_scalar_roughness(ref["a87_z0"], ref["a87_usr"])
    assert_port(z0T, ref["a87_z0T"])
    assert_port(z0q, ref["a87_z0q"])
    assert_port(R, ref["a87_R_star"])


def test_reference_covers_all_regimes(ref):
    R = ref["a87_R_star"]
    assert (R <= R_SMOOTH).any()
    assert ((R > R_SMOOTH) & (R < R_ROUGH)).any()
    assert ((R >= R_ROUGH) & (R <= 1000)).any()
    assert (R > 1000).any()  # beyond the stated range (AUDIT F3)


def test_form_drag(ref):
    assert_port(lupkes_form_drag(ref["fd_A"], ref["fd_z0w"]), ref["fd_CD10nf"])


def test_form_drag_limits():
    cd = lupkes_form_drag(np.array([0.0, 0.009, np.nan, 1.0]), 2e-4)
    assert cd[0] == 0.0 and cd[1] == 0.0 and np.isnan(cd[2])
    assert 0.0 < cd[3] < 1e-11  # A is clipped to 1 - 1e-6, so not exactly zero at A = 1
