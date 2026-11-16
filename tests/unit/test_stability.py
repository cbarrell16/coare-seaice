"""stability: psi functions equal v2 meteo.psiu_26 / meteo.psit_26."""

import numpy as np
import pytest

from coare_seaice.stability import psit_26, psiu_26
from conftest import assert_port, load_ref


@pytest.fixture(scope="module")
def ref():
    return load_ref("thermo_stability")


def test_psiu_26(ref):
    assert_port(psiu_26(ref["zeta"]), ref["psiu_26"])


def test_psit_26(ref):
    assert_port(psit_26(ref["zeta"]), ref["psit_26"])


@pytest.mark.parametrize("f", [psiu_26, psit_26])
def test_shape_and_nan(f):
    assert f(0.3).shape == ()
    assert f(np.zeros((2, 3))).shape == (2, 3)
    assert np.isnan(f(np.nan))


def test_psiu_neutral_is_zero():
    assert psiu_26(0.0) == pytest.approx(0.0, abs=1e-12)


def test_close_to_pycoare_where_expected(ref):
    # pycoare uses exact fractions and slightly different constants in places. This
    # documents the size of that known difference (not a property of v2).
    from pycoare import util

    z = ref["zeta"]
    ok = np.isfinite(z)
    assert np.max(np.abs(psiu_26(z[ok]) - util.psiu_26(z[ok]))) < 1e-3
    mild = ok & (np.abs(z) <= 1.0)
    assert np.max(np.abs(psit_26(z[mild]) - util.psit_26(z[mild]))) < 1e-2
