"""Task F8: measurement heights may be arrays and any real number type."""

import numpy as np
import pytest

from coare_seaice.legacy import coare35vn_ice


def _inputs(n, seed=3):
    rng = np.random.default_rng(seed)
    args = (rng.uniform(2, 20, n), rng.uniform(-25, 2, n), rng.uniform(60, 98, n),
            rng.uniform(-1.8, 3, n))  # fmt: skip
    kw = dict(P=rng.uniform(970, 1030, n), lat=rng.uniform(65, 80, n),
              Rs=rng.uniform(0, 100, n), Rl=rng.uniform(180, 300, n),
              ice_conc=rng.uniform(0, 1, n), verbose=False)  # fmt: skip
    return args, kw


@pytest.mark.parametrize("method", ["air_freezing", "seb_steady"])
def test_per_sample_heights_equal_a_loop_over_scalars(method):
    n = 12
    args, kw = _inputs(n)
    rng = np.random.default_rng(9)
    zu, zt, zq = rng.uniform(20, 90, n), rng.uniform(15, 90, n), rng.uniform(15, 90, n)
    A = coare35vn_ice(*args, zu=zu, zt=zt, zq=zq, t_ice_method=method, **kw)
    for i in range(n):
        one = {k: (v[i : i + 1] if isinstance(v, np.ndarray) else v) for k, v in kw.items()}
        B = coare35vn_ice(*(a[i : i + 1] for a in args), zu=float(zu[i]), zt=float(zt[i]),
                          zq=float(zq[i]), t_ice_method=method, **one)  # fmt: skip
        np.testing.assert_allclose(A[i], B[0], rtol=1e-10, atol=1e-12)


@pytest.mark.parametrize("kind", [np.float32, np.float64, np.int64, int, float])
def test_any_real_number_type(kind):
    args, kw = _inputs(8)
    ref = coare35vn_ice(*args, zu=36, zt=22, zq=22, t_ice_method="air_freezing", **kw)
    out = coare35vn_ice(*args, zu=kind(36), zt=kind(22), zq=kind(22),
                        t_ice_method="air_freezing", **kw)  # fmt: skip
    np.testing.assert_array_equal(out, ref)


def test_bad_heights_are_rejected():
    args, kw = _inputs(4)
    with pytest.raises(ValueError, match="zt must be positive"):
        coare35vn_ice(*args, zu=10, zt=0.0, zq=10, **kw)
    with pytest.raises(ValueError, match="wrong length"):
        coare35vn_ice(*args, zu=np.array([10.0, 12.0]), zt=10, zq=10, **kw)
