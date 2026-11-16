"""Task T8: COARE 3.6 is selectable and wave inputs pass through."""

import numpy as np
import pytest

import coare_seaice as cs

rng = np.random.default_rng(0)
N = 120
ARGS = (rng.uniform(1, 22, N), rng.uniform(-25, 3, N), rng.uniform(60, 98, N),
        rng.uniform(970, 1030, N), rng.uniform(-1.8, 4, N), rng.uniform(0, 1, N))  # fmt: skip
KW = dict(zu=18, zt=15, zq=15, lat=72, lw_down=rng.uniform(180, 320, N),
          sw_down=rng.uniform(0, 200, N))  # fmt: skip


def test_versions_agree_without_wave_inputs():
    r35 = cs.ice_fluxes(*ARGS, **KW)
    r36 = cs.ice_fluxes(*ARGS, coare="3.6", **KW)
    # pycoare's two versions differ a little (measured: 0.8 W m-2, 1.2e-4 N m-2 at most)
    np.testing.assert_allclose(r36.water.hsb, r35.water.hsb, atol=1.0)
    np.testing.assert_allclose(r36.water.hlb, r35.water.hlb, atol=1.0)
    np.testing.assert_allclose(r36.water.tau, r35.water.tau, atol=3e-4)
    for name in ("tau", "hsb", "hlb"):  # the ice tile does not depend on the COARE version
        np.testing.assert_array_equal(getattr(r36.ice, name), getattr(r35.ice, name))
    assert r36.meta["coare"] == "3.6" and r35.meta["coare"] == "3.5"


@pytest.mark.parametrize("version", ["3.5", "3.6"])
def test_wave_inputs_pass_through(version):
    plain = cs.ice_fluxes(*ARGS, coare=version, **KW)
    waves = cs.ice_fluxes(*ARGS, coare=version, cp=8.0, sigH=np.full(N, 2.0), **KW)
    assert np.median(np.abs(waves.water.tau / plain.water.tau - 1.0)) > 0.02
    np.testing.assert_array_equal(waves.ice.tau, plain.ice.tau)
    # the form drag uses the open-water roughness, so it responds too
    assert np.nanmax(np.abs(waves.blended.tau_form - plain.blended.tau_form)) > 0


def test_salinity_is_used_by_3_6_only_and_bad_version_is_rejected():
    a = cs.ice_fluxes(*ARGS, coare="3.6", **KW)
    b = cs.ice_fluxes(*ARGS, coare="3.6", salinity=30.0, **KW)
    assert np.max(np.abs(a.water.hlb - b.water.hlb)) > 0.1
    with pytest.raises(ValueError, match="version must be"):
        cs.ice_fluxes(*ARGS, coare="4.0", **KW)
