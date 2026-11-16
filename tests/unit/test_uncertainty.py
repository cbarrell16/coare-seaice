"""Task U1: perturbed-input ensembles."""

import numpy as np
import pytest

from coare_seaice.uncertainty import ensemble, log_uniform, normal, uniform

N = 8
INPUTS = dict(u=np.linspace(3, 15, N), t=np.full(N, -14.0), rh=np.full(N, 82.0), p=1000.0,
              sst=-1.5, sic=np.linspace(0.1, 0.95, N))  # fmt: skip
KW = dict(zu=20, zt=15, zq=15, lat=75, lw_down=220.0, sw_down=0.0)


def test_zero_perturbation_returns_the_central_value():
    ens = ensemble(INPUTS, dict(t=normal(0.0), u=normal(0.0, relative=True)), n=5, **KW)
    for name, pct in ens.percentiles.items():
        central = getattr(ens.central.blended, name.split(".")[1])
        for row in pct:
            np.testing.assert_array_equal(row, central)


def test_reproducible_with_a_seed():
    pert = dict(t=normal(0.3), sic=normal(0.1, per_sample=True), z0_ice=log_uniform(1e-4, 1e-2))
    a = ensemble(INPUTS, pert, n=20, seed=4, z0_ice=5e-4, **KW)
    b = ensemble(INPUTS, pert, n=20, seed=4, z0_ice=5e-4, **KW)
    c = ensemble(INPUTS, pert, n=20, seed=5, z0_ice=5e-4, **KW)
    np.testing.assert_array_equal(a.percentiles["blended.hsb"], b.percentiles["blended.hsb"])
    assert not np.array_equal(a.percentiles["blended.hsb"], c.percentiles["blended.hsb"])


def test_spread_brackets_the_central_value_and_grows_with_the_error():
    small = ensemble(INPUTS, dict(t=normal(0.1)), n=60, seed=1, **KW)
    large = ensemble(INPUTS, dict(t=normal(1.0)), n=60, seed=1, **KW)
    lo, mid, hi = large.percentiles["blended.hsb"]
    assert np.all(lo < large.central.blended.hsb) and np.all(large.central.blended.hsb < hi)
    width = lambda e: e.percentiles["blended.hsb"][2] - e.percentiles["blended.hsb"][0]  # noqa: E731
    assert np.all(width(large) > 5 * width(small))
    assert mid.shape == (N,)


def test_settings_other_outputs_and_members():
    ens = ensemble(INPUTS, dict(z0_ice=uniform(1e-4, 1e-2)), n=15, seed=2, z0_ice=5e-4,
                   outputs=("tiles.thick_ice.tau", "blended.tau"), levels=(25, 75),
                   keep_members=True, **KW)  # fmt: skip
    assert ens.percentiles["tiles.thick_ice.tau"].shape == (2, N)
    assert ens.members["blended.tau"].shape == (15, N)
    assert np.all(np.ptp(ens.members["tiles.thick_ice.tau"], axis=0) > 0)


def test_perturbed_ice_concentration_stays_a_fraction_and_bad_names_are_rejected():
    ens = ensemble(INPUTS, dict(sic=normal(2.0, per_sample=True)), n=10, seed=0, **KW)
    assert np.isfinite(ens.percentiles["blended.hsb"]).all()
    with pytest.raises(ValueError, match="nothing to perturb"):
        ensemble(INPUTS, dict(z0_ice=uniform(1e-4, 1e-2)), n=2, **KW)
    with pytest.raises(ValueError, match="missing"):
        ensemble(dict(u=1.0), {}, n=2, **KW)
