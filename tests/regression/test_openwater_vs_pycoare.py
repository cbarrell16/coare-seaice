"""openwater.open_water equals plain pycoare, and its extra quantities are consistent."""

import numpy as np
import pytest
from pycoare import coare_35, coare_36

from coare_seaice import thermo
from coare_seaice.openwater import open_water

# pycoare's COARE 3.6 raises a harmless RuntimeWarning for sea temperatures below 1 C
# ((ts - 1) ** 0.82 of a negative number, of which it keeps the real part).
pytestmark = pytest.mark.filterwarnings("ignore:invalid value encountered in power")


def _forcing():
    rng = np.random.default_rng(7)
    n = 200
    return dict(
        u=rng.uniform(0.5, 25, n), t=rng.uniform(-30, 5, n), rh=rng.uniform(50, 100, n),
        ts=rng.uniform(-1.9, 6, n), p=rng.uniform(960, 1040, n), lat=rng.uniform(60, 82, n),
        rs=rng.uniform(0, 200, n), rl=rng.uniform(150, 320, n),
    )  # fmt: skip


@pytest.mark.parametrize("version,cls", [("3.5", coare_35), ("3.6", coare_36)])
@pytest.mark.parametrize("heights", [(36.0, 22.0, 22.0), (10.0, 10.0, 10.0)])
def test_fluxes_equal_plain_pycoare(version, cls, heights):
    f = _forcing()
    zu, zt, zq = heights
    ow = open_water(**f, zu=zu, zt=zt, zq=zq, version=version)
    c = cls(f["u"], t=f["t"], rh=f["rh"].copy(), ts=f["ts"], p=f["p"], lat=f["lat"],
            rs=f["rs"], rl=f["rl"], zu=zu, zt=zt, zq=zq)  # fmt: skip
    np.testing.assert_array_equal(ow.tau, c.fluxes.tau)
    np.testing.assert_array_equal(ow.hsb, c.fluxes.hsb)
    np.testing.assert_array_equal(ow.hlb, c.fluxes.hlb)
    np.testing.assert_array_equal(ow.usr, c.velocities.usr)
    np.testing.assert_array_equal(ow.zo, c.stability_parameters.zo)
    np.testing.assert_array_equal(ow.u_n_rf, c.velocities.u_n_rf)


@pytest.mark.parametrize("version", ["3.5", "3.6"])
def test_recomputed_quantities_consistent_with_pycoare(version):
    f = _forcing()
    ow = open_water(**f, zu=36.0, zt=22.0, zq=22.0, version=version)
    # tau = rhoa u*^2 / gf in pycoare, so this recovers pycoare's internal air density.
    np.testing.assert_allclose(ow.tau * ow.gust_factor / ow.usr**2, ow.rhoa, rtol=1e-13)
    # hlb = -rhoa Le u* q* gives pycoare's latent heat.
    np.testing.assert_allclose(-ow.hlb / (ow.rhoa * ow.usr * ow.qsr), ow.le, rtol=1e-12)
    # pycoare's own util (not used by the package) for the humidities.
    from pycoare import util

    np.testing.assert_allclose(ow.q_air, util.qair(f["t"], f["p"], f["rh"].copy()) / 1000,
                               rtol=1e-14)  # fmt: skip
    np.testing.assert_allclose(ow.q_sea, util.qsea(f["ts"], f["p"]) / 1000, rtol=1e-14)


def test_units_kg_per_kg():
    f = _forcing()
    ow = open_water(**f, zu=10.0, zt=10.0, zq=10.0)
    assert np.all(ow.q_air < 0.02) and np.all(ow.q_rf < 0.02) and np.all(np.abs(ow.qsr) < 1e-3)
    # dqer = wetc dter with wetc = 0.622 Le qs / (R (ts + 273.16)^2), all in kg/kg
    wetc = 0.622 * ow.le * ow.q_sea / (thermo.RGAS * (f["ts"] + thermo.TDK) ** 2)
    np.testing.assert_allclose(ow.dqer, wetc * ow.dter, rtol=1e-12)


def test_inputs_not_modified():
    f = _forcing()
    before = {k: v.copy() for k, v in f.items()}
    open_water(**f)
    for k in f:
        np.testing.assert_array_equal(f[k], before[k])


def test_bad_version():
    with pytest.raises(ValueError):
        open_water(5.0, -5.0, 80.0, 0.0, version="3.0")
