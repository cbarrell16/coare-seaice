"""thermo: each function equals its v2 counterpart (meteo.py / ice_surface_temperature.py)."""

import numpy as np
import pytest

from coare_seaice import thermo
from conftest import assert_port, load_ref


@pytest.fixture(scope="module")
def ref():
    return load_ref("thermo_stability")


def test_esat_water(ref):
    assert_port(thermo.esat_water(ref["t"], ref["p"]), ref["esat_water"])


def test_qsat_sea(ref):
    assert_port(thermo.qsat_sea(ref["t"], ref["p"]), ref["qsat_sea"])


def test_qsat_sea_salinity_35_equals_default(ref):
    assert_port(thermo.qsat_sea(ref["t"], ref["p"], 35.0), ref["qsat_sea"])


def test_specific_humidity_and_vapour_pressure(ref):
    assert_port(thermo.specific_humidity(ref["t"], ref["p"], ref["rh"]), ref["q_air"])
    assert_port(thermo.vapour_pressure(ref["t"], ref["p"], ref["rh"]), ref["e_air"])


def test_gravity(ref):
    assert_port(thermo.gravity(ref["lat"]), ref["gravity"])


def test_ice_saturation(ref):
    assert_port(thermo.esat_ice(ref["t"], ref["p"]), ref["esat_ice"])
    assert_port(thermo.qsat_ice(ref["t"], ref["p"]), ref["qsat_ice"])
    assert_port(thermo.dqsat_ice_dt(ref["t"], ref["p"]), ref["dqsat_ice_dt"])


def test_air_density_and_latent_heat_equal_v2_columns(inputs):
    # v2 output columns 25 (Le) and 26 (rhoa) are computed from the inputs alone.
    A = load_ref("legacy_openwater_z36")["A"]
    q = thermo.specific_humidity(inputs["t"], inputs["P"], inputs["rh"])
    assert_port(thermo.air_density(inputs["t"], inputs["P"], q), A[:, 26])
    assert_port(thermo.latent_heat_vaporisation(inputs["ts"]), A[:, 25])


def test_inputs_not_modified():
    rh = np.array([80.0, 90.0])
    thermo.specific_humidity(np.array([-5.0, 0.0]), np.array([1000.0, 1000.0]), rh)
    assert rh.tolist() == [80.0, 90.0]


def test_air_viscosity_matches_pycoare_and_audit_values():
    # Task F2.
    from pycoare import coare_35

    from coare_seaice.thermo import air_viscosity

    t = np.array([-30.0, -20.0, -10.0, 0.0, 5.0])
    np.testing.assert_allclose(air_viscosity(0.0), 1.326e-5, rtol=1e-12)
    np.testing.assert_allclose(air_viscosity(-10.0), 1.24e-5, rtol=2e-3)
    np.testing.assert_allclose(air_viscosity(-20.0), 1.157e-5, rtol=2e-3)
    # pycoare keeps its value on a private object, so this part is skipped if that moves.
    c = coare_35(np.full(5, 5.0), t=t.copy(), rh=np.full(5, 80.0), ts=np.full(5, 1.0))
    visa = getattr(getattr(c, "_bulk_loop_inputs", None), "visa", None)
    if visa is not None:
        np.testing.assert_allclose(air_viscosity(t), visa, rtol=1e-12)
