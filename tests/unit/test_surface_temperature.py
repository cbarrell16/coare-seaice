"""surface_temperature.simple / seb and the LW fallback equal their v2 counterparts."""

import numpy as np
import pytest

from coare_seaice.radiation import lw_down_konig_langlo
from coare_seaice.surface_temperature import (
    air_capped,
    air_freezing,
    effective_conductance,
    empirical,
    seb_terms,
    solve_seb_steady,
    surface_albedo,
)
from conftest import assert_port, load_ref

ARRAY_ARGS = ("u", "t", "Q_air", "P", "grav", "zu", "zt", "zq", "rhoa", "ta")


@pytest.fixture(scope="module")
def ref():
    return load_ref("thermo_stability")


@pytest.fixture(scope="module")
def seb():
    return load_ref("seb")


def test_simple_methods(ref):
    assert_port(air_freezing(ref["t"]), ref["air_freezing"])
    assert_port(air_capped(ref["t"]), ref["air_capped"])
    assert_port(empirical(ref["t"]), ref["empirical"])
    assert_port(empirical(ref["t"], offset=3.5), ref["empirical_3p5"])


def test_empirical_func():
    t = np.array([-10.0, -1.0])
    assert_port(empirical(t, func=lambda x, k: k * np.ones_like(x), k=1.5), [-8.5, 0.0])


def test_albedo_conductance(ref):
    assert_port(surface_albedo(ref["h_snow"]), ref["surface_albedo"])
    assert_port(effective_conductance(ref["h_ice"], ref["h_snow"]), ref["effective_conductance"])
    assert_port(
        effective_conductance(ref["h_ice"], ref["h_snow"], 1.9, 0.25),
        ref["effective_conductance_k"],
    )


def test_lw_down_konig_langlo(ref):
    assert_port(lw_down_konig_langlo(ref["t"], ref["cloud"]), ref["lw_down_konig_langlo"])


def test_seb_terms(seb, inputs):
    lw, sw, fc = seb_terms(inputs["t_ice"], inputs["Rl"], inputs["Rs"], seb["albedo"],
                           seb["G"], -1.8)  # fmt: skip
    assert_port(lw, seb["lw_net"])
    assert_port(sw, seb["sw_net"])
    assert_port(fc, seb["F_cond"])


def test_solve_seb_steady(seb, inputs):
    args = {k: seb[f"in_{k}"] for k in ARRAY_ARGS}
    args.update(z0_ice=5e-4, cpa=1004.67, von=0.4, fdg=1.0, ice_scheme="blended_a87",
                scalar_ratio=0.2, nits=10)  # fmt: skip
    Ts, res = solve_seb_steady(args, inputs["Rl"], inputs["Rs"], seb["albedo"], seb["G"],
                               -1.8, inputs["t"])  # fmt: skip
    assert_port(Ts, seb["Ts"])
    assert_port(res, seb["res"])
