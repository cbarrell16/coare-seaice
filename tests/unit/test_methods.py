"""Task T2: surface temperature methods behind one interface."""

from dataclasses import dataclass

import numpy as np
import pytest

from coare_seaice import thermo
from coare_seaice.most import ice_surface_layer
from coare_seaice.surface_temperature import (
    SEB,
    AirCapped,
    AirFreezing,
    Empirical,
    Prescribed,
    SurfaceTemperature,
    TileContext,
    effective_conductance,
    method_from_name,
    simple,
    solve_seb_steady,
    surface_albedo,
)


@pytest.fixture(scope="module")
def ctx(inputs):
    n = inputs["u"].size
    t, P = inputs["t"], inputs["P"]
    q = thermo.specific_humidity(t, P, inputs["rh"])
    args = dict(
        u=inputs["u"], t=t, Q_air=q, z0_ice=5e-4, P=P, grav=thermo.gravity(inputs["lat"]),
        zu=np.full(n, 36.0), zt=np.full(n, 22.0), zq=np.full(n, 22.0),
        rhoa=thermo.air_density(t, P, q), cpa=thermo.CPA, von=thermo.VON, fdg=thermo.FDG,
        ta=t + thermo.TDK, ice_scheme="blended_a87", scalar_ratio=0.2,
    )  # fmt: skip
    h_snow = np.full(n, 0.2)
    return TileContext(
        t_air=t, zt=args["zt"], p=P, q_air=q, lw_down=inputs["Rl"], sw_down=inputs["Rs"],
        sic=inputs["sic"], flux=lambda Ts: ice_surface_layer(t_ice=Ts, **args),
        h_ice=np.full(n, 1.0), h_snow=h_snow, k_snow=0.31, albedo=surface_albedo(h_snow),
        t_base=-1.8,
    ), args  # fmt: skip


def test_simple_methods_match_the_functions(ctx):
    c, _ = ctx
    np.testing.assert_array_equal(AirFreezing().solve(c).t_surface, simple.air_freezing(c.t_air))
    np.testing.assert_array_equal(AirCapped().solve(c).t_surface, simple.air_capped(c.t_air))
    np.testing.assert_array_equal(
        Empirical(offset=3.0).solve(c).t_surface, simple.empirical(c.t_air, 3.0)
    )


def test_seb_matches_the_solver(ctx):
    c, args = ctx
    G = effective_conductance(c.h_ice, c.h_snow, 2.03, 0.31)
    Ts, res = solve_seb_steady(args, c.lw_down, c.sw_down, c.albedo, G, -1.8, c.t_air)
    out = SEB().solve(c)
    np.testing.assert_array_equal(out.t_surface, Ts)
    np.testing.assert_array_equal(out.seb_residual, res)


def test_prescribed_and_fallback(ctx):
    c, _ = ctx
    T = np.where(np.arange(c.n) % 3 == 0, np.nan, -12.0)
    plain = Prescribed(T).solve(c)
    assert np.isnan(plain.t_surface[0]) and not plain.filled.any()
    filled = Prescribed(T, fallback=AirFreezing()).solve(c)
    np.testing.assert_array_equal(filled.filled, np.isnan(T))
    np.testing.assert_array_equal(filled.t_surface[::3], simple.air_freezing(c.t_air)[::3])
    assert np.all(Prescribed(-5.0).solve(c).t_surface == -5.0)
    with pytest.raises(ValueError):
        Prescribed(np.zeros(3)).solve(c)


def test_names():
    assert isinstance(method_from_name("air_freezing"), AirFreezing)
    assert isinstance(method_from_name("seb_steady", h_ice=0.3), SEB)
    with pytest.raises(ValueError):
        method_from_name("nope")


def test_a_new_method_needs_no_change_elsewhere(ctx):
    @dataclass
    class TwoBelowAir:  # written by a user, outside the package
        energy_balance = False

        def solve(self, c):
            return SurfaceTemperature(np.minimum(c.t_air - 2.0, 0.0))

    c, _ = ctx
    out = Prescribed(np.full(c.n, np.nan), fallback=TwoBelowAir()).solve(c)
    np.testing.assert_array_equal(out.t_surface, np.minimum(c.t_air - 2.0, 0.0))
