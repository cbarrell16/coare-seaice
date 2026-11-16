"""Task T7: roughness options."""

import numpy as np
import pytest

import coare_seaice as cs
from coare_seaice import thermo
from coare_seaice.roughness import (
    ICE_ROUGHNESS,
    cdn10_from_z0,
    resolve_ice_roughness,
    z0_andreas2010_winter,
    z0_from_cdn10,
)

N = 6
ARGS = (np.linspace(3.0, 18.0, N), np.full(N, -15.0), np.full(N, 80.0), np.full(N, 1000.0),
        np.full(N, -1.5), np.full(N, 0.7))  # fmt: skip
KW = dict(zu=10, zt=10, zq=10, lat=75, lw_down=220.0, sw_down=0.0)


def test_drag_coefficient_and_roughness_length_conversions():
    assert np.isclose(cdn10_from_z0(5e-4), 1.631e-3, rtol=1e-3)  # audit value
    assert np.isclose(z0_from_cdn10(cdn10_from_z0(3e-3)), 3e-3)
    # Renfrew et al. (2019) Table 1: z0_ice = 3e-3 m is "equivalent to CDN10_ice = 2.4e-3"
    assert np.isclose(cdn10_from_z0(3e-3), 2.4e-3, rtol=0.02)


def test_andreas_2010_eq_5_1():
    nu = 1.3e-5
    # hand values: 0.135 nu / u* + 2.30e-4 tanh^3(13 u*)
    np.testing.assert_allclose(z0_andreas2010_winter(0.05, nu),
                               0.135 * nu / 0.05 + 2.30e-4 * np.tanh(0.65) ** 3)  # fmt: skip
    assert np.isclose(z0_andreas2010_winter(0.6, nu), 2.30e-4, rtol=0.02)  # plateau
    assert np.isclose(cdn10_from_z0(2.30e-4), 1.4e-3, rtol=0.02)


def test_named_settings():
    assert resolve_ice_roughness("metum_smooth") == 0.5e-3
    assert np.isclose(cdn10_from_z0(resolve_ice_roughness("iceland_sea")), 3.4e-3)
    assert resolve_ice_roughness("andreas2010_winter") is z0_andreas2010_winter
    assert all("source" in v for v in ICE_ROUGHNESS.values())
    with pytest.raises(ValueError, match="unknown ice roughness"):
        resolve_ice_roughness("bumpy")
    with pytest.raises(ValueError, match="positive"):
        resolve_ice_roughness(-1.0)


def test_array_equals_loop_over_constants():
    z0 = np.geomspace(1e-4, 1e-2, N)
    res = cs.ice_fluxes(*ARGS, z0_ice=z0, **KW)
    for i in range(N):
        one = cs.ice_fluxes(*(a[i] for a in ARGS), z0_ice=float(z0[i]), **KW)
        np.testing.assert_allclose(res.ice.tau[i], one.ice.tau[0], rtol=1e-12)
        np.testing.assert_allclose(res.blended.hsb[i], one.blended.hsb[0], rtol=1e-12)
    np.testing.assert_array_equal(res.ice.z0, z0)


def test_rougher_ice_gives_more_stress_and_named_equals_value():
    smooth = cs.ice_fluxes(*ARGS, z0_ice="fram_strait_large_flat_floes", **KW)
    rough = cs.ice_fluxes(*ARGS, z0_ice="iceland_sea", **KW)
    assert np.all(rough.ice.tau > smooth.ice.tau)
    by_value = cs.ice_fluxes(*ARGS, z0_ice=float(z0_from_cdn10(3.4e-3)), **KW)
    np.testing.assert_array_equal(rough.ice.tau, by_value.ice.tau)


def test_friction_velocity_dependent_roughness_is_self_consistent():
    res = cs.ice_fluxes(*ARGS, z0_ice="andreas2010_winter", **KW)
    nu = thermo.air_viscosity(-15.0)
    np.testing.assert_allclose(res.ice.z0, z0_andreas2010_winter(res.ice.usr, nu), rtol=2e-3)
    assert not res.ice.qc["not_converged"].any()
    mine = cs.ice_fluxes(*ARGS, z0_ice=lambda usr, nu: 1e-3 + 0 * usr, **KW)  # a user function
    const = cs.ice_fluxes(*ARGS, z0_ice=1e-3, **KW)
    np.testing.assert_allclose(mine.ice.hsb, const.ice.hsb, rtol=1e-12)


def test_prescribed_scalar_roughness_on_an_ice_tile():
    tiles = cs.two_tile(z0=0.5e-3, z0t=0.1e-3, z0q=0.1e-3)  # Renfrew et al. (2019) Table 1
    res = cs.ice_fluxes(*ARGS, tiles=tiles, **KW)
    np.testing.assert_array_equal(res.ice.z0t, 0.1e-3)
    k = 0.4
    np.testing.assert_allclose(
        res.ice.chn10, k**2 / (np.log(10 / 0.5e-3) * np.log(10 / 0.1e-3)), rtol=1e-12
    )
    with pytest.raises(ValueError, match="both z0t and z0q"):
        cs.two_tile(z0t=1e-4)


def test_water_tile_with_prescribed_roughness():
    # Given COARE's own roughness lengths, the generic solver is close to COARE.
    ref = cs.ice_fluxes(*ARGS, jcool=0, **KW)
    w = ref.water
    tiles = cs.two_tile(water=dict(z0=w.z0, z0t=w.z0t, z0q=w.z0q))
    res = cs.ice_fluxes(*ARGS, tiles=tiles, jcool=0, **KW)
    # (COARE adds gustiness to the wind over water: 12 % at 3 m/s here, 2 % at 9 m/s.)
    strong = ARGS[0] >= 9.0
    for name in ("hsb", "hlb", "tau"):
        np.testing.assert_allclose(getattr(res.water, name)[strong], getattr(w, name)[strong],
                                   rtol=0.02, err_msg=name)  # fmt: skip
        assert np.all(getattr(res.water, name) <= getattr(w, name))
    np.testing.assert_allclose(res.blended.tau_form[strong], ref.blended.tau_form[strong],
                               rtol=0.02)  # fmt: skip
    np.testing.assert_array_equal(res.ice.hsb, ref.ice.hsb)
    # halving the scalar roughness lowers the heat flux; this is what anchoring does
    low = cs.two_tile(water=dict(z0=w.z0, z0t=w.z0t / 100, z0q=w.z0q / 100))
    assert np.all(cs.ice_fluxes(*ARGS, tiles=low, jcool=0, **KW).water.hsb < res.water.hsb)
    with pytest.raises(ValueError, match="needs z0, z0t and z0q"):
        cs.two_tile(water=dict(z0=1e-4))
