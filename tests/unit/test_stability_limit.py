"""Task F1: the named stability limit and the bulk-Richardson-number flag."""

import numpy as np

from coare_seaice import thermo
from coare_seaice.legacy import coare35vn_ice
from coare_seaice.most import (
    RIB_CRITICAL,
    SurfaceLayerOptions,
    bulk_richardson,
    ice_surface_layer,
)


def _solve(u, t, t_ice, options=None):
    u, t = np.atleast_1d(u).astype(float), np.atleast_1d(t).astype(float)
    n = u.size
    P = np.full(n, 1000.0)
    q = thermo.specific_humidity(t, P, np.full(n, 80.0))
    kw = dict(
        u=u, t=t, Q_air=q, t_ice=np.full(n, float(t_ice)), z0_ice=5e-4, P=P,
        grav=np.full(n, 9.82), zu=np.full(n, 10.0), zt=np.full(n, 10.0), zq=np.full(n, 10.0),
        rhoa=thermo.air_density(t, P, q), cpa=thermo.CPA, von=thermo.VON, fdg=thermo.FDG,
        ta=t + thermo.TDK, ice_scheme="blended_a87", scalar_ratio=0.2,
    )  # fmt: skip
    if options is not None:
        kw["options"] = options
    return ice_surface_layer(**kw)


def test_rib_formula_grachev_eq8():
    # g z / theta_v * d_theta / U^2 for dry air: 9.8 * 10 / 250 * 2 / 16 = 0.049
    rib = bulk_richardson(252.0, 250.0, 0.0, 0.0, 4.0, 10.0, 250.0, 9.8)
    assert np.isclose(rib, 0.049)


def test_flags_on_both_sides_of_the_limit():
    # windy and near neutral / calm and strongly stable
    fl = _solve(u=[10.0, 1.0], t=[-10.0, -10.0], t_ice=-20.0)
    assert fl["Rib"][0] < RIB_CRITICAL < fl["Rib"][1]
    assert list(fl["qc_rib"]) == [False, True]
    assert list(fl["zeta_limited"]) == [False, True]
    assert abs(fl["zeta"][0]) < 10.0 and np.isclose(fl["zeta"][1], 10.0)


def test_unstable_side_is_not_flagged_by_rib():
    fl = _solve(u=[5.0], t=[-25.0], t_ice=-5.0)
    assert fl["Rib"][0] < 0 and not fl["qc_rib"][0] and not fl["zeta_limited"][0]


def test_limit_is_a_parameter_and_default_is_v2():
    a = _solve(u=[1.0], t=[-10.0], t_ice=-20.0)
    b = _solve(u=[1.0], t=[-10.0], t_ice=-20.0, options=SurfaceLayerOptions(zeta_max=10.0))
    c = _solve(u=[1.0], t=[-10.0], t_ice=-20.0, options=SurfaceLayerOptions(zeta_max=2.0))
    assert a["hsb"][0] == b["hsb"][0]
    assert np.isclose(c["zeta"][0], 2.0) and c["hsb"][0] != a["hsb"][0]


def test_legacy_returns_diagnostics(inputs):
    kw = dict(P=inputs["P"], lat=inputs["lat"], Rs=inputs["Rs"], Rl=inputs["Rl"], zu=36,
              zt=22, zq=22, ice_conc=inputs["sic"], t_ice_method="air_freezing",
              verbose=False)  # fmt: skip
    args = (inputs["u"], inputs["t"], inputs["rh"], inputs["ts"])
    A = coare35vn_ice(*args, **kw)
    B, diag = coare35vn_ice(*args, return_diagnostics=True, **kw)
    np.testing.assert_array_equal(A, B)
    assert set(diag) >= {"zeta", "Rib", "zeta_limited", "qc_rib"}
    assert diag["qc_rib"].dtype == bool and diag["qc_rib"].shape == (A.shape[0],)


def test_viscosity_option_f2():
    # Colder air is less viscous, so R* is larger and z0T smaller than with 1.5e-5.
    from coare_seaice.most import V2_OPTIONS

    new = _solve(u=[8.0], t=[-20.0], t_ice=-15.0)
    old = _solve(u=[8.0], t=[-20.0], t_ice=-15.0, options=V2_OPTIONS)
    nu = thermo.air_viscosity(-20.0)
    assert np.isclose(new["R_star"][0] / old["R_star"][0], 1.5e-5 / nu, rtol=2e-2)
    assert new["z0T"][0] < old["z0T"][0]
    assert 0.0 < abs(new["hsb"][0] / old["hsb"][0] - 1.0) < 0.05
