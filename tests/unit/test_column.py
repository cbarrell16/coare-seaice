"""surface_temperature.column.IceColumn equals v2 IceColumn step by step."""

import numpy as np
import pytest

from coare_seaice.surface_temperature import IceColumn
from conftest import assert_port, load_ref

CONFIGS = {
    "default": dict(),
    "nosnow_fresh": dict(h_ice=0.6, h_snow=0.0, n_ice=5, salinity_dependent=False),
    "thin_k": dict(h_ice=0.3, h_snow=0.05, n_ice=6, n_snow=2, k_snow=0.25, s_ice=8.0),
}
FORCING = ("theta_a", "q_a", "p", "lw_down", "sw_down", "h_s", "h_l")


@pytest.fixture(scope="module")
def ref():
    return load_ref("column")


@pytest.mark.parametrize("name", sorted(CONFIGS))
def test_steps(ref, name):
    col = IceColumn(v2_compat=True, **CONFIGS[name])
    assert_port(col.steady_profile(-15.0), ref[f"{name}__steady_profile"])
    col.initialise(-12.0)
    f = {k: ref[f"steps__in_{k}"] for k in FORCING}
    out = np.array([col.step(600.0, *(f[k][i] for k in FORCING)) for i in range(40)])
    assert_port(out[:, 0], ref[f"{name}__step_Ts"])
    assert_port(out[:, 1], ref[f"{name}__step_F_cond"])
    assert_port(out[:, 2], ref[f"{name}__step_melt"])
    assert_port(col.T, ref[f"{name}__step_T"][-1])
    col.initialise(-20.0)
    r = col.step(3 * 3600.0, *(f[k][5] for k in FORCING), albedo=0.7)
    assert_port(np.array(r), ref[f"{name}__step_albedo"])


@pytest.mark.parametrize("name", sorted(CONFIGS))
def test_reference_reaches_melt(ref, name):
    assert (ref[f"{name}__step_melt"] > 0).any()  # the melt branch is exercised


def test_run(ref):
    f = {k: ref[f"run__in_{k}"] for k in FORCING}
    res = IceColumn(v2_compat=True).run(ref["run__in_time_s"], *(f[k] for k in FORCING),
                          ref["run__in_Ts_init"], sic=ref["run__in_sic"])  # fmt: skip
    for k in ("Ts", "F_cond", "melt"):
        assert_port(res[k], ref[f"run__{k}"])
    np.testing.assert_array_equal(res["reinit"], ref["run__reinit"])
    res = IceColumn(h_snow=0.0, v2_compat=True).run(ref["run__in_time_s"], *(f[k] for k in FORCING),
                                    ref["run__in_Ts_init"], albedo=0.6)  # fmt: skip
    for k in ("Ts", "F_cond", "melt"):
        assert_port(res[k], ref[f"run_nosic__{k}"])
    np.testing.assert_array_equal(res["reinit"], ref["run_nosic__reinit"])


def test_conductive_flux_when_surface_is_warmer_than_base():
    # Task F4: surface at -1 C, base at -1.8 C. The flux is finite and downward (negative).
    col = IceColumn(h_ice=1.0, h_snow=0.2)
    col.initialise(-1.0)
    forcing = dict(theta_a=-1.0, q_a=3.4e-3, p=1000.0, lw_down=305.0, sw_down=0.0, h_s=1e4,
                   h_l=0.0)  # fmt: skip
    Ts, F_cond, _ = col.step(600.0, **forcing)
    assert abs(Ts + 1.0) < 0.05
    assert np.isfinite(F_cond) and -5.0 < F_cond < 0.0  # about G dT = 0.88 * (-0.8) W m-2


def test_legacy_diffusion_f_cond_is_the_column_flux(inputs):
    from coare_seaice.legacy import OUTPUT_COLUMNS, coare35vn_ice

    n = inputs["u"].size
    t = np.linspace(-6.0, 3.0, n)  # warms through the base temperature
    kw = dict(P=1000.0, lat=75.0, Rs=0.0, Rl=300.0, zu=10, zt=10, zq=10, ice_conc=np.ones(n),
              t_ice_method="diffusion", time_s=600.0 * np.arange(n), verbose=False)  # fmt: skip
    args = (np.full(n, 6.0), t, np.full(n, 90.0), np.full(n, -1.8))
    i_f, i_t = OUTPUT_COLUMNS.index("F_cond"), OUTPUT_COLUMNS.index("T_ice")
    new = coare35vn_ice(*args, **kw)
    old = coare35vn_ice(*args, _v2_compat=True, **kw)
    warm = new[:, i_t] > -1.8
    assert warm.sum() > 10
    assert np.abs(new[warm, i_f]).max() < 50.0  # finite and of physical size
    assert np.abs(old[warm, i_f]).max() > 1e3  # the v2 defect


PINNED = dict(theta_a=-15.0, q_a=1.0e-3, p=1000.0, lw_down=250.0, sw_down=0.0, h_s=1e7, h_l=0.0)


@pytest.mark.parametrize("name", sorted(CONFIGS))
def test_no_drift_after_initialisation(name):
    # Task F5: with the surface held at -15 C the initial profile is already steady.
    col = IceColumn(**CONFIGS[name])
    col.initialise(-15.0)
    T0, F0 = col.T.copy(), col.steady_flux(-15.0)
    for _ in range(144):
        _, F_cond, _ = col.step(600.0, **PINNED)
    assert np.max(np.abs(col.T - T0)) < 1e-4
    assert abs(F_cond - F0) < 1e-3 * abs(F0)


def test_v2_initial_profile_drifted():
    col = IceColumn(v2_compat=True)  # salinity-dependent column, fresh-ice initial profile
    col.initialise(-15.0)
    T0 = col.T.copy()
    for _ in range(144):
        col.step(600.0, **PINNED)
    assert np.max(np.abs(col.T - T0)) > 1e-2


def test_steady_profile_fresh_ice_matches_analytic_slab():
    col = IceColumn(h_ice=1.0, h_snow=0.2, salinity_dependent=False)
    G = 1.0 / (1.0 / 2.03 + 0.2 / 0.31)
    assert np.isclose(col.steady_flux(-20.0), G * (-1.8 + 20.0), rtol=1e-12)
    old = IceColumn(h_ice=1.0, h_snow=0.2, salinity_dependent=False, v2_compat=True)
    np.testing.assert_allclose(col.steady_profile(-20.0), old.steady_profile(-20.0), rtol=1e-12)


def test_melt_is_reported_only_for_a_melting_surface():
    # Task F5: melt is reset at each Picard pass, so a melt flux can only be reported
    # together with a surface at 0 C. (In v2 a value from an earlier pass could carry
    # over; no such case occurred in 20 000 random steps, so the effect is very rare.)
    rng = np.random.default_rng(5)
    melted = 0
    for _ in range(300):
        col = IceColumn(h_ice=rng.uniform(0.1, 2.0), h_snow=rng.choice([0.0, 0.05, 0.2]))
        col.initialise(rng.uniform(-3.0, -0.05))
        Ts, _, melt = col.step(
            600.0, theta_a=rng.uniform(-3.0, 6.0), q_a=rng.uniform(0.002, 0.006), p=1005.0,
            lw_down=rng.uniform(250.0, 340.0), sw_down=rng.uniform(0.0, 400.0),
            h_s=rng.uniform(5.0, 30.0), h_l=rng.uniform(1e4, 6e4),
        )  # fmt: skip
        assert melt >= 0.0 and (melt == 0.0 or Ts == 0.0)
        melted += melt > 0
    assert 20 < melted < 280  # both branches are exercised


def test_diffusion_with_thickness_array_warns_and_uses_median(inputs):
    from coare_seaice.legacy import coare35vn_ice

    n = 60
    kw = dict(P=1000.0, lat=75.0, Rs=0.0, Rl=230.0, zu=10, zt=10, zq=10, ice_conc=np.ones(n),
              t_ice_method="diffusion", time_s=600.0 * np.arange(n), verbose=False)  # fmt: skip
    args = (np.full(n, 6.0), np.full(n, -15.0), np.full(n, 85.0), np.full(n, -1.8))
    h = np.linspace(0.5, 1.5, n)
    with pytest.warns(UserWarning, match="one slab"):
        a = coare35vn_ice(*args, h_ice=h, **kw)
    b = coare35vn_ice(*args, h_ice=float(np.median(h)), **kw)
    np.testing.assert_allclose(a, b, rtol=1e-12, equal_nan=True)
