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
    col = IceColumn(**CONFIGS[name])
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
    res = IceColumn().run(ref["run__in_time_s"], *(f[k] for k in FORCING),
                          ref["run__in_Ts_init"], sic=ref["run__in_sic"])  # fmt: skip
    for k in ("Ts", "F_cond", "melt"):
        assert_port(res[k], ref[f"run__{k}"])
    np.testing.assert_array_equal(res["reinit"], ref["run__reinit"])
    res = IceColumn(h_snow=0.0).run(ref["run__in_time_s"], *(f[k] for k in FORCING),
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
