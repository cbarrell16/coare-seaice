"""most.ice_surface_layer equals v2 _ice_fluxes on stored inputs."""

import pytest

from coare_seaice.most import ice_surface_layer
from conftest import assert_port, load_ref

ARRAY_ARGS = ("u", "t", "Q_air", "P", "grav", "zu", "zt", "zq", "rhoa", "ta")
FLUX_KEYS = ("tau", "hsb", "hlb", "Cdn", "Chn", "Cen", "R_star", "z0T", "z0q", "h_s", "h_l",
             "Qs", "L")  # fmt: skip
CASES = {"a87": "blended_a87", "a87_z10": "blended_a87", "simple": "simple"}


@pytest.fixture(scope="module")
def ref():
    return load_ref("most")


@pytest.mark.parametrize("case", sorted(CASES))
def test_ice_surface_layer(ref, case):
    args = {k: ref[f"{case}__in_{k}"] for k in ARRAY_ARGS}
    fl = ice_surface_layer(
        t_ice=ref["t_ice"], z0_ice=float(ref[f"{case}__z0_ice"]), cpa=1004.67, von=0.4,
        fdg=1.0, ice_scheme=CASES[case], scalar_ratio=0.2, nits=10, **args,
    )  # fmt: skip
    for k in FLUX_KEYS:
        assert_port(fl[k], ref[f"{case}__{k}"])
