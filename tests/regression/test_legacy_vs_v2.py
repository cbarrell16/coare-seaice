"""legacy.coare35vn_ice against stored v2 outputs (tests/fixtures/v2_reference).

The comparison runs with ``_v2_compat=True``, which switches off the audit fixes that
change results (release 0.3.0), so it still proves that the port equals v2 and that
those fixes are the only differences. The default behaviour is covered by the unit tests
of each fix.

Open-water columns 0-37 come from pycoare, which differs slightly from the NOAA COARE 3.5
loop in v2 (WGS84 gravity; exact 1/3 and 2/3 in the gustiness and stability functions
where the NOAA code has 0.333 / 0.3333 / 0.6667 and the constants 14.28, 8.525). So they
are compared with rtol = 2e-3 plus an absolute tolerance per column, chosen to be well
below what matters physically or can be measured. Ice columns 38-55 are compared with
rtol = 1e-3; only the blends and form drag depend on the open-water part.

Columns whose rtol alone is not met, and why (measured on the fixture, see report):
    16 dter, 17 dqer  cool-skin depression: relative differences up to 8e-3, at most
                      1.4e-4 K / 4.6e-8 kg/kg; the cool skin amplifies the small flux
                      differences.
    20 Trf            2 m temperature: relative difference up to 5e-2 where Trf is near
                      0 C; absolute difference at most 1.5e-3 K.
    48 hsb_blend,     relative difference up to 4e-2 where the open-water and ice
    49 hlb_blend      contributions nearly cancel; absolute at most 0.09 W m-2.
    54 SEB_res        for 'seb_steady' this is the converged residual, ~1e-10 W m-2 (zero
                      to rounding); different numpy builds differ there by ~1e-13 W m-2.
                      It is identical to v2 with the numpy used to make the fixtures.
"""

from pathlib import Path

import numpy as np
import pytest

from coare_seaice import thermo
from coare_seaice.legacy import OUTPUT_COLUMNS, _ice_columns, coare35vn_ice
from coare_seaice.most import V2_OPTIONS
from conftest import REF_DIR, assert_port, load_ref, resolve_kwargs

RTOL_OW = 2e-3
RTOL_ICE = 1e-3

# Absolute tolerance per column (units of that column).
ATOL = {
    "usr": 1e-5, "tau": 1e-5, "hsb": 0.05, "hlb": 0.05, "hbb": 0.05, "hsbb": 0.05,
    "hlwebb": 1e-3, "tsr": 1e-5, "qsr": 1e-6, "zot": 1e-9, "zoq": 1e-9,
    "Cd": 1e-7, "Ch": 1e-7, "Ce": 1e-7, "L": 1.0, "zet": 1e-4,
    "dter": 5e-4, "dqer": 2e-7, "tkt": 1e-6,
    "Urf": 1e-3, "Trf": 5e-3, "Qrf": 1e-3, "RHrf": 0.05, "UrfN": 1e-3, "Rnl": 1e-3,
    "Le": 0.0, "rhoa": 0.0, "UN": 1e-3, "U10": 1e-3, "U10N": 1e-3,
    "Cdn_10": 1e-5, "Chn_10": 1e-5, "Cen_10": 1e-5, "RF": 1e-3, "Evap": 1e-5,
    "Qs": 1e-4, "Q10": 1e-3, "RH10": 0.05,
    "tau_blend": 1e-5, "hsb_blend": 0.05, "hlb_blend": 0.05, "SEB_res": 1e-6,
    "tau_form": 1e-6,
}  # fmt: skip

RUNS = sorted(p.stem[len("legacy_") :] for p in Path(REF_DIR).glob("legacy_*.npz"))


def _v2_array(name):
    d = load_ref(f"legacy_{name}")
    if "A" in d:
        return d["A"], d["kwargs"]
    ow = load_ref(f"legacy_openwater_{d['openwater']}")["A"]
    return np.column_stack([ow, d["ice"]]), d["kwargs"]


def _run(inputs, kwargs_json):
    kw = resolve_kwargs(kwargs_json, inputs)
    return coare35vn_ice(
        inputs["u"], inputs["t"], inputs["rh"], inputs["ts"], _v2_compat=True, **kw
    )


@pytest.mark.parametrize("name", RUNS)
def test_legacy_matches_v2(inputs, name):
    A, kw = _v2_array(name)
    B = _run(inputs, kw)
    assert B.shape == A.shape
    for i in range(A.shape[1]):
        col = OUTPUT_COLUMNS[i]
        rtol = RTOL_OW if i < 38 else RTOL_ICE
        np.testing.assert_allclose(
            B[:, i], A[:, i], rtol=rtol, atol=ATOL.get(col, 0.0), equal_nan=True,
            err_msg=f"run {name}, column {i} ({col})",
        )  # fmt: skip


def test_runs_cover_every_method():
    import json

    methods = {json.loads(str(load_ref(f"legacy_{r}")["kwargs"])).get("t_ice_method") for r in RUNS}
    assert set(methods) >= {
        "air_freezing", "air_capped", "empirical", "prescribed", "seb_steady", "diffusion",
    }  # fmt: skip


@pytest.mark.parametrize("name", [r for r in RUNS if not r.startswith("openwater")])
def test_ice_columns_exact_given_v2_open_water(inputs, name):
    """With v2's own open-water values fed in, the ice port reproduces v2 to rounding."""
    A, kw_json = _v2_array(name)
    kw = resolve_kwargs(kw_json, inputs)
    n = A.shape[0]
    t, rh, P = inputs["t"], inputs["rh"], kw["P"]
    h_ice = np.asarray(kw.get("h_ice", 1.0)) * np.ones(n)
    h_snow = np.asarray(kw.get("h_snow", 0.2)) * np.ones(n)
    sic = kw["ice_conc"]
    method = kw.get("t_ice_method", "seb_steady")
    if kw.get("t_ice") is not None:
        method = "prescribed"
    # v2's open-water roughness from its Cdn_10 column: Cdn_10 = 1e3 k^2 / ln(10/zo)^2.
    zo_v2 = 10.0 / np.exp(0.4 / np.sqrt(A[:, 30] / 1000.0))
    cols, _ = _ice_columns(
        u=inputs["u"], t=t, rh=rh, P=P, Rs=kw["Rs"], Rl=kw["Rl"],
        grav=thermo.gravity(kw["lat"]),
        zu=kw["zu"] * np.ones(n), zt=kw["zt"] * np.ones(n), zq=kw["zq"] * np.ones(n),
        q_air=thermo.specific_humidity(t, P, rh), rhoa=A[:, 26], zo_w=zo_v2,
        u10n_w=A[:, 29], tau_w=A[:, 1], hsb_w=A[:, 2], hlb_w=A[:, 3],
        A_ice=np.where(np.isfinite(sic), np.clip(sic, 0, 1), np.nan),
        z0_ice=5e-4, ice_scheme=kw.get("ice_scheme", "blended_a87"), scalar_ratio=0.2,
        t_ice_method=method, t_ice=kw.get("t_ice"),
        t_ice_fallback=kw.get("t_ice_fallback", "seb_steady"),
        h_ice=h_ice, h_snow=h_snow, k_snow=kw.get("k_snow", 0.31),
        albedo=np.where(h_snow > 0.02, 0.85, 0.65), t_base=-1.8,
        t_ice_offset=kw.get("t_ice_offset", 2.0), time_s=kw.get("time_s"),
        diffusion_kwargs=None, rh_reference=kw.get("rh_reference", "water"), verbose=False,
        options=V2_OPTIONS, v2_compat=True,
    )  # fmt: skip
    for j, col in enumerate(cols):
        name_j = OUTPUT_COLUMNS[38 + j]
        if name_j == "SEB_res":  # converged residual ~1e-10 W m-2 for the SEB methods
            np.testing.assert_allclose(col, A[:, 38 + j], rtol=1e-10, atol=1e-6, equal_nan=True)
        else:
            assert_port(col, A[:, 38 + j], rtol=1e-10)
