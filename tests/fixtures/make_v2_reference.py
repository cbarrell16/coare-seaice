"""Generate the v2 reference fixtures in ``tests/fixtures/v2_reference/``.

Run once, by hand, with the v2 code (coare35vn_ice.py 2.0.0, ice_surface_temperature.py,
meteo.py, util.py) in a directory given on the command line:

    python tests/fixtures/make_v2_reference.py /path/to/v2 [output_dir]

It is not run by CI. It imports only numpy and the v2 modules (not coare_seaice), so the
references are independent of the code they test.

Contents (all synthetic, fixed seed):
    inputs.npz              300 polar cases: forcing, SIC (with exact 0, 1 and NaN),
                            prescribed T_ice (with NaN), thickness arrays, time axis
    thermo_stability.npz    saturation, humidity, gravity, psi functions, LW fallback,
                            albedo, conductance, simple T_ice methods
    roughness_form_drag.npz Andreas (1987) scalar roughness grid; form drag grid
    most.npz                ice surface-layer solver on cases x surface temperatures
    seb.npz                 seb_terms and the steady SEB solve
    column.npz              IceColumn: steady profile, step sequences, run()
    legacy_openwater_*.npz  v2 columns 0-37, per height set
    legacy_<run>.npz        v2 columns 38-55 (or all 56 when the open-water part has its
                            own inputs) plus the keyword arguments of the call as JSON
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np

SEED = 20261113
N = 300


def make_inputs(rng):
    u = rng.uniform(0.5, 25.0, N)
    t = rng.uniform(-30.0, 5.0, N)
    rh = rng.uniform(50.0, 100.0, N)
    ts = rng.uniform(-1.9, 6.0, N)
    P = rng.uniform(960.0, 1040.0, N)
    lat = rng.uniform(60.0, 82.0, N)
    sic = rng.uniform(0.0, 1.0, N)
    sic[0:15] = 0.0
    sic[15:30] = 1.0
    sic[30:36] = np.nan
    sic[36:42] = rng.uniform(0.0, 0.01, 6)  # below the form-drag threshold
    Rs = rng.uniform(0.0, 200.0, N)
    Rl = rng.uniform(150.0, 320.0, N)
    t_ice = np.minimum(t + rng.uniform(-2.0, 4.0, N), 0.0)
    t_ice[rng.choice(N, 40, replace=False)] = np.nan
    h_ice = rng.uniform(0.1, 3.0, N)
    h_snow = rng.uniform(0.0, 0.5, N)
    h_snow[::10] = 0.01  # bare-ice albedo
    dt = np.full(N, 600.0)
    dt[[50, 180]] = 4 * 3600.0  # gaps longer than max_gap_s
    time_s = np.cumsum(dt)
    rain = rng.uniform(0.0, 5.0, N)
    return dict(
        u=u, t=t, rh=rh, ts=ts, P=P, lat=lat, sic=sic, Rs=Rs, Rl=Rl, t_ice=t_ice,
        h_ice=h_ice, h_snow=h_snow, time_s=time_s, rain=rain,
    )  # fmt: skip


def thermo_stability(meteo, ist):
    rng = np.random.default_rng(SEED + 1)
    t = np.linspace(-45.0, 10.0, 221)
    p = rng.uniform(960.0, 1040.0, t.size)
    rh = rng.uniform(40.0, 100.0, t.size)
    qa, em = meteo.qair(t, p, rh)
    lat = np.linspace(-90.0, 90.0, 181)
    pos = np.logspace(-9, 3, 300)
    zeta = np.concatenate([-pos[::-1], [0.0, np.nan], pos])
    cloud = np.linspace(-0.2, 1.2, t.size)
    h = np.linspace(-0.1, 3.0, 63)
    hs = np.linspace(-0.05, 0.6, 63)
    return dict(
        t=t, p=p, rh=rh, lat=lat, zeta=zeta, cloud=cloud, h_ice=h, h_snow=hs,
        esat_water=meteo.qsat(t, p),
        qsat_sea=meteo.qsea(t, p) / 1000.0,
        q_air=qa / 1000.0,
        e_air=em,
        gravity=meteo.grv(lat),
        esat_ice=ist.esat_ice(t, p),
        qsat_ice=ist.qsat_ice(t, p),
        dqsat_ice_dt=ist.dqsat_ice_dt(t, p),
        psiu_26=meteo.psiu_26(zeta),
        psit_26=meteo.psit_26(zeta),
        lw_down_konig_langlo=ist.lw_down_konig_langlo(t, cloud),
        surface_albedo=ist.surface_albedo(hs),
        effective_conductance=ist.effective_conductance(h, hs),
        effective_conductance_k=ist.effective_conductance(h, hs, 1.9, 0.25),
        air_freezing=ist.t_ice_air_freezing(t),
        air_capped=ist.t_ice_air_capped(t),
        empirical=ist.t_ice_empirical(t),
        empirical_3p5=ist.t_ice_empirical(t, offset=3.5),
    )  # fmt: skip


def roughness_form_drag(v2):
    z0 = np.array([1e-5, 1e-4, 5e-4, 1e-3, 1e-2])
    usr = np.logspace(-3, 0.5, 60)
    Z0, US = (a.ravel() for a in np.meshgrid(z0, usr, indexing="ij"))
    Z0 = np.append(Z0, [5e-4, np.nan])
    US = np.append(US, [np.nan, 0.3])
    z0T, z0q, R = v2._a87_scalar_roughness(Z0, US)
    A = np.array([np.nan, -0.1, 0.0, 0.005, 0.0099, 0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9,
                  0.99, 0.999999, 1.0, 1.2])  # fmt: skip
    z0w = np.array([1e-5, 2e-4, 1e-3, 0.6])
    AA, ZW = (a.ravel() for a in np.meshgrid(A, z0w, indexing="ij"))
    return dict(
        a87_z0=Z0, a87_usr=US, a87_z0T=z0T, a87_z0q=z0q, a87_R_star=R,
        fd_A=AA, fd_z0w=ZW, fd_CD10nf=v2._lupkes_form_drag(AA, ZW),
    )  # fmt: skip


def flux_args(v2, meteo, inp, idx, zu, zt, zq, ice_scheme="blended_a87", z0_ice=5e-4):
    """The arguments v2 passes to _ice_fluxes, computed as v2's main function does."""
    u, t, rh, P, lat = (inp[k][idx] for k in ("u", "t", "rh", "P", "lat"))
    n = u.size
    Q = meteo.qair(t, P, rh)[0] / 1000.0
    rhoa = P * 100.0 / (287.1 * (t + 273.16) * (1 + 0.61 * Q))
    return dict(
        u=u, t=t, Q_air=Q, z0_ice=z0_ice, P=P, grav=meteo.grv(lat),
        zu=zu * np.ones(n), zt=zt * np.ones(n), zq=zq * np.ones(n), rhoa=rhoa,
        cpa=1004.67, von=0.4, fdg=1.0, ta=t + 273.16, ice_scheme=ice_scheme,
        scalar_ratio=v2.SCALAR_RATIO, nits=10,
    )  # fmt: skip


ARRAY_ARGS = ("u", "t", "Q_air", "P", "grav", "zu", "zt", "zq", "rhoa", "ta")
FLUX_KEYS = ("tau", "hsb", "hlb", "Cdn", "Chn", "Cen", "R_star", "z0T", "z0q", "h_s", "h_l",
             "Qs", "L")  # fmt: skip


def most(v2, meteo, inp):
    out = {}
    t_sfc = np.array([-35.0, -20.0, -10.0, -1.8, 0.0, 2.0])
    base = np.arange(60)
    idx = np.repeat(base, t_sfc.size)
    ts_grid = np.tile(t_sfc, base.size)
    ts_grid[::37] = np.nan
    cases = {
        "a87": dict(zu=36, zt=22, zq=22),
        "a87_z10": dict(zu=10, zt=10, zq=10),
        "simple": dict(zu=36, zt=22, zq=22, ice_scheme="simple", z0_ice=1e-3),
    }
    out["t_ice"] = ts_grid
    for name, kw in cases.items():
        args = flux_args(v2, meteo, inp, idx, **kw)
        fl = v2._ice_fluxes(t_ice=ts_grid, **args)
        for k in ARRAY_ARGS:
            out[f"{name}__in_{k}"] = args[k]
        out[f"{name}__z0_ice"] = np.array(args["z0_ice"])
        for k in FLUX_KEYS:
            out[f"{name}__{k}"] = fl[k]
    return out


def seb(v2, meteo, ist, inp):
    idx = np.arange(N)
    args = flux_args(v2, meteo, inp, idx, 36, 22, 22)
    G = ist.effective_conductance(inp["h_ice"], inp["h_snow"])
    alb = ist.surface_albedo(inp["h_snow"])
    Ts, res = v2._solve_seb_steady(args, inp["Rl"], inp["Rs"], alb, G, -1.8, inp["t"])
    lw, sw, fc = v2._seb_terms(inp["t_ice"], inp["Rl"], inp["Rs"], alb, G, -1.8)
    out = {f"in_{k}": args[k] for k in ARRAY_ARGS}
    out.update(G=G, albedo=alb, Ts=Ts, res=res, lw_net=lw, sw_net=sw, F_cond=fc)
    return out


def _forcing(n, warm_from=None):
    i = np.arange(n)
    theta = -15.0 + 8.0 * np.sin(2 * np.pi * i / 24.0)
    q = 0.8 * 0.622 * 6.1115 * np.exp(22.452 * np.minimum(theta, 0) / (272.55 + theta)) / 1000
    f = dict(
        theta_a=theta, q_a=q, p=np.full(n, 1005.0),
        lw_down=220.0 + 60.0 * np.cos(2 * np.pi * i / 17.0),
        sw_down=np.maximum(0.0, 150.0 * np.sin(2 * np.pi * i / 24.0)),
        h_s=np.full(n, 15.0) + 5.0 * np.sin(i / 3.0),
        h_l=np.full(n, 4.0e4) + 1.0e4 * np.cos(i / 4.0),
    )  # fmt: skip
    if warm_from is not None:  # warm, moist, sunny spell that melts the surface
        f["theta_a"][warm_from:] = 6.0
        f["q_a"][warm_from:] = 0.006
        f["lw_down"][warm_from:] = 330.0
        f["sw_down"][warm_from:] = 500.0
    return f


def column(ist):
    out = {}
    configs = {
        "default": dict(),
        "nosnow_fresh": dict(h_ice=0.6, h_snow=0.0, n_ice=5, salinity_dependent=False),
        "thin_k": dict(h_ice=0.3, h_snow=0.05, n_ice=6, n_snow=2, k_snow=0.25, s_ice=8.0),
    }
    n = 40
    for name, kw in configs.items():
        col = ist.IceColumn(**kw)
        out[f"{name}__steady_profile"] = col.steady_profile(-15.0)
        col.initialise(-12.0)
        f = _forcing(n, warm_from=30)
        for k, v in f.items():
            out[f"steps__in_{k}"] = v
        Ts, Fc, melt, prof = [], [], [], []
        for i in range(n):
            r = col.step(600.0, *(f[k][i] for k in f))
            Ts.append(r[0])
            Fc.append(r[1])
            melt.append(r[2])
            prof.append(col.T.copy())
        out[f"{name}__step_Ts"] = np.array(Ts)
        out[f"{name}__step_F_cond"] = np.array(Fc)
        out[f"{name}__step_melt"] = np.array(melt)
        out[f"{name}__step_T"] = np.array(prof)
        # one step with an explicit albedo and a long time step
        col.initialise(-20.0)
        r = col.step(3 * 3600.0, *(f[k][5] for k in f), albedo=0.7)
        out[f"{name}__step_albedo"] = np.array(r)
    # run(): gaps, SIC below reset and NaN, NaN forcing
    m = 80
    f = _forcing(m, warm_from=70)
    time_s = np.cumsum(np.full(m, 600.0))
    time_s[30:] += 4 * 3600.0
    sic = np.full(m, 0.9)
    sic[45:48] = 0.05
    sic[55] = np.nan
    ts_init = -10.0 + 3.0 * np.sin(np.arange(m) / 5.0)
    f["lw_down"][60] = np.nan
    run_in = dict(time_s=time_s, Ts_init=ts_init, sic=sic, **f)
    for k, v in run_in.items():
        out[f"run__in_{k}"] = v
    col = ist.IceColumn()
    res = col.run(time_s, *(f[k] for k in f), ts_init, sic=sic)
    for k, v in res.items():
        out[f"run__{k}"] = v
    res = ist.IceColumn(h_snow=0.0).run(time_s, *(f[k] for k in f), ts_init, albedo=0.6)
    for k, v in res.items():
        out[f"run_nosic__{k}"] = v
    return out


# Whole-function runs of v2 coare35vn_ice. Values "@name" refer to arrays in inputs.npz.
Z36 = dict(zu=36, zt=22, zq=22)
Z10 = dict(zu=10, zt=10, zq=10)
COMMON = dict(P="@P", lat="@lat", Rs="@Rs", Rl="@Rl", verbose=False)
ICE = dict(ice_conc="@sic")
RUNS = {
    "air_freezing": dict(**Z36, **ICE, t_ice_method="air_freezing"),
    "air_capped": dict(**Z36, **ICE, t_ice_method="air_capped"),
    "empirical": dict(**Z36, **ICE, t_ice_method="empirical", t_ice_offset=3.0),
    "prescribed": dict(**Z36, **ICE, t_ice_method="prescribed", t_ice="@t_ice"),
    "prescribed_fb_air": dict(
        **Z36, **ICE, t_ice_method="prescribed", t_ice="@t_ice", t_ice_fallback="air_freezing"
    ),
    "seb_steady": dict(**Z36, **ICE, t_ice_method="seb_steady"),
    "seb_steady_arrays": dict(
        **Z36, **ICE, t_ice_method="seb_steady", h_ice="@h_ice", h_snow="@h_snow", k_snow=0.25
    ),
    "diffusion": dict(**Z36, **ICE, t_ice_method="diffusion", time_s="@time_s"),
    "simple_seb": dict(**Z36, **ICE, ice_scheme="simple", t_ice_method="seb_steady"),
    "rh_ice": dict(**Z36, **ICE, t_ice_method="seb_steady", rh_reference="ice"),
    "z10_air_freezing": dict(**Z10, **ICE, t_ice_method="air_freezing"),
    "z10_seb_steady": dict(**Z10, **ICE, t_ice_method="seb_steady"),
    "z10_rain_jcool0": dict(**Z10, **ICE, t_ice_method="air_freezing", rain="@rain", jcool=0),
}
OPENWATER = {"z36": dict(**Z36), "z10": dict(**Z10)}


def resolve(kw, inp):
    return {
        k: (inp[v[1:]] if isinstance(v, str) and v.startswith("@") else v) for k, v in kw.items()
    }


def main(v2_dir, out_dir):
    sys.path.insert(0, str(Path(v2_dir).resolve()))
    import coare35vn_ice as v2
    import ice_surface_temperature as ist
    import meteo

    assert v2.__version__ == "2.0.0", v2.__version__
    warnings.simplefilter("ignore", RuntimeWarning)  # v2 meteo.psit_26 warns harmlessly
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    inp = make_inputs(rng)
    np.savez_compressed(out / "inputs.npz", **inp)
    np.savez_compressed(out / "thermo_stability.npz", **thermo_stability(meteo, ist))
    np.savez_compressed(out / "roughness_form_drag.npz", **roughness_form_drag(v2))
    np.savez_compressed(out / "most.npz", **most(v2, meteo, inp))
    np.savez_compressed(out / "seb.npz", **seb(v2, meteo, ist, inp))
    np.savez_compressed(out / "column.npz", **column(ist))

    def call(kw):
        k = {**COMMON, **kw}
        return v2.coare35vn_ice(inp["u"], inp["t"], inp["rh"], inp["ts"], **resolve(k, inp))

    ow = {}
    for name, kw in OPENWATER.items():
        A = call(kw)
        assert A.shape == (N, 38)
        ow[name] = A
        np.savez_compressed(
            out / f"legacy_openwater_{name}.npz", A=A, kwargs=json.dumps({**COMMON, **kw})
        )
    for name, kw in RUNS.items():
        A = call(kw)
        assert A.shape == (N, 56)
        hset = "z10" if kw["zu"] == 10 else "z36"
        own_ow = "rain" in kw or "jcool" in kw
        if not own_ow:
            assert np.array_equal(A[:, :38], ow[hset], equal_nan=True), name
        data = dict(A=A) if own_ow else dict(ice=A[:, 38:], openwater=hset)
        np.savez_compressed(out / f"legacy_{name}.npz", kwargs=json.dumps({**COMMON, **kw}), **data)
    for f in sorted(out.glob("*.npz")):
        print(f"{f.name:40s} {f.stat().st_size / 1024:7.1f} kB")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    default_out = Path(__file__).resolve().parent / "v2_reference"
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else default_out)
