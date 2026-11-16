"""Open-water fluxes: a thin adapter over the public ``pycoare`` interface.

Everything here comes from the public attribute groups of a ``pycoare.coare_35`` or
``pycoare.coare_36`` instance (``fluxes``, ``velocities``, ``temperatures``,
``humidities``, ``stability_parameters``, ``transfer_coefficients``,
``stability_functions``). A few input-derived quantities that pycoare does not expose
(air density, air and surface humidity, latent heat) are recomputed with the same
formulas in :mod:`coare_seaice.thermo`.

Units: humidities in kg/kg (pycoare reports some in g/kg; they are converted here).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from pycoare import coare_35, coare_36

from . import thermo

_VERSIONS = {"3.5": coare_35, "3.6": coare_36}


@dataclass(frozen=True)
class OpenWater:
    """Open-water results. Arrays have the shape of ``u``; heat fluxes positive upward.

    Reference-height values (``*_rf``) are at ``zrf``.
    """

    version: str
    zrf: float
    # fluxes (N m-2, W m-2, mm h-1)
    tau: np.ndarray
    hsb: np.ndarray
    hlb: np.ndarray
    hbb: np.ndarray
    hsbb: np.ndarray
    hlwebb: np.ndarray
    evap: np.ndarray
    rnl: np.ndarray
    rf: np.ndarray  # rain heat flux; NaN when no rain is given (pycoare behaviour)
    # scaling parameters and stability
    usr: np.ndarray  # m/s
    tsr: np.ndarray  # K
    qsr: np.ndarray  # kg/kg
    obukhov_length: np.ndarray  # m
    zet: np.ndarray  # zu / L
    # roughness (m)
    zo: np.ndarray
    zot: np.ndarray
    zoq: np.ndarray
    # cool skin
    dter: np.ndarray  # K
    dqer: np.ndarray  # kg/kg
    tkt: np.ndarray  # m
    # transfer coefficients (not x 1e3)
    cd: np.ndarray
    ch: np.ndarray
    ce: np.ndarray
    cdn_rf: np.ndarray
    chn_rf: np.ndarray
    cen_rf: np.ndarray
    # winds (m/s)
    ut: np.ndarray  # wind including gustiness
    gust_factor: np.ndarray  # ut / du
    u_n: np.ndarray  # neutral wind at zu
    u_rf: np.ndarray
    u_n_rf: np.ndarray
    # temperature and humidity at zrf
    t_rf: np.ndarray  # deg C
    q_rf: np.ndarray  # kg/kg
    rh_rf: np.ndarray  # %
    # recomputed from the inputs with the COARE formulas (thermo)
    rhoa: np.ndarray  # kg m-3
    q_air: np.ndarray  # kg/kg, air specific humidity at zq (RH over liquid water)
    q_sea: np.ndarray  # kg/kg, saturation humidity at the sea surface
    le: np.ndarray  # J kg-1, latent heat of vaporisation at ts
    jcool: int
    coare: object  # the pycoare instance itself


def open_water(
    u,
    t,
    rh,
    ts,
    p=1015.0,
    rs=150.0,
    rl=370.0,
    zu=18.0,
    zt=18.0,
    zq=18.0,
    lat=45.0,
    zi=600.0,
    rain=None,
    cp=None,
    sigH=None,
    jcool=1,
    zrf=10.0,
    version="3.5",
    salinity=None,
):
    """Run pycoare for open water and collect the quantities the ice code needs.

    Inputs follow pycoare (u m/s at zu; t deg C at zt; rh % at zq; ts deg C; p mb;
    rs, rl downwelling W m-2; heights m; lat deg; zi m; rain mm/h). ``zrf`` is the single
    reference height for the ``*_rf`` outputs. ``version`` is "3.5" or "3.6";
    ``salinity`` (psu) is passed to COARE 3.6 only.
    """
    if version not in _VERSIONS:
        raise ValueError(f"version must be one of {sorted(_VERSIONS)}")

    # pycoare 0.4.3 changes some inputs in place (util.qair divides a float64 `rh` array
    # by 100), so pycoare is given its own copies and the arrays used here stay intact.
    def _copy(x):
        return None if x is None else np.array(x, dtype=float)

    kwargs = dict(
        t=t, rh=rh, zu=zu, zt=zt, zq=zq, zrf=zrf, ts=ts, p=p, lat=lat, zi=zi,
        rs=rs, rl=rl, rain=rain, cp=cp, sigH=sigH, jcool=jcool,
    )  # fmt: skip
    kwargs = {k: (v if k == "jcool" else _copy(v)) for k, v in kwargs.items()}
    if version == "3.6" and salinity is not None:
        kwargs["ss"] = salinity
    c = _VERSIONS[version](_copy(u), **kwargs)

    jc = 0 if jcool == 0 else 1
    q_air = thermo.specific_humidity(t, p, rh)
    if version == "3.6":
        q_sea = thermo.qsat_sea(ts, p, 35.0 if salinity is None else salinity)
    else:
        q_sea = thermo.qsat_sea(ts, p)

    f, v, sp = c.fluxes, c.velocities, c.stability_parameters
    tc, tm, hm = c.transfer_coefficients, c.temperatures, c.humidities
    return OpenWater(
        version=version,
        zrf=float(np.asarray(zrf, dtype=float).ravel()[0]),
        tau=f.tau,
        hsb=f.hsb,
        hlb=f.hlb,
        hbb=f.hbb,
        hsbb=f.hsbb,
        hlwebb=f.hlwebb,
        evap=f.evap,
        rnl=f.rnl,
        rf=f.rf,
        usr=v.usr,
        tsr=sp.tsr,
        qsr=sp.qsr,
        obukhov_length=sp.obukL,
        zet=sp.zet,
        zo=sp.zo,
        zot=sp.zot,
        zoq=sp.zoq,
        dter=tm.dter,
        dqer=hm.dqer / 1000.0,  # pycoare reports g/kg
        tkt=sp.tkt,
        cd=tc.cd,
        ch=tc.ch,
        ce=tc.ce,
        cdn_rf=tc.cdn_rf,
        chn_rf=tc.chn_rf,
        cen_rf=tc.cen_rf,
        ut=v.ut,
        gust_factor=v.gf,
        u_n=v.u_n,
        u_rf=v.u_rf,
        u_n_rf=v.u_n_rf,
        t_rf=tm.t_rf,
        q_rf=hm.q_rf / 1000.0,  # pycoare reports g/kg
        rh_rf=hm.rh_rf,
        rhoa=thermo.air_density(t, p, q_air),
        q_air=q_air,
        q_sea=q_sea,
        le=thermo.latent_heat_vaporisation(ts),
        jcool=jc,
        coare=c,
    )
