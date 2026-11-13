"""Smoke tests: the package imports and its one hard dependency is usable."""

import numpy as np


def test_package_imports():
    import coare_seaice

    assert isinstance(coare_seaice.__version__, str)


def test_pycoare_open_water_runs():
    from pycoare import coare_35

    c = coare_35(
        [8.0], t=[-10.0], rh=[80.0], ts=[-1.0], zu=10, zt=10, zq=10, p=[1000.0], lat=[75.0]
    )
    assert np.isfinite(np.asarray(c.fluxes.hsb)).all()
    assert np.asarray(c.fluxes.hsb)[0] > 0  # cold air over warmer water: upward heat flux
