"""Shared test helpers: loading the stored v2 reference fixtures."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

REF_DIR = Path(__file__).parent / "fixtures" / "v2_reference"

# Tolerance for "ported function equals v2 on identical inputs".
RTOL_PORT = 1e-12


def load_ref(name):
    """Load ``tests/fixtures/v2_reference/<name>.npz`` as a dict of arrays."""
    with np.load(REF_DIR / f"{name}.npz") as f:
        return {k: f[k] for k in f.files}


def resolve_kwargs(kwargs_json, inputs):
    """Turn a stored kwargs JSON string into call arguments ("@name" -> inputs[name])."""
    kw = json.loads(str(kwargs_json))
    return {
        k: (inputs[v[1:]] if isinstance(v, str) and v.startswith("@") else v) for k, v in kw.items()
    }


def assert_port(actual, expected, rtol=RTOL_PORT):
    """Equal to v2 to relative tolerance ``rtol``.

    An absolute floor of ``rtol`` times the largest magnitude in ``expected`` is added, so
    that values that cancel to near zero (e.g. psi at zeta ~ 0) do not fail on last-bit
    differences between numpy builds. With the numpy used to make the fixtures the
    results are bit-identical.
    """
    expected = np.asarray(expected, dtype=float)
    scale = np.nanmax(np.abs(expected)) if np.isfinite(expected).any() else 0.0
    np.testing.assert_allclose(actual, expected, rtol=rtol, atol=rtol * scale, equal_nan=True)


@pytest.fixture(scope="session")
def inputs():
    return load_ref("inputs")
