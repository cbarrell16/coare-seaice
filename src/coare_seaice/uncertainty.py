"""Uncertainty by perturbed-input ensembles.

:func:`ensemble` reruns :func:`coare_seaice.ice_fluxes` many times with inputs or
settings perturbed, and returns percentiles of chosen outputs. It measures how input
uncertainty spreads into the fluxes. It says nothing about errors in the method itself.

Perturbations are made with :func:`normal`, :func:`uniform` and :func:`log_uniform`:

    ens = ensemble(
        dict(u=u, t=t, rh=rh, p=p, sst=sst, sic=sic),
        dict(t=normal(0.2), u=normal(0.05, relative=True), sic=normal(0.1),
             z0_ice=log_uniform(1e-4, 1e-2)),
        n=200, seed=1, zu=36, zt=22, zq=22, lat=lat, lw_down=lw, sw_down=sw,
    )
    ens.percentiles["blended.hsb"]      # shape (3, N): 5th, 50th, 95th percentile
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np

from .api import ice_fluxes

__all__ = ["EnsembleResult", "ensemble", "log_uniform", "normal", "uniform"]

DEFAULT_OUTPUTS = ("blended.tau", "blended.hsb", "blended.hlb")


def normal(sd, relative=False, per_sample=False):
    """Gaussian error with standard deviation ``sd`` added to the value.

    relative=True: ``sd`` is a fraction of the value (0.05 = 5 %). per_sample=False (the
    default) draws one error per ensemble member and applies it to every sample, as for
    a calibration offset; per_sample=True draws an independent error for every sample,
    as for random noise.
    """

    def draw(rng, value):
        value = np.asarray(value, dtype=float)
        e = rng.normal(0.0, 1.0, value.shape if per_sample else ()) * sd
        return value * (1.0 + e) if relative else value + e

    return draw


def uniform(low, high):
    """Replace the value by one draw per member from a uniform distribution."""
    return lambda rng, value: rng.uniform(low, high)


def log_uniform(low, high):
    """Replace the value by one draw per member, uniform in the logarithm (for roughness)."""
    return lambda rng, value: float(np.exp(rng.uniform(np.log(low), np.log(high))))


@dataclass
class EnsembleResult:
    """``central`` is the unperturbed FluxResult; ``percentiles[name]`` has shape
    (len(levels), N); ``members[name]`` has shape (n, N) if members were kept;
    ``n_valid[name]`` (shape N) counts the members that gave a finite value."""

    central: object
    levels: tuple
    percentiles: dict
    members: dict | None
    n: int
    seed: int
    n_valid: dict | None = None


def _get(result, path):
    obj = result
    for part in path.split("."):
        obj = obj[part] if isinstance(obj, dict) else getattr(obj, part)
    return np.asarray(obj, dtype=float)


def ensemble(inputs, perturbations, n=100, seed=0, levels=(5, 50, 95),
             outputs=DEFAULT_OUTPUTS, keep_members=False, **settings):  # fmt: skip
    """Percentiles of flux outputs under perturbed inputs.

    inputs : dict with u, t, rh, p, sst, sic (the positional inputs of ``ice_fluxes``).
    perturbations : dict of name -> perturbation (see :func:`normal`, :func:`uniform`,
        :func:`log_uniform`), where name is an input or a keyword setting such as
        ``z0_ice``, ``lw_down`` or ``zu``. A perturbed ``sic`` is limited to 0-1 and a
        perturbed ``rh`` and ``u`` to be non-negative, and perturbed heights to at least
    0.1 m; nothing else is limited.
    n, seed : ensemble size and random seed; the same seed gives the same result.
    levels : percentiles to return. outputs : dotted names in the result, for example
        "blended.hsb" or "tiles.thick_ice.t_surface".
    Other keywords go to ``ice_fluxes`` unchanged.

    NaN members are ignored in the percentiles (a sample is NaN only if every member
    is); ``n_valid`` reports how many members each percentile is based on.
    """
    order = ("u", "t", "rh", "p", "sst", "sic")
    missing = [k for k in order if k not in inputs]
    if missing:
        raise ValueError(f"ensemble: inputs is missing {missing}")
    unknown = [k for k in perturbations if k not in inputs and k not in settings]
    if unknown:
        raise ValueError(f"ensemble: nothing to perturb for {unknown}; give its central value")
    rng = np.random.default_rng(seed)
    central = ice_fluxes(*(inputs[k] for k in order), **settings)
    stack = {name: np.empty((n, central.n)) for name in outputs}
    for i in range(n):
        inp, kw = dict(inputs), dict(settings)
        for name, draw in perturbations.items():
            target = inp if name in inp else kw
            target[name] = draw(rng, target[name])
        inp["sic"] = np.clip(inp["sic"], 0.0, 1.0)
        inp["rh"] = np.maximum(inp["rh"], 0.0)
        inp["u"] = np.maximum(inp["u"], 0.0)  # a speed, not a component
        for z in ("zu", "zt", "zq"):
            if z in perturbations:
                kw[z] = np.maximum(kw[z], 0.1)
        member = ice_fluxes(*(inp[k] for k in order), **kw)
        for name in outputs:
            stack[name][i] = _get(member, name)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)  # all-NaN samples
        pct = {name: np.nanpercentile(v, levels, axis=0) for name, v in stack.items()}
    n_valid = {name: np.isfinite(v).sum(axis=0) for name, v in stack.items()}
    return EnsembleResult(central=central, levels=tuple(levels), percentiles=pct,
                          members=stack if keep_members else None, n=n, seed=seed,
                          n_valid=n_valid)  # fmt: skip
