# coare-seaice

Turbulent air-sea-ice fluxes over sea ice and the marginal ice zone, as an extension of
the COARE bulk flux algorithm.

**Status: early development. Not yet released. Interfaces and results will change.**

## What it does

Given standard surface-layer measurements (wind, temperature, humidity, pressure), sea
surface temperature and sea-ice concentration, it returns area-mean momentum, sensible
heat and latent heat fluxes over a mixed surface of water and ice.

- **Open water:** COARE 3.5 or 3.6 through [`pycoare`](https://github.com/pyCOARE/coare).
- **Momentum over ice:** form drag after Lüpkes et al. (2012) with the Elvidge et al.
  (2016) parameters.
- **Heat and moisture over ice:** the Andreas (1987) scalar roughness model.
- **Ice surface temperature:** a simple air-temperature assumption, a prescribed value, or
  a surface energy balance with conduction through ice and snow.
- **Surface tiles:** open water, thick ice and thin ice, blended by area.

## Install (development)

```bash
git clone https://github.com/cbarrell16/coare-seaice
cd coare-seaice
pip install -e ".[dev,verify]"
pytest
```

Python 3.10 or later.

## Documentation

- `docs/design/` — software design, verification study design and paper outline.
- `docs/caveats.md` — scientific caveats that apply to every result.
- `TASKS.md` — the development backlog.
- `CHANGELOG.md` — changes that alter results.

## Scope and limits

Cold season only. The accuracy of the fluxes is limited by the sea-ice concentration, the
ice surface temperature, the ice roughness and, for the energy-balance method, the
downwelling longwave radiation. See `docs/caveats.md`.

## Citing

See `CITATION.cff`. Authors: Chris Barrell and Ian Renfrew (University of East Anglia).

## AI assistance

Parts of the code, tests and documentation were drafted with Claude large language models
(Anthropic), under the direction and review of the authors. See `docs/ai_assistance.md`.

## Licence

MIT. See `LICENSE`.
