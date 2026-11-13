#!/usr/bin/env bash
# Setup script for Claude Code on the web (cloud environment).
# Needs only PyPI and GitHub, which the default network allowlist permits.
# Data downloads and campaign runs happen on the UEA machine, not here.
set -euo pipefail
python -m pip install --upgrade pip
python -m pip install -e ".[dev,verify,docs]"
pre-commit install || true
python -c "import coare_seaice, pycoare; print('coare-seaice', coare_seaice.__version__)"
