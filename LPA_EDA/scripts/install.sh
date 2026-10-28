#!/usr/bin/env bash
set -euo pipefail
cd /workspace/APL-ML---Project/LPA_EDA
export PIP_CACHE_DIR=/workspace/.cache/lpa-eda/pip
export MPLCONFIGDIR=/workspace/.cache/lpa-eda/matplotlib
export XDG_CACHE_HOME=/workspace/.cache/lpa-eda
mkdir -p "$PIP_CACHE_DIR" "$MPLCONFIGDIR" /workspace/.venvs
if [ ! -x /workspace/.venvs/lpa-eda/bin/python ]; then
  python3 -m venv --system-site-packages /workspace/.venvs/lpa-eda
fi
/workspace/.venvs/lpa-eda/bin/python -m pip install --requirement requirements.txt
/workspace/.venvs/lpa-eda/bin/python -m pip check
