#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash install_cg2all_local.sh /extra/tmp/carbonara-pseudoWaxsis/cg2all_local

BASE="${1:-/extra/tmp/carbonara-pseudoWaxsis/cg2all_local}"

MICROMAMBA="$BASE/bin/micromamba"
MAMBA_ROOT_PREFIX="$BASE/micromamba"
CG2ALL_ENV="$MAMBA_ROOT_PREFIX/envs/cg2all"
CG2ALL_BIN="$CG2ALL_ENV/bin/convert_cg2all"

TMPDIR_LOCAL="$BASE/tmp"
PIP_CACHE_DIR_LOCAL="$BASE/pip-cache"

mkdir -p "$BASE/bin" "$MAMBA_ROOT_PREFIX" "$TMPDIR_LOCAL" "$PIP_CACHE_DIR_LOCAL"

export MAMBA_ROOT_PREFIX
export TMPDIR="$TMPDIR_LOCAL"
export PIP_CACHE_DIR="$PIP_CACHE_DIR_LOCAL"

echo ">>> BASE              : $BASE"
echo ">>> MICROMAMBA        : $MICROMAMBA"
echo ">>> MAMBA_ROOT_PREFIX : $MAMBA_ROOT_PREFIX"
echo ">>> CG2ALL_ENV        : $CG2ALL_ENV"
echo ">>> TMPDIR            : $TMPDIR"
echo ">>> PIP_CACHE_DIR     : $PIP_CACHE_DIR"

# -------------------------------------------------------------------
# 1. Install micromamba locally
# -------------------------------------------------------------------
if [[ ! -x "$MICROMAMBA" ]]; then
    echo
    echo ">>> Installing micromamba locally..."

    MM_TAR="$TMPDIR_LOCAL/micromamba.tar.bz2"

    if command -v wget >/dev/null 2>&1; then
        wget -qO "$MM_TAR" https://micro.mamba.pm/api/micromamba/linux-64/latest
    elif command -v curl >/dev/null 2>&1; then
        curl -L -o "$MM_TAR" https://micro.mamba.pm/api/micromamba/linux-64/latest
    else
        echo "ERROR: need wget or curl to download micromamba."
        exit 1
    fi

    tar -xjf "$MM_TAR" -C "$BASE" bin/micromamba
else
    echo
    echo ">>> micromamba already present: $MICROMAMBA"
fi

"$MICROMAMBA" --help >/dev/null

# -------------------------------------------------------------------
# 2. Create CG2ALL environment
# -------------------------------------------------------------------
if [[ ! -x "$CG2ALL_BIN" ]]; then
    echo
    echo ">>> Creating cg2all micromamba environment..."

    "$MICROMAMBA" create -y -p "$CG2ALL_ENV" python=3.10 pip

    "$MICROMAMBA" run -p "$CG2ALL_ENV" python -m pip install --upgrade pip setuptools wheel

    echo
    echo ">>> Installing pinned compatibility stack..."

    "$MICROMAMBA" run -p "$CG2ALL_ENV" python -m pip install --no-cache-dir "numpy==1.26.4"
    "$MICROMAMBA" run -p "$CG2ALL_ENV" python -m pip install --no-cache-dir "torch==1.13.1"
    "$MICROMAMBA" run -p "$CG2ALL_ENV" python -m pip install --no-cache-dir "dgl==0.9.1" -f https://data.dgl.ai/wheels/repo.html
    "$MICROMAMBA" run -p "$CG2ALL_ENV" python -m pip install --no-cache-dir "e3nn==0.4.4"
    "$MICROMAMBA" run -p "$CG2ALL_ENV" python -m pip install --no-cache-dir "ml-collections==0.1.1"

    echo
    echo ">>> Installing CG2ALL-related dependencies..."

    "$MICROMAMBA" run -p "$CG2ALL_ENV" python -m pip install --no-cache-dir git+https://github.com/huhlim/mdtraj.git
    "$MICROMAMBA" run -p "$CG2ALL_ENV" python -m pip install --no-cache-dir git+https://github.com/huhlim/SE3Transformer.git

    echo
    echo ">>> Installing cg2all..."

    "$MICROMAMBA" run -p "$CG2ALL_ENV" python -m pip install --no-cache-dir --no-deps \
        git+https://github.com/huhlim/cg2all@a00b8816736c08852944f147e39164d0f5e1834e
else
    echo
    echo ">>> cg2all already installed: $CG2ALL_BIN"
fi

# -------------------------------------------------------------------
# 3. Verify
# -------------------------------------------------------------------
echo
echo ">>> Testing convert_cg2all..."

"$MICROMAMBA" run -p "$CG2ALL_ENV" convert_cg2all --help >/dev/null

echo
echo ">>> CG2ALL environment is ready."
echo
echo ">>> Use this executable form:"
echo "    $MICROMAMBA run -p $CG2ALL_ENV convert_cg2all"
echo
echo ">>> Python in CG2ALL env:"
echo "    $MICROMAMBA run -p $CG2ALL_ENV python"
