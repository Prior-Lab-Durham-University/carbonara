#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")" && pwd)
EXT="$ROOT/external"
REQ="$ROOT/requirements.txt"

mkdir -p "$EXT"

echo ">>> Project root: $ROOT"
echo ">>> External dir : $EXT"

# ----------------------------
# Python dependencies
# ----------------------------
if [[ -f "$REQ" ]]; then
    echo ">>> Installing Python requirements..."
    python3 -m pip install -r "$REQ"
else
    echo "WARNING: requirements.txt not found at $REQ"
fi


# ----------------------------
# pyFoXS setup (required)
# ----------------------------
PYFOXS_TARGET="$EXT/pyFoXS"
PYFOXS_GIT_URL="https://github.com/biocatiit/pyFoXS.git"

echo ">>> Installing pyFoXS dependencies..."
python3 -m pip install numpy scipy numba

if [[ -d "$PYFOXS_TARGET/.git" ]]; then
    echo ">>> pyFoXS already present, pulling latest changes..."
    git -C "$PYFOXS_TARGET" pull --ff-only
else
    echo ">>> Cloning pyFoXS from $PYFOXS_GIT_URL"
    rm -rf "$PYFOXS_TARGET"
    git clone "$PYFOXS_GIT_URL" "$PYFOXS_TARGET"
fi

if [[ ! -f "$PYFOXS_TARGET/pyFoXS/foxs.py" ]]; then
    echo "ERROR: pyFoXS installation failed."
    echo "Expected file not found:"
    echo "  $PYFOXS_TARGET/pyFoXS/foxs.py"
    exit 1
fi

echo ">>> pyFoXS ready at: $PYFOXS_TARGET/pyFoXS/foxs.py"


echo
echo ">>> Setup complete."
echo ">>> Next steps:"
echo "    1) Ensure MODELLER is installed and licensed"
echo "    2) Ensure predictStructureQvary is built"
echo "    3) Launch from notebook or shell"
