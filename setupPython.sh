#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------------------
# Carbonara Python setup
#
# This script is intended to be the supported install route for the
# Python-side Carbonara tooling. It deliberately does more than
# `pip install -r requirements.txt` because some dependencies need
# ordered/special handling:
#
#   - biobox needs numpy available at build time, so it is installed
#     separately with --no-build-isolation.
#   - pyFoXS is cloned from GitHub and needs torch at runtime.
#
# Usage:
#   source .venv/bin/activate
#   bash setupPython.sh
#
# Optional for quota-limited machines:
#   export PIP_INSTALL_OPTS="--no-cache-dir --no-compile"
#   bash setupPython.sh
# ---------------------------------------------------------------------

ROOT=$(cd "$(dirname "$0")" && pwd)
EXT="$ROOT/external"
REQ="$ROOT/requirements.txt"

# Use the active Python environment by default.
# Users can override explicitly:
#   PYTHON=/path/to/python bash setupPython.sh
PYTHON="${PYTHON:-python}"

# Extra pip options, e.g. "--no-cache-dir --no-compile"
PIP_INSTALL_OPTS="${PIP_INSTALL_OPTS:-}"

mkdir -p "$EXT"

echo ">>> Project root : $ROOT"
echo ">>> External dir : $EXT"
echo ">>> Python       : $($PYTHON -c 'import sys; print(sys.executable)')"
echo ">>> Python prefix: $($PYTHON -c 'import sys; print(sys.prefix)')"

# ----------------------------
# Bootstrap build dependencies
# ----------------------------
echo
echo ">>> Installing bootstrap build dependencies..."
"$PYTHON" -m pip install $PIP_INSTALL_OPTS --upgrade pip setuptools wheel
"$PYTHON" -m pip install $PIP_INSTALL_OPTS numpy Cython

# ----------------------------
# Python dependencies
# ----------------------------
if [[ -f "$REQ" ]]; then
    echo
    echo ">>> Installing Python requirements..."
    "$PYTHON" -m pip install $PIP_INSTALL_OPTS -r "$REQ"
else
    echo "WARNING: requirements.txt not found at $REQ"
fi

# ----------------------------
# biobox special handling
# ----------------------------
echo
echo ">>> Installing biobox..."
echo ">>> Note: biobox is installed separately because it needs numpy at build time."
"$PYTHON" -m pip install $PIP_INSTALL_OPTS --no-build-isolation biobox

# ----------------------------
# pyFoXS setup
# ----------------------------
PYFOXS_TARGET="$EXT/pyFoXS"
PYFOXS_GIT_URL="https://github.com/biocatiit/pyFoXS.git"

echo
echo ">>> Ensuring pyFoXS runtime dependencies..."
"$PYTHON" -m pip install $PIP_INSTALL_OPTS numpy scipy numba torch

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

# Create a pyfoxs wrapper inside the active environment's script/bin dir.
BIN_DIR="$($PYTHON -c 'import sysconfig; print(sysconfig.get_path("scripts"))')"
PYFOXS_WRAPPER="$BIN_DIR/pyfoxs"

echo ">>> Creating pyfoxs wrapper at: $PYFOXS_WRAPPER"

cat > "$PYFOXS_WRAPPER" <<EOF
#!/usr/bin/env bash
exec "$PYTHON" "$PYFOXS_TARGET/pyFoXS/foxs.py" "\$@"
EOF

chmod +x "$PYFOXS_WRAPPER"

# Convenience wrapper to run Carbonara scripts with the configured Python.
CARBONARA_PYTHON_WRAPPER="$BIN_DIR/carbonara-python"

echo ">>> Creating carbonara-python wrapper at: $CARBONARA_PYTHON_WRAPPER"

cat > "$CARBONARA_PYTHON_WRAPPER" <<EOF
#!/usr/bin/env bash
exec "$PYTHON" "\$@"
EOF

chmod +x "$CARBONARA_PYTHON_WRAPPER"

# ----------------------------
# Import smoke test
# ----------------------------
echo
echo ">>> Running import smoke test..."

"$PYTHON" - <<'PY'
import importlib
import sys

print("Python:", sys.executable)
print("Prefix:", sys.prefix)

checks = [
    ("numpy", "numpy"),
    ("scipy", "scipy"),
    ("matplotlib", "matplotlib"),
    ("pandas", "pandas"),
    ("Bio", "Bio"),
    ("IPython", "IPython"),
    ("ipywidgets", "ipywidgets"),
    ("tqdm", "tqdm"),
    ("plotly", "plotly"),
    ("networkx", "networkx"),
    ("sklearn", "sklearn"),
    ("numba", "numba"),
    ("Cython", "Cython"),
    ("mdtraj", "mdtraj"),
    ("openmm", "openmm"),
    ("pdbfixer", "pdbfixer"),
    ("py3Dmol", "py3Dmol"),
    ("torch", "torch"),
    ("biobox", "biobox"),
]

bad = []

for label, module_name in checks:
    try:
        mod = importlib.import_module(module_name)
        print(f"OK  {label:20s} {getattr(mod, '__version__', 'version unknown')}")
    except Exception as exc:
        print(f"BAD {label:20s} {type(exc).__name__}: {exc}")
        bad.append(label)

# Exact imports used by Carbonara/all-atom scripts.
exact_checks = [
    ("pdbfixer.PDBFixer", "from pdbfixer import PDBFixer"),
    ("openmm.app.PDBFile", "from openmm.app import PDBFile"),
    ("CarbonaraDataTools", "import CarbonaraDataTools"),
]

for label, statement in exact_checks:
    try:
        exec(statement, {})
        print(f"OK  {label:20s}")
    except Exception as exc:
        print(f"BAD {label:20s} {type(exc).__name__}: {exc}")
        bad.append(label)

if bad:
    raise SystemExit("Import smoke test failed for: " + ", ".join(bad))

print(">>> Import smoke test passed.")
PY

# ----------------------------
# pyFoXS smoke test
# ----------------------------
echo
echo ">>> Running pyFoXS smoke test..."

"$PYFOXS_WRAPPER" --help >/dev/null || {
    echo "ERROR: pyfoxs wrapper failed."
    echo "Tried: $PYFOXS_WRAPPER --help"
    exit 1
}

echo ">>> pyFoXS smoke test passed."

echo
echo ">>> Setup complete."
echo ">>> Python             : $($PYTHON -c 'import sys; print(sys.executable)')"
echo ">>> pyFoXS script      : $PYFOXS_TARGET/pyFoXS/foxs.py"
echo ">>> pyfoxs wrapper     : $PYFOXS_WRAPPER"
echo ">>> carbonara-python   : $CARBONARA_PYTHON_WRAPPER"
echo
echo ">>> Next steps:"
echo "    1) Ensure predictStructureQvary is built:"
echo "         mkdir -p build && cd build && cmake .. && make -j4"
echo "    2) In notebooks, run scripts using the kernel Python:"
echo "         import sys"
echo "         !{sys.executable} setup_carbonara_allAtom.py ..."
echo "    3) If using Modeller backmapping, ensure MODELLER is installed and licensed."
