#!/usr/bin/env bash
# ======================================================================
#  06_adflow.sh  —  Build + install ADflow (REAL build, no eigenvalue code)
#  Requires Steps 1-5 (env + OpenMPI + PETSc + CGNS + py wrappers).
#
#  Clones the repo, drops in the macOS config, runs `make` (real), then
#  `pip install .` into the 'adflow' conda env.
#
#  NEEDS the file  adflow_config_friend.mk  next to this script.
#  Self-contained otherwise:  chmod +x 06_adflow.sh && ./06_adflow.sh
#  ~10-20 min. Safe to re-run.
# ======================================================================
set -eo pipefail
log()  { printf '\n\033[1;34m==>\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[warn]\033[0m %s\n' "$*"; }
die()  { printf '\033[1;31m[error]\033[0m %s\n' "$*" >&2; exit 1; }

# ----------------------------------------------------------------------
#  ADflow repo — public mdolab upstream (no auth, no SSH key needed).
#  Override at runtime if ever needed:  ADFLOW_REPO=<url> ./06_adflow.sh
# ----------------------------------------------------------------------
ADFLOW_REPO="${ADFLOW_REPO:-https://github.com/mdolab/adflow.git}"

CFG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CFG_FILE="$CFG_DIR/adflow_config_friend.mk"
[[ -f "$CFG_FILE" ]] || die "missing $CFG_FILE — keep it next to this script."

export MINIFORGE_HOME="${MINIFORGE_HOME:-$HOME/miniforge3}"
[[ -f "$MINIFORGE_HOME/etc/profile.d/conda.sh" ]] || die "run 01_setup.sh first."
source "$MINIFORGE_HOME/etc/profile.d/conda.sh"
conda activate adflow || die "conda env 'adflow' missing — run 01_setup.sh first."

# --- env the build needs (self-contained defaults) --------------------
export PACKAGES="${PACKAGES:-$HOME/packages}"
export MPI_INSTALL_DIR="${MPI_INSTALL_DIR:-$PACKAGES/openmpi-5.0.9/opt-gfortran}"
export PETSC_DIR="${PETSC_DIR_REAL:-$PACKAGES/petsc-3.21.0}"
export PETSC_ARCH="real-debug"
export CGNS_HOME="${CGNS_HOME:-$PACKAGES/CGNS-4.5.0/opt-gfortran}"
export PATH="$MPI_INSTALL_DIR/bin:$CGNS_HOME/bin:$PATH"; hash -r

[[ -x "$MPI_INSTALL_DIR/bin/mpifort" ]]                || die "OpenMPI missing — run 02."
[[ -f "$PETSC_DIR/$PETSC_ARCH/lib/libpetsc.dylib" ]]   || die "real PETSc missing — run 03."
[[ -f "$CGNS_HOME/lib/libcgns.dylib" ]]                || die "CGNS missing — run 04."
python -c "import mpi4py, baseclasses" 2>/dev/null || die "py wrappers missing — run 05."

# --- clone -----------------------------------------------------------
cd "$PACKAGES"
if [[ ! -d adflow ]]; then
  log "Cloning ADflow from $ADFLOW_REPO …"
  git clone "$ADFLOW_REPO" adflow || die "git clone failed — check repo URL / SSH access."
fi
cd adflow

# --- config + build --------------------------------------------------
log "Installing macOS config.mk (real, no SLEPc/complexify)…"
cp "$CFG_FILE" config/config.mk

log "Building ADflow (real)… this compiles the Fortran core, ~10-20 min."
make

[[ -f adflow/libadflow.so ]] || die "adflow/libadflow.so not produced — build failed."
log "libadflow.so built."

# --- python install --------------------------------------------------
log "Installing the adflow python package…"
pip install .

# --- verify ----------------------------------------------------------
log "What libadflow.so links against:"
otool -L adflow/libadflow.so | grep -iE "cgns|petsc|mpi|metis|superlu" | sed 's/^/    /'
log "Import check:"
python -c "from adflow import ADFLOW; print('  ADFLOW import OK')"

log "Step 6 complete — ADflow built and importable."
echo; echo "NEXT:  ./07_smoketest.sh"
