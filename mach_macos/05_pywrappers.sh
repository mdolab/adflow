#!/usr/bin/env bash
# ======================================================================
#  05_pywrappers.sh  —  Python packages ADflow needs at runtime
#  Requires Steps 1-4 (env + OpenMPI + real PETSc + CGNS).
#
#  Installs into the 'adflow' conda env:
#    * mpi4py 3.1.6        (built against YOUR OpenMPI)
#    * mdolab-baseclasses  (pure-python; ADflow imports it)
#
#  (No petsc4py — the non-eigenvalue ADflow does not import it; it links
#   the PETSc C library directly in the Fortran build.)
#
#  Self-contained:  chmod +x 05_pywrappers.sh && ./05_pywrappers.sh
#  ~2-4 min. Safe to re-run.
# ======================================================================
set -eo pipefail
log()  { printf '\n\033[1;34m==>\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[warn]\033[0m %s\n' "$*"; }
die()  { printf '\033[1;31m[error]\033[0m %s\n' "$*" >&2; exit 1; }

export MINIFORGE_HOME="${MINIFORGE_HOME:-$HOME/miniforge3}"
[[ -f "$MINIFORGE_HOME/etc/profile.d/conda.sh" ]] || die "run 01_setup.sh first."
source "$MINIFORGE_HOME/etc/profile.d/conda.sh"
conda activate adflow || die "conda env 'adflow' missing — run 01_setup.sh first."

export PACKAGES="${PACKAGES:-$HOME/packages}"
export MPI_INSTALL_DIR="${MPI_INSTALL_DIR:-$PACKAGES/openmpi-5.0.9/opt-gfortran}"
export PATH="$MPI_INSTALL_DIR/bin:$PATH"; hash -r
[[ -x "$MPI_INSTALL_DIR/bin/mpicc" ]] || die "OpenMPI missing — run 02_openmpi.sh."

# --- mpi4py (against our OpenMPI) -------------------------------------
if python -c "import mpi4py" 2>/dev/null; then
  log "mpi4py already installed ($(python -c 'import mpi4py;print(mpi4py.__version__)'))."
else
  # mpi4py 3.1.6's bundled mpidistutils calls new_compiler(dry_run=...),
  # which setuptools>=74 removed -> "unexpected keyword argument 'dry_run'".
  # Pin setuptools<74 (the working machine has 73.0.1) and build WITHOUT
  # pip's build isolation so that pinned setuptools is actually used.
  log "Pinning setuptools<74 for the mpi4py build…"
  pip install "setuptools<74" wheel
  log "Installing mpi4py 3.1.6 (compiles against $(command -v mpicc))…"
  MPICC="$MPI_INSTALL_DIR/bin/mpicc" \
    pip install "mpi4py==3.1.6" --no-cache --no-build-isolation
fi

# --- mdolab-baseclasses ----------------------------------------------
if python -c "import baseclasses" 2>/dev/null; then
  log "baseclasses already installed."
else
  log "Installing mdolab-baseclasses 1.8.5…"
  pip install "mdolab-baseclasses==1.8.5"
fi

# --- verify -----------------------------------------------------------
log "Import check:"
python - <<'PY'
import mpi4py, baseclasses
print("  mpi4py     ", mpi4py.__version__)
print("  baseclasses", baseclasses.__version__)
PY
log "Parallel MPI check (2 ranks):"
"$MPI_INSTALL_DIR/bin/mpirun" -n 2 python -c \
  "from mpi4py import MPI; print('  rank', MPI.COMM_WORLD.rank, 'of', MPI.COMM_WORLD.size)"

log "Step 5 complete — mpi4py + baseclasses installed."
echo; echo "NEXT:  ./06_adflow.sh"
