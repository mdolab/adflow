#!/usr/bin/env bash
# ======================================================================
#  07_smoketest.sh  —  Final end-to-end check of the ADflow install
#  Requires Steps 1-6.
#
#  Self-contained:  chmod +x 07_smoketest.sh && ./07_smoketest.sh
# ======================================================================
set -eo pipefail
log()  { printf '\n\033[1;34m==>\033[0m %s\n' "$*"; }
die()  { printf '\033[1;31m[error]\033[0m %s\n' "$*" >&2; exit 1; }

export MINIFORGE_HOME="${MINIFORGE_HOME:-$HOME/miniforge3}"
source "$MINIFORGE_HOME/etc/profile.d/conda.sh"
conda activate adflow || die "conda env 'adflow' missing."
export PACKAGES="${PACKAGES:-$HOME/packages}"
export MPI_INSTALL_DIR="${MPI_INSTALL_DIR:-$PACKAGES/openmpi-5.0.9/opt-gfortran}"
export PATH="$MPI_INSTALL_DIR/bin:$PATH"; hash -r

log "1) Versions"
python - <<'PY'
import mpi4py, baseclasses, adflow
print("   python     ", __import__("sys").version.split()[0])
print("   mpi4py     ", mpi4py.__version__)
print("   baseclasses", baseclasses.__version__)
print("   adflow     ", getattr(adflow, "__version__", "(installed)"))
PY

log "2) ADFLOW class imports"
python -c "from adflow import ADFLOW; print('   OK')"

log "3) Parallel run under MPI (2 ranks)"
mpirun -n 2 python -c "from mpi4py import MPI; c=MPI.COMM_WORLD; print('   rank', c.rank, 'of', c.size)"

log "4) Native libs libadflow.so resolves"
otool -L "$(python -c 'import adflow,os;print(os.path.join(os.path.dirname(adflow.__file__),"libadflow.so"))')" \
  | grep -iE "cgns|petsc|mpi" | sed 's/^/   /'

log "ALL CHECKS PASSED — ADflow is ready on this Mac."
echo
echo "To use it in a new terminal:  bash -l   then   conda activate adflow"
