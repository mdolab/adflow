#!/usr/bin/env bash
# ======================================================================
#  02_openmpi.sh  —  Build OpenMPI 5.0.9 from source
#  Step 2 of 10.  Requires Step 1 done (conda env 'adflow' + toolchain).
#
#  Installs to:  $HOME/packages/openmpi-5.0.9/opt-gfortran
#  Compilers  :  conda clang / clang++ / gfortran  (NOT Apple/Homebrew)
#
#  Self-contained: activates the env itself. Just:
#    chmod +x 02_openmpi.sh
#    ./02_openmpi.sh
#
#  Takes ~10-20 min. Safe to re-run (skips the build if already installed).
# ======================================================================
set -eo pipefail

log()  { printf '\n\033[1;34m==>\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[warn]\033[0m %s\n' "$*"; }
die()  { printf '\033[1;31m[error]\033[0m %s\n' "$*" >&2; exit 1; }

# --- activate the conda toolchain -------------------------------------
export MINIFORGE_HOME="${MINIFORGE_HOME:-$HOME/miniforge3}"
[[ -f "$MINIFORGE_HOME/etc/profile.d/conda.sh" ]] || \
  die "miniforge not found at $MINIFORGE_HOME — run 01_setup.sh first."
source "$MINIFORGE_HOME/etc/profile.d/conda.sh"
conda activate adflow || die "conda env 'adflow' missing — run 01_setup.sh first."
command -v gfortran >/dev/null || die "gfortran not on PATH — is the 'adflow' env active?"

# --- paths (match the MACH bashrc block; self-contained defaults) ------
export PACKAGES="${PACKAGES:-$HOME/packages}"
export MPI_INSTALL_DIR="${MPI_INSTALL_DIR:-$PACKAGES/openmpi-5.0.9/opt-gfortran}"
VER=5.0.9
TARBALL="openmpi-${VER}.tar.bz2"
URL="https://download.open-mpi.org/release/open-mpi/v5.0/${TARBALL}"
NPROC="$(sysctl -n hw.ncpu)"
mkdir -p "$PACKAGES"

# --- compilers: exactly what the working mac used ----------------------
export CC=clang CXX=clang++ FC=gfortran F77=gfortran
log "Toolchain:"
echo "    CC=$(command -v clang)"
echo "    FC=$(command -v gfortran)"
echo "    prefix=$MPI_INSTALL_DIR"

# --- build (skip if already installed) --------------------------------
if [[ -x "$MPI_INSTALL_DIR/bin/mpicc" && -x "$MPI_INSTALL_DIR/bin/mpifort" ]]; then
  log "OpenMPI already installed at $MPI_INSTALL_DIR — skipping build."
else
  cd "$PACKAGES"
  if [[ ! -f "$TARBALL" ]]; then
    log "Downloading OpenMPI $VER…"
    curl -fL -o "$TARBALL" "$URL" || wget -O "$TARBALL" "$URL" || die "download failed."
  fi
  log "Unpacking…"
  rm -rf "openmpi-${VER}"
  tar xjf "$TARBALL"
  cd "openmpi-${VER}"

  log "Configuring (this is the slow part)…"
  ./configure --prefix="$MPI_INSTALL_DIR"

  log "Building with -j$NPROC…"
  make -j"$NPROC"

  log "Installing…"
  make install
fi

# --- put it on PATH for verification ----------------------------------
export PATH="$MPI_INSTALL_DIR/bin:$PATH"
hash -r

# --- verify -----------------------------------------------------------
log "Verifying wrapper compilers point at the conda toolchain:"
"$MPI_INSTALL_DIR/bin/mpicc"   --showme:command
"$MPI_INSTALL_DIR/bin/mpifort" --showme:command
echo
"$MPI_INSTALL_DIR/bin/ompi_info" | grep -iE "Open MPI:|C compiler:|Fort compiler:" | head -3

log "Compiling + running a 2-rank MPI test…"
cat > /tmp/mpitest.c <<'EOF'
#include <mpi.h>
#include <stdio.h>
int main(int argc, char** argv){
  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  printf("  hello from rank %d of %d\n", rank, size);
  MPI_Finalize();
  return 0;
}
EOF
"$MPI_INSTALL_DIR/bin/mpicc" /tmp/mpitest.c -o /tmp/mpitest
"$MPI_INSTALL_DIR/bin/mpirun" -n 2 /tmp/mpitest

log "Step 2 complete — OpenMPI $VER built and working."
cat <<'MSG'

--------------------------------------------------------------------
NEXT:
  Open a new terminal (or `source ~/.bashrc`), then:
    conda activate adflow
    which mpicc          # -> ~/packages/openmpi-5.0.9/opt-gfortran/bin/mpicc
  Then run:  ./03_petsc.sh
--------------------------------------------------------------------
MSG
