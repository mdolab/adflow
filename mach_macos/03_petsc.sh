#!/usr/bin/env bash
# ======================================================================
#  03_petsc.sh  —  Build PETSc 3.21.0  (REAL, arch real-debug)
#  Step 3 of the sequence.  Requires Steps 1-2 (env + OpenMPI).
#
#  Installs to:  $HOME/packages/petsc-3.21.0   (arch: real-debug)
#  Downloads+builds: metis, parmetis, superlu_dist, mumps, scalapack, hypre
#  BLAS/LAPACK:  the conda OpenBLAS in the adflow env.
#  MPI:          the OpenMPI you built in Step 2 (via --with-mpi-dir).
#
#  Self-contained. Just:
#    chmod +x 03_petsc.sh
#    ./03_petsc.sh
#
#  ~20-40 min. Safe to re-run (skips if libpetsc is already built).
# ======================================================================
set -eo pipefail

log()  { printf '\n\033[1;34m==>\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[warn]\033[0m %s\n' "$*"; }
die()  { printf '\033[1;31m[error]\033[0m %s\n' "$*" >&2; exit 1; }

# --- activate the conda toolchain -------------------------------------
export MINIFORGE_HOME="${MINIFORGE_HOME:-$HOME/miniforge3}"
[[ -f "$MINIFORGE_HOME/etc/profile.d/conda.sh" ]] || die "run 01_setup.sh first."
source "$MINIFORGE_HOME/etc/profile.d/conda.sh"
conda activate adflow || die "conda env 'adflow' missing — run 01_setup.sh first."

# --- paths (self-contained defaults, match the MACH bashrc block) ------
export PACKAGES="${PACKAGES:-$HOME/packages}"
export MPI_INSTALL_DIR="${MPI_INSTALL_DIR:-$PACKAGES/openmpi-5.0.9/opt-gfortran}"
export PETSC_DIR="${PETSC_DIR_REAL:-$PACKAGES/petsc-3.21.0}"
export PETSC_ARCH="real-debug"
VER=3.21.0
TARBALL="petsc-${VER}.tar.gz"
URL="https://web.cels.anl.gov/projects/petsc/download/release-snapshots/${TARBALL}"
BLAS_LIB="$CONDA_PREFIX/lib/libopenblas.dylib"

[[ -x "$MPI_INSTALL_DIR/bin/mpicc" ]] || die "OpenMPI not found — run 02_openmpi.sh first."
[[ -f "$BLAS_LIB" ]] || die "OpenBLAS not found at $BLAS_LIB — is the adflow env active?"
export PATH="$MPI_INSTALL_DIR/bin:$PATH"; hash -r

log "Build plan:"
echo "    PETSC_DIR  = $PETSC_DIR"
echo "    PETSC_ARCH = $PETSC_ARCH"
echo "    MPI        = $MPI_INSTALL_DIR  (mpicc: $(command -v mpicc))"
echo "    BLAS       = $BLAS_LIB"

# --- skip if already built --------------------------------------------
if [[ -f "$PETSC_DIR/$PETSC_ARCH/lib/libpetsc.dylib" ]]; then
  log "PETSc real already built at $PETSC_DIR/$PETSC_ARCH — skipping configure/make."
else
  # --- fetch + unpack -------------------------------------------------
  cd "$PACKAGES"
  if [[ ! -d "petsc-${VER}" ]]; then
    if [[ ! -f "$TARBALL" ]]; then
      log "Downloading PETSc $VER…"
      curl -fL -o "$TARBALL" "$URL" || wget -O "$TARBALL" "$URL" || die "download failed."
    fi
    log "Unpacking…"
    tar xzf "$TARBALL"
  fi
  cd "petsc-${VER}"

  # IMPORTANT: do NOT export CC/CXX/FC here — PETSc uses the MPI wrappers
  # from --with-mpi-dir, which already wrap the conda clang/gfortran.
  log "Configuring PETSc (real-debug)…"
  ./configure \
    --PETSC_ARCH=real-debug \
    --with-scalar-type=real \
    --with-debugging=1 \
    --with-mpi-dir="$MPI_INSTALL_DIR" \
    --with-blaslapack-lib="$BLAS_LIB" \
    --with-shared-libraries=yes \
    --with-fortran-bindings=1 \
    --with-cxx-dialect=C++11 \
    --download-metis=yes \
    --download-parmetis=yes \
    --download-superlu_dist=yes \
    --download-mumps=yes \
    --download-scalapack=yes \
    --download-hypre=yes

  log "Building PETSc (this compiles metis/mumps/superlu_dist/hypre — slow)…"
  make PETSC_DIR="$PWD" PETSC_ARCH=real-debug all
fi

# --- verify -----------------------------------------------------------
[[ -f "$PETSC_DIR/$PETSC_ARCH/lib/libpetsc.dylib" ]] || die "libpetsc.dylib not found — build failed."
log "PETSc library present: $PETSC_DIR/$PETSC_ARCH/lib/libpetsc.dylib"

log "Running PETSc's own check (small MPI solve)…"
if make -C "$PETSC_DIR" PETSC_DIR="$PETSC_DIR" PETSC_ARCH=real-debug check; then
  log "PETSc check passed."
else
  warn "PETSc 'make check' reported an issue — the library built OK, we'll"
  warn "confirm functionality later via petsc4py/adflow. Note it and continue."
fi

log "Step 3 complete — real PETSc built."
cat <<'MSG'

--------------------------------------------------------------------
NEXT:  ./04_petsc_complex.sh   (complex PETSc, for the LST/eigen build)
       Similar length; builds into ~/packages/petsc_complex.
--------------------------------------------------------------------
MSG
