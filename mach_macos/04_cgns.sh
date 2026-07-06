#!/usr/bin/env bash
# ======================================================================
#  04_cgns.sh  —  Build CGNS 4.5.0  (ADF format, no HDF5, Fortran bindings)
#  Requires Steps 1-3 (env + OpenMPI + PETSc).
#
#  Installs to:  $HOME/packages/CGNS-4.5.0/opt-gfortran
#  Compilers  :  conda clang (C) + gfortran (Fortran)   [serial, not MPI]
#
#  Self-contained:  chmod +x 04_cgns.sh && ./04_cgns.sh
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
export CGNS_HOME="${CGNS_HOME:-$PACKAGES/CGNS-4.5.0/opt-gfortran}"
VER=4.5.0
TARBALL="v${VER}.tar.gz"
URL="https://github.com/CGNS/CGNS/archive/refs/tags/${TARBALL}"
NPROC="$(sysctl -n hw.ncpu)"
mkdir -p "$PACKAGES"

# keep cmake<4 (PETSc's world was built with it; harmless backstop here too)
export CMAKE_POLICY_VERSION_MINIMUM=3.5

if [[ -f "$CGNS_HOME/lib/libcgns.dylib" ]]; then
  log "CGNS already installed at $CGNS_HOME — skipping."
else
  cd "$PACKAGES"
  if [[ ! -d "CGNS-${VER}" ]]; then
    [[ -f "$TARBALL" ]] || { log "Downloading CGNS $VER…"; curl -fL -o "$TARBALL" "$URL" || wget -O "$TARBALL" "$URL" || die "download failed."; }
    tar xzf "$TARBALL"
  fi
  cd "CGNS-${VER}"
  rm -f CMakeCache.txt
  log "Configuring CGNS (ADF, no HDF5, Fortran on)…"
  cmake \
    -D CGNS_ENABLE_FORTRAN=ON \
    -D CMAKE_INSTALL_PREFIX="$CGNS_HOME" \
    -D CGNS_ENABLE_64BIT=OFF \
    -D CGNS_ENABLE_HDF5=OFF \
    -D CGNS_BUILD_CGNSTOOLS=OFF \
    -D CMAKE_C_COMPILER="$(command -v clang)" \
    -D CMAKE_Fortran_COMPILER="$(command -v gfortran)" \
    -D CMAKE_C_FLAGS="-fPIC" \
    -D CMAKE_Fortran_FLAGS="-fPIC" .
  log "Building + installing…"
  make -j"$NPROC" install
fi

[[ -f "$CGNS_HOME/lib/libcgns.dylib" ]] || die "libcgns.dylib not found — build failed."
log "CGNS present: $CGNS_HOME/lib/libcgns.dylib"
ls "$CGNS_HOME/include/" | grep -i cgnslib >/dev/null && log "cgnslib.h header present."
log "Step 4 complete — CGNS built."
echo; echo "NEXT:  ./05_pywrappers.sh"
