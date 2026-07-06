#!/usr/bin/env bash
# ======================================================================
#  01_setup.sh  —  MACH-Aero / ADflow bootstrap for a FRESH macOS machine
#  Step 1 of 10.
#
#  Sets up, from nothing:
#    * Xcode Command Line Tools   (git, make, clang, curl)
#    * Homebrew (+ wget)
#    * miniforge3                 (conda)
#    * conda env 'adflow'         (python 3.11 + clang/gfortran + OpenBLAS
#                                  + cmake/make/swig + numpy/scipy)
#    * the MACH environment block written into ~/.bashrc
#
#  Auto-detects Apple Silicon (arm64) vs Intel (x86_64). Safe to re-run.
#
#  USAGE:
#    chmod +x 01_setup.sh
#    ./01_setup.sh
#  Then open a NEW terminal and:  conda activate adflow
# ======================================================================
set -eo pipefail    # NOTE: no `-u` — conda's activation scripts use unset vars

log()  { printf '\n\033[1;34m==>\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[warn]\033[0m %s\n' "$*"; }
die()  { printf '\033[1;31m[error]\033[0m %s\n' "$*" >&2; exit 1; }

ARCH="$(uname -m)"
log "Detected architecture: $ARCH"
[[ "$ARCH" == "arm64" || "$ARCH" == "x86_64" ]] || die "Unsupported arch: $ARCH"

# ----------------------------------------------------------------------
# 1. Xcode Command Line Tools
# ----------------------------------------------------------------------
if ! xcode-select -p >/dev/null 2>&1; then
  log "Installing Xcode Command Line Tools — a GUI dialog will pop up."
  xcode-select --install || true
  die "Finish the Xcode CLT install in the pop-up, then RE-RUN this script."
fi
log "Xcode CLT present: $(xcode-select -p)"

# ----------------------------------------------------------------------
# 2. Homebrew (+ wget)
# ----------------------------------------------------------------------
if [[ "$ARCH" == "arm64" ]]; then BREW_PREFIX=/opt/homebrew; else BREW_PREFIX=/usr/local; fi
if ! command -v brew >/dev/null 2>&1 && [[ ! -x "$BREW_PREFIX/bin/brew" ]]; then
  log "Installing Homebrew…  (may ask for your macOS password)"
  NONINTERACTIVE=1 /bin/bash -c \
    "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
fi
eval "$("$BREW_PREFIX/bin/brew" shellenv bash)"
command -v wget >/dev/null 2>&1 || brew install wget || \
  warn "wget install failed — later steps can fall back to curl."
log "Homebrew ready: $(brew --version | head -1)"

# ----------------------------------------------------------------------
# 3. miniforge3  (conda)
# ----------------------------------------------------------------------
export MINIFORGE_HOME="$HOME/miniforge3"
if [[ ! -d "$MINIFORGE_HOME" ]]; then
  log "Installing miniforge3 to $MINIFORGE_HOME…"
  if [[ "$ARCH" == "arm64" ]]; then MF=Miniforge3-MacOSX-arm64.sh; else MF=Miniforge3-MacOSX-x86_64.sh; fi
  curl -fsSL -o "/tmp/$MF" \
    "https://github.com/conda-forge/miniforge/releases/latest/download/$MF"
  bash "/tmp/$MF" -b -p "$MINIFORGE_HOME"
fi
# make `conda` usable inside this script
source "$MINIFORGE_HOME/etc/profile.d/conda.sh"
log "conda ready: $(conda --version)"

# ----------------------------------------------------------------------
# 4. conda env 'adflow' + toolchain
# ----------------------------------------------------------------------
if ! conda env list | grep -qE '^adflow[[:space:]]'; then
  log "Creating conda env 'adflow' (python 3.11)…"
  conda create -y -n adflow python=3.11
fi
conda activate adflow

log "Installing compiler toolchain + build tools from conda-forge…"
conda install -y -c conda-forge \
  c-compiler cxx-compiler fortran-compiler \
  openblas libopenblas cmake make swig pkg-config

log "Installing pinned numpy / scipy (numpy<2 is required for f2py builds)…"
pip install "numpy==1.26.4" "scipy==1.15.3"

# ----------------------------------------------------------------------
# 5. Write the MACH environment block into ~/.bashrc  (idempotent)
# ----------------------------------------------------------------------
BASHRC="$HOME/.bashrc"; touch "$BASHRC"
BEGIN="# >>> MACH-Aero env >>>"; END="# <<< MACH-Aero env <<<"
if grep -qF "$BEGIN" "$BASHRC"; then
  log "Refreshing existing MACH block in ~/.bashrc"
  awk -v b="$BEGIN" -v e="$END" '$0==b{s=1} !s{print} $0==e{s=0}' \
    "$BASHRC" > "$BASHRC.tmp" && mv "$BASHRC.tmp" "$BASHRC"
fi
cat >> "$BASHRC" <<'EOF'
# >>> MACH-Aero env >>>
export MINIFORGE_HOME="$HOME/miniforge3"
[ -f "$MINIFORGE_HOME/etc/profile.d/conda.sh" ] && . "$MINIFORGE_HOME/etc/profile.d/conda.sh"
adflow-env() { conda activate adflow; }

export PACKAGES="$HOME/packages"

# OpenMPI
export MPI_INSTALL_DIR="$PACKAGES/openmpi-5.0.9/opt-gfortran"

# PETSc (real + complex)
export PETSC_DIR_REAL="$PACKAGES/petsc-3.21.0";     export PETSC_ARCH_REAL="real-debug"
export PETSC_DIR_COMPLEX="$PACKAGES/petsc_complex"; export PETSC_ARCH_COMPLEX="complex-debug"

# SLEPc (real + complex)
export SLEPC_DIR_REAL="$PACKAGES/slepc-3.21.0";     export SLEPC_ARCH_REAL="real-debug"
export SLEPC_DIR_COMPLEX="$PACKAGES/slepc_complex"; export SLEPC_ARCH_COMPLEX="complex-debug"

# default to real for tools expecting single values
export PETSC_DIR="$PETSC_DIR_REAL";  export PETSC_ARCH="$PETSC_ARCH_REAL"
export SLEPC_DIR="$SLEPC_DIR_REAL";  export SLEPC_ARCH="$SLEPC_ARCH_REAL"

# CGNS
export CGNS_HOME="$PACKAGES/CGNS-4.5.0/opt-gfortran"
export CGNS_INCLUDE_FLAGS="-I$CGNS_HOME/include"
export CGNS_LINKER_FLAGS="-L$CGNS_HOME/lib -lcgns"

# complexify
export COMPLEXIFY_DIR="$PACKAGES/complexify/opt-gfortran"
export COMPLEXIFY_INCLUDE_FLAGS="-I$COMPLEXIFY_DIR/include"
export COMPLEXIFY_LINKER_FLAGS="-L$COMPLEXIFY_DIR/lib -lcomplexify"

export PATH="$MPI_INSTALL_DIR/bin:$CGNS_HOME/bin:$PATH"
# <<< MACH-Aero env <<<
EOF
log "Wrote MACH block to ~/.bashrc"

# Make sure Terminal (a login shell) actually sources ~/.bashrc, and brew.
PROFILE="$HOME/.bash_profile"; touch "$PROFILE"
grep -qF '.bashrc' "$PROFILE" || \
  printf '\n[ -f "$HOME/.bashrc" ] && . "$HOME/.bashrc"\n' >> "$PROFILE"
grep -qF 'brew shellenv' "$PROFILE" || \
  printf '\neval "$(%s/bin/brew shellenv bash)"\n' "$BREW_PREFIX" >> "$PROFILE"

# ----------------------------------------------------------------------
# 6. Report
# ----------------------------------------------------------------------
log "Step 1 complete. Verification:"
echo "    conda    : $(command -v conda)"
echo "    python   : $(python --version 2>&1)  ->  $(command -v python)"
echo "    clang    : $(command -v clang)"
echo "    gfortran : $(command -v gfortran)"
echo "    cmake    : $(command -v cmake)"
cat <<'MSG'

--------------------------------------------------------------------
NEXT:
  1. Close this terminal and open a NEW one.
  2. If your prompt is zsh (default on macOS), start bash first:   bash -l
  3. Activate the env:                                             conda activate adflow
  4. Confirm you see the packages env, then run:                   ./02_openmpi.sh
--------------------------------------------------------------------
MSG
