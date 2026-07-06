# ----------------------------------------------------------------------
# ADflow config — macOS (Apple Silicon), REAL build only.
# No SLEPc, no complexify (this machine runs the non-eigenvalue ADflow).
# Ready to run:  cp this to  adflow/config/config.mk
# ----------------------------------------------------------------------

# ------- Define a possible parallel make (use PMAKE = make otherwise)--
PMAKE = make -j 4

# ------- Define the MPI Compilers--------------------------------------
FF90 = mpifort
CC   = mpicc

# ------- Define Precision Flags ---------------------------------------
FF90_INTEGER_PRECISION_FLAG =
FF90_REAL_PRECISION_FLAG    =
CC_INTEGER_PRECISION_FLAG   =
CC_REAL_PRECISION_FLAG      =

# ------- Define CGNS Inlcude and linker flags -------------------------
CGNS_INCLUDE_FLAGS=-I$(CGNS_HOME)/include
CGNS_LINKER_FLAGS=-Wl,-rpath,$(CGNS_HOME)/lib -L$(CGNS_HOME)/lib -lcgns

# ------- Define Compiler Flags ----------------------------------------
FF77_FLAGS = -fPIC -fdefault-real-8 -fdefault-double-8 -march=native
FF90_FLAGS = $(FF77_FLAGS) -std=f2008
FFXX_OPT_FLAGS = -O3 -ffast-math
C_FLAGS   = -fPIC -O

# ------- Define Archiver  and Flags -----------------------------------
AR       = ar
AR_FLAGS = -rvs

# ------- Define Linker Flags ------------------------------------------
LINKER       = $(FF90)
LINKER_FLAGS = -undefined dynamic_lookup

# ------- Define Petsc Info --- Should not need to modify this -----
include ${PETSC_DIR}/lib/petsc/conf/variables # PETSc 3.6
PETSC_INCLUDE_FLAGS=${PETSC_CC_INCLUDES} -I$(PETSC_DIR)
PETSC_LINKER_FLAGS=${PETSC_LIB}

# Combine flags from above -- don't modify here
FF90_PRECISION_FLAGS = $(FF90_INTEGER_PRECISION_FLAG)$(FF90_REAL_PRECISION_FLAG)
CC_PRECISION_FLAGS   = $(CC_INTEGER_PRECISION_FLAG) $(CC_REAL_PRECISION_FLAG)

# Define potentially different python, python-config and f2py executables:
PYTHON = python
# Pin to conda 3.11 (a bare python3-config can resolve to Homebrew python@3.14).
PYTHON-CONFIG = python3.11-config
F2PY = f2py
