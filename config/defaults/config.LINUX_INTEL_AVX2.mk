# ----------------------------------------------------------------------
# Config file for Intel ifort
# ----------------------------------------------------------------------

# ------- Define a possible parallel make (use PMAKE = make otherwise)--
PMAKE = make -j 4

# ------- Define the MPI Compilers--------------------------------------
ifdef I_MPI_ROOT # Using Intel MPI
  # Note that ";" is there to avoid make shell optimization, otherwise the shell command may fail
  ICC_EXISTS := $(shell command -v icc;)

  ifdef ICC_EXISTS
    # icc only exists on older Intel versions
    # Assume that we want to use the old compilers
    FF90 = mpiifort
    CC = mpiicc
  else
    # Use the new compilers
    FF90 = mpiifx
    CC = mpiicx
  endif

else # Using HPE MPI
  FF90 = ifort -lmpi
  CC   = icc -lmpi
endif

# ------- Define Precision Flags ---------------------------------------
# Options for Integer precision flags: -DUSE_LONG_INT
# Options for Real precision flags: -DUSE_SINGLE_PRECISION, -DUSE_QUADRUPLE_PRECISION
# Default (nothing specified) is 4 byte integers and 8 byte reals

FF90_INTEGER_PRECISION_FLAG =
FF90_REAL_PRECISION_FLAG    =
CC_INTEGER_PRECISION_FLAG   =
CC_REAL_PRECISION_FLAG      =

# ------- Define CGNS Inlcude and linker flags -------------------------
# Define the CGNS include directory and linking flags for the CGNS library.
# We are assuming that HDF5 came from PETSc so it is included in ${PETSC_LIB}.
# Otherwise you will have to specify the HDF5 library.
CGNS_INCLUDE_FLAGS=-I$(CGNS_HOME)/include
CGNS_LINKER_FLAGS=-L$(CGNS_HOME)/lib -lcgns

# ------- Define complexify inlcude and linker flags -------------------------
COMPLEXIFY_INCLUDE_FLAGS=-I$(COMPLEXIFY_DIR)/include
COMPLEXIFY_LINKER_FLAGS=-L$(COMPLEXIFY_DIR)/lib -lcomplexify

# ------- Define Compiler Flags ----------------------------------------
FF77_FLAGS   = -fPIC -r8
FF90_FLAGS = $(FF77_FLAGS) -std08
FFXX_OPT_FLAGS = -O2 -xCORE-AVX2
C_FLAGS      = -fPIC -O -xCORE-AVX2

# ------- Define Archiver  and Flags -----------------------------------
AR       = ar
AR_FLAGS = -rvs

# ------- Define Linker Flags ------------------------------------------
LINKER       = $(FF90)
LINKER_FLAGS = -nofor-main

# ------- Define Petsc Info ---
include ${PETSC_DIR}/lib/petsc/conf/variables
PETSC_INCLUDE_FLAGS=${PETSC_CC_INCLUDES} -I$(PETSC_DIR)
PETSC_LINKER_FLAGS=${PETSC_LIB}

# Combine flags from above -- don't modify here
FF90_PRECISION_FLAGS = $(FF90_INTEGER_PRECISION_FLAG)$(FF90_REAL_PRECISION_FLAG)
CC_PRECISION_FLAGS   = $(CC_INTEGER_PRECISION_FLAG) $(CC_REAL_PRECISION_FLAG)

# Define potentially different python, python-config and f2py executables:
PYTHON = python
PYTHON-CONFIG = python3-config # use python-config for python 2
F2PY = f2py
