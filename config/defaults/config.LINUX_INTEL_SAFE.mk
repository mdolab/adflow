# ----------------------------------------------------------------------
# Config file for Intel ifort -- safe version. ifort is TERRIBLE for
# generating WRONG code with O2+. There is a 50-50 chance O2 code is
# simply wrong. O2 and above MUST BE ASSUMED WRONG until it can be
# proven otherwise. This is a safe config file at O1 which does not
# do vectorization which will generally give correct code.
# ----------------------------------------------------------------------

# ------- Define a possible parallel make (use PMAKE = make otherwise)--
PMAKE = make -j 4

# ------- Define the MPI Compilers--------------------------------------
FF90 = mpiifort
CC   = mpiicc

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

# ------- Define Compiler Flags ----------------------------------------
FF77_FLAGS   = -DHAS_ISNAN -fPIC -r8 -O1 -g
FF90_FLAGS = $(FF77_FLAGS) -std08
C_FLAGS      = -DHAS_ISNAN -O -fPIC

# ------- Define Archiver  and Flags -----------------------------------
AR       = ar
AR_FLAGS = -rvs

# ------- Define Linker Flags ------------------------------------------
LINKER       = $(FF90)
LINKER_FLAGS = -nofor_main

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
