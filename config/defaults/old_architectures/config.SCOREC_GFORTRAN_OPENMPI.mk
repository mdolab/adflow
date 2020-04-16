# ----------------------------------------------------------------------
# Config file for RPI SCOREC clusters using GCC 4.7.0 and OpenMPI
# ----------------------------------------------------------------------

# ------- Define a possible parallel make (use PMAKE = make otherwise)--
PMAKE = make -j 4

# ------- Define the MPI Compilers--------------------------------------
FF90 = mpif90
CC   = mpicc

# ------- Define Precision Flags ---------------------------------------
# Options for Integer precision flags: -DUSE_LONG_INT
# Options for Real precision flags: -DUSE_SINGLE_PRECISION, -DUSE_QUADRUPLE_PRECISION
# Default (nothing specified) is 4 byte integers and 8 byte reals

FF90_INTEGER_PRECISION_FLAG =
FF90_REAL_PRECISION_FLAG    = 
CC_INTEGER_PRECISION_FLAG   =
CC_REAL_PRECISION_FLAG      = 

# ------- Define CGNS Inlcude and linker flags -------------------------
CGNS_INCLUDE_FLAGS = -I${CGNS_INCLUDE_DIR}
CGNS_LINKER_FLAGS = -Wl,-rpath,${CGNS_LIB_DIR} -L${CGNS_LIB_DIR} -lcgns

# ------- Define Compiler Flags ----------------------------------------

FF90_FLAGS  = -DHAS_ISNAN -fPIC -fdefault-real-8 -O2 -fdefault-double-8 -g
C_FLAGS     = -DHAS_ISNAN -O2 -fPIC -g

# ------- Define Archiver  and Flags -----------------------------------
AR       = ar
AR_FLAGS = -rvs

# ------- Define Linker Flags ------------------------------------------
LINKER       = $(FF90)
LINKER_FLAGS = 

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
