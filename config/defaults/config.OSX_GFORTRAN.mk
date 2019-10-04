# ----------------------------------------------------------------------
# Config file for Intel ifort  with OpenMPI
# ----------------------------------------------------------------------

# ------- Define a possible parallel make (use PMAKE = make otherwise)--
PMAKE = make -j 4

# ------- Define the MPI Compilers--------------------------------------
FF90 = mpif90
CC   = mpicc

#-------- Define Exec Suffix for executable ----------------------------
EXEC_SUFFIX = _linux_gfortran_openmpi

# ------- Define Precision Flags ---------------------------------------
# Options for Integer precision flags: -DUSE_LONG_INT
# Options for Real precision flags: -DUSE_SINGLE_PRECISION, -DUSE_QUADRUPLE_PRECISION
# Default (nothing specified) is 4 byte integers and 8 byte reals

FF90_INTEGER_PRECISION_FLAG =
FF90_REAL_PRECISION_FLAG    =
CC_INTEGER_PRECISION_FLAG   =
CC_REAL_PRECISION_FLAG      =

# ------- Define CGNS Inlcude and linker flags -------------------------
# Define the CNGS include directory and linking flags for CGNSlib. We
# can use 3.2.x OR CGNS 3.3+. You must define which version is being
# employed as shown below. We are assuming that HDF5 came from PETSc
# so it is included in ${PETSC_LIB}. Otherwise you will have to
# specify the HDF5 library.

# ----------- CGNS 3.2.x ------------------
CGNS_VERSION_FLAG=
CGNS_INCLUDE_FLAGS=-I$(CGNS_HOME)/include
CGNS_LINKER_FLAGS=-L$(CGNS_HOME)/lib -lcgns

# # ----------- CGNS 3.3.x ------------------
# CGNS_VERSION_FLAG=-DUSECGNSMODULE
# CGNS_INCLUDE_FLAGS=-I$(HOME)/packages/CGNS/src
# CGNS_LINKER_FLAGS=-L$(HOME)/packages/CGNS/src/lib -lcgns

# ------- Define Compiler Flags ----------------------------------------
FF90_FLAGS = -DHAS_ISNAN  -fPIC -fdefault-real-8 -fdefault-double-8 -g  -O3 -march=native -ffast-math
C_FLAGS   = -DHAS_ISNAN  -O -fPIC -g

# ------- Define Archiver  and Flags -----------------------------------
AR       = ar
AR_FLAGS = -rvs

# ------- Define Linker Flags ------------------------------------------
LINKER       = $(FF90)
LINKER_FLAGS =  -undefined dynamic_lookup

# ------- Define Petsc Info --- Should not need to modify this -----
include ${PETSC_DIR}/lib/petsc/conf/variables # PETSc 3.6
#include ${PETSC_DIR}/conf/variables # PETSc 3.5
PETSC_INCLUDE_FLAGS=${PETSC_CC_INCLUDES} -I$(PETSC_DIR)
PETSC_LINKER_FLAGS=${PETSC_LIB}

# # Combine flags from above -- don't modify here
# FF90_FLAGS = $(FF90_GEN_FLAGS) $(FF90_OPT_FLAGS) $(FF90_DEBUG_FLAGS)
# CC_FLAGS   = $(CC_GEN_FLAGS)   $(CC_OPT_FLAGS)   $(CC_DEBUG_FLAGS)
FF90_PRECISION_FLAGS = $(FF90_INTEGER_PRECISION_FLAG)$(FF90_REAL_PRECISION_FLAG)
CC_PRECISION_FLAGS   = $(CC_INTEGER_PRECISION_FLAG) $(CC_REAL_PRECISION_FLAG)

# Define potentially different python, python-config and f2py executables:
PYTHON = python
PYTHON-CONFIG = python-config
F2PY = f2py
