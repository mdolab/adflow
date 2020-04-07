# ----------------------------------------------------------------------
# Config file for Intel ifort  with OpenMPI
# ----------------------------------------------------------------------

# ------- Define a possible parallel make (use PMAKE = make otherwise)--
PMAKE = make -j 4

# ------- Define the MPI Compilers--------------------------------------
FF90 = mpif90
CC   = mpicc

#-------- Define Exec Suffix for executable ----------------------------
EXEC_SUFFIX = _linux_intel_openmpi_nyx

# ------- Define Precision Flags ---------------------------------------
# Options for Integer precision flags: -DUSE_LONG_INT
# Options for Real precision flags: -DUSE_SINGLE_PRECISION, -DUSE_QUADRUPLE_PRECISION
# Default (nothing specified) is 4 byte integers and 8 byte reals

FF90_INTEGER_PRECISION_FLAG =
FF90_REAL_PRECISION_FLAG    = 
CC_INTEGER_PRECISION_FLAG   =
CC_REAL_PRECISION_FLAG      = 

# ------- Define CGNS Inlcude and linker flags -------------------------
CGNS_INCLUDE_FLAGS = -I/nobackup/kenway/opt/include 
CGNS_LINKER_FLAGS  = -L/nobackup/kenway/opt/lib/ -lcgns 

# ------- Define Compiler Flags ----------------------------------------

FF90_GEN_FLAGS = -DHAS_ISNAN 
CC_GEN_FLAGS   = -DHAS_ISNAN  

FF90_OPT_FLAGS   =  -fPIC -r8 -O1  -xHost #fp-model precise
CC_OPT_FLAGS     = -O -fPIC

FF90_DEBUG_FLAGS = #-check bounds -check all
CC_DEBUG_FLAGS   = #-g -Wall -pedantic -DDEBUG_MODE

# ------- Define Archiver  and Flags -----------------------------------
AR       = ar
AR_FLAGS = -rvs

# ------- Define Linker Flags ------------------------------------------
LINKER       = $(FF90)
LINKER_FLAGS = -nofor_main

# ------- Define Petsc Info --- Should not need to modify this -----
include ${PETSC_DIR}/lib/petsc/conf/variables # PETSc 3.6
#include ${PETSC_DIR}/conf/variables # PETSc 3.5
PETSC_INCLUDE_FLAGS=${PETSC_CC_INCLUDES} -I$(PETSC_DIR)
PETSC_LINKER_FLAGS=${PETSC_LIB}

# Combine flags from above -- don't modify here
FF90_FLAGS = $(FF90_GEN_FLAGS) $(FF90_OPT_FLAGS) $(FF90_DEBUG_FLAGS)
CC_FLAGS   = $(CC_GEN_FLAGS)   $(CC_OPT_FLAGS)   $(CC_DEBUG_FLAGS)
FF90_PRECISION_FLAGS = $(FF90_INTEGER_PRECISION_FLAG)$(FF90_REAL_PRECISION_FLAG)
CC_PRECISION_FLAGS   = $(CC_INTEGER_PRECISION_FLAG) $(CC_REAL_PRECISION_FLAG)

# Define potentially different python, python-config and f2py executables:
PYTHON = python
PYTHON-CONFIG = python3-config # use python-config for python 2
F2PY = f2py
