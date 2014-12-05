# ----------------------------------------------------------------------
# Config file for RPI SCOREC clusters using GCC 4.7.0 and OpenMPI
# ----------------------------------------------------------------------

# ------- Define a possible parallel make (use PMAKE = make otherwise)--
PMAKE = make -j 4

# ------- Define the MPI Compilers--------------------------------------
FF90 = mpif90
CC   = mpicc

#-------- Define Exec Suffix for executable ----------------------------
EXEC_SUFFIX =

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

FF90_GEN_FLAGS = -DHAS_ISNAN 
CC_GEN_FLAGS   = -DHAS_ISNAN  

FF90_OPT_FLAGS   = -fPIC -fdefault-real-8 -O2 -fdefault-double-8 -g
CC_OPT_FLAGS     = -O2 -fPIC -g

FF90_DEBUG_FLAGS = #-check bounds -check all
CC_DEBUG_FLAGS   = #-g -Wall -pedantic -DDEBUG_MODE

# ------- Define Archiver  and Flags -----------------------------------
AR       = ar
AR_FLAGS = -rvs

# ------- Defin pV3 and pvm3 linker flags ------------------------------
#PV3_FLAGS          = -DUSE_PV3
#PV3_LINKER_FLAGS   = -L/usr/local/pV3/clients/LINUX-INTEL -lpV3
#PVM3_LINKER_FLAGS  = -L/usr/local/pvm3/lib/LINUX -lgpvm3 -lpvm3
#PV3_INT_SRC_DIR    = src/pv3Interface

# ------- Define Linker Flags ------------------------------------------
LINKER       = $(FF90)
LINKER_FLAGS = 

# ------- Define Petsc Info --- Should not need to modify this -----
include ${PETSC_DIR}/conf/variables
PETSC_INCLUDE_FLAGS=${PETSC_CC_INCLUDES} -I$(PETSC_DIR)
PETSC_LINKER_FLAGS=${PETSC_LIB}

# Combine flags from above -- don't modify here
FF90_FLAGS = $(FF90_GEN_FLAGS) $(FF90_OPT_FLAGS) $(FF90_DEBUG_FLAGS)
CC_FLAGS   = $(CC_GEN_FLAGS)   $(CC_OPT_FLAGS)   $(CC_DEBUG_FLAGS)
FF90_PRECISION_FLAGS = $(FF90_INTEGER_PRECISION_FLAG)$(FF90_REAL_PRECISION_FLAG)
CC_PRECISION_FLAGS   = $(CC_INTEGER_PRECISION_FLAG) $(CC_REAL_PRECISION_FLAG)

# Define potentially different python, python-config and f2py executables:
PYTHON = python
PYTHON-CONFIG = python-config
F2PY = f2py