#      ******************************************************************
#      *                                                                *
#      * File:          config.LINUX_INTEL_MPICH.mk                     *
#      * Author:        Edwin van der Weide, Andre C. Marta             *
#      * Starting date: 12-07-2002                                      *
#      * Last modified: 11-29-2006                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Defines the compiler settings and other commands  *
#      *              to have "make" function correctly. This file      *
#      *              defines the settings for a parallel Linux         *
#      *              executable in combination with MPICH. Assumed is  *
#      *              that mpif90 and mpicc are based on the intel and  *
#      *              gcc compilers respectively.                       *
#      *                                                                *
#      ******************************************************************

#      ==================================================================

#      ******************************************************************
#      *                                                                *
#      * Possibly overrule the make command to allow for parallel make. *
#      *                                                                *
#      ******************************************************************

SHELL      = /bin/bash

RM         = /bin/rm -rf
MV         = /bin/mv -f
SYM_LINK   = ln -sf
#MAKE       = make
MAKE_CLEAN = make
MAKE = make -j 8

#      ******************************************************************
#      *                                                                *
#      * Suffix of the executable such that binaries of multiple        *
#      * platforms can be stored in the bin directory.                  *
#      *                                                                *
#      ******************************************************************

EXEC_SUFFIX = _linux_intel_openmpi_scinet

#      ******************************************************************
#      *                                                                *
#      * F90 and C compiler definitions.                                *
#      *                                                                *
#      ******************************************************************

#FF90 = /usr/local/mpich-intel/bin/mpif90
#CC   = /usr/local/mpich-intel/bin/mpicc
FF90 = mpif90 -fPIC
CC   = mpicc -fPIC


#      ******************************************************************
#      *                                                                *
#      * CGNS include and linker flags.                                 *
#      *                                                                *
#      ******************************************************************
CGNS_INCLUDE_FLAGS = -I$(HOME)/packages/cgnslib_2.5
CGNS_LINKER_FLAGS  = -L$(HOME)/packages/cgnslib_2.5/LINUX -lcgns
#CGNS_INCLUDE_FLAGS = -I$(HOME)/packages/cgnslib_2.5_debug
#CGNS_LINKER_FLAGS  = -L$(HOME)/packages/cgnslib_2.5_debug/LINUX -lcgns
#CGNS_INCLUDE_FLAGS = -I$(HOME)/packages/cgnslib_2.4
#CGNS_LINKER_FLAGS  = -L$(HOME)/packages/cgnslib_2.4/LINUX -lcgns
#CGNS_INCLUDE_FLAGS = -I/usr/local/include
#CGNS_LINKER_FLAGS  = -L/usr/local/lib64 -lcgns.intel
#CGNS_LINKER_FLAGS  = -L/usr/local/lib -lcgns.intel

#      ******************************************************************
#      *                                                                *
#      * Precision flags. When nothing is specified 4 byte integer and  *
#      * and 8 byte real types are used. The flag -DUSE_LONG_INT sets   *
#      * the 8 byte integer to the standard integer type. The flag      *
#      * -DUSE_SINGLE_PRECISION sets the 4 byte real to the standard    *
#      * real type and -DUSE_QUADRUPLE_PRECISION will make 16 byte      *
#      * reals the default type. The latter option may not work on all  *
#      * platforms.                                                     *
#      *                                                                *
#      ******************************************************************

#FF90_INTEGER_PRECISION_FLAG = -DUSE_LONG_INT
#CC_INTEGER_PRECISION_FLAG   = -DUSE_LONG_INT

#FF90_REAL_PRECISION_FLAG = -DUSE_SINGLE_PRECISION
#FF90_REAL_PRECISION_FLAG = -DUSE_QUADRUPLE_PRECISION
#CC_REAL_PRECISION_FLAG   = -DUSE_SINGLE_PRECISION
#CC_REAL_PRECISION_FLAG   = -DUSE_QUADRUPLE_PRECISION

FF90_PRECISION_FLAGS = $(FF90_INTEGER_PRECISION_FLAG) \
		       $(FF90_REAL_PRECISION_FLAG)
CC_PRECISION_FLAGS   = $(CC_INTEGER_PRECISION_FLAG) \
		       $(CC_REAL_PRECISION_FLAG)

#      ******************************************************************
#      *                                                                *
#      * Compiler flags. It is assumed that mpif90 is based on the      *
#      * intel compiler ifort and mpicc on the gcc gnu compiler.        *
#      *                                                                *
#      ******************************************************************

COMMAND_SEARCH_PATH_MODULES = -I

FF90_GEN_FLAGS = -DHAS_ISNAN -DUSE_PETSC_3
CC_GEN_FLAGS   =

#FF90_OPTFLAGS   = -O3 -ipo -ipo_obj
FF90_OPTFLAGS   = -r8 -O1 #-g -check all#-O1 #-r8#-O1 -r8 #-check all# -O2 -r8 -fpic# -check all #-tpp7 -xW -unroll -ip
#CC_OPTFLAGS     = -O3 -fexpensive-optimizations -frerun-cse-after-loop \
#		  -fthread-jumps -funroll-loops -finline-functions
CC_OPTFLAGS     = -O -g

#FF90_DEBUGFLAGS = -g -C -implicitnone -ftrapuv -debug extended \
#		  -traceback -DDEBUG_MODE
#FF90_DEBUGFLAGS = -g -implicitnone -DDEBUG_MODE
#CC_DEBUGFLAGS   = -g -Wall -pedantic -DDEBUG_MODE

FF90_FLAGS = $(FF90_GEN_FLAGS) $(FF90_OPTFLAGS) $(FF90_DEBUGFLAGS)
CC_FLAGS   = $(CC_GEN_FLAGS)   $(CC_OPTFLAGS)   $(CC_DEBUGFLAGS)

#      ******************************************************************
#      *                                                                *
#      * pV3 and pvm3 linker flags.                                     *
#      *                                                                *
#      ******************************************************************

#PV3_FLAGS          = -DUSE_PV3
#PV3_LINKER_FLAGS   = -L/usr/local/pV3/clients/LINUX-INTEL -lpV3
#PVM3_LINKER_FLAGS  = -L/usr/local/pvm3/lib/LINUX -lgpvm3 -lpvm3
#PV3_INT_SRC_DIR    = src/pv3Interface

#      ******************************************************************
#      *                                                                *
#      * PETSc include and linker flags. The flag -DUSE_NO_PETSC        *
#      * indicates that the PETSc libraries should not be included in   *
#      * the built process. These are only needed when running ADjoint  *
#      * problems, that use PETSc to solve the system of equations.     *
#      * If -DUSE_NO_PETSC is uncommented then the PETSc include and    *
#      * linker flags should be commented out.                          *
#      *                                                                *
#      ******************************************************************

#PETSC_ROOT_DIR      = /usr/local/petsc-2.3.1-p15
#PETSC_INCLUDE_FLAGS = -DUSE_NO_PETSC
#PETSC_INCLUDE_FLAGS = -I$(PETSC_ROOT_DIR) \
#		      -I$(PETSC_ROOT_DIR)/bmake/<platform> \
#		      -I$(PETSC_ROOT_DIR)/include
#PETSC_LINKER_FLAGS  = -L$(PETSC_ROOT_DIR)/lib/<platform> \
#		      -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc

X11_DIR = /usr/lib64

include ${PETSC_DIR}/conf/variables

PETSC_INCLUDE_FLAGS=${PETSC_CC_INCLUDES} -I$(PETSC_DIR)
PETSC_LINKER_FLAGS=${PETSC_LIB}

# PETSC_DIR = /home/kenway/packages/petsc-3.0.0-p10
# PETSC_ARCH = linux-gnu-c-debug
# PETSC_INCLUDE_FLAGS = -I$(PETSC_DIR) -I$(PETSC_DIR)/bmake/$(PETSC_ARCH) -I$(PETSC_DIR)/include -I$(PETSC_DIR)/include/mpiuni -I$(PETSC_DIR)/$(PETSC_ARCH)/include
# PETSC_LINKER_FLAGS  = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc -lpetsccontrib -lpetscsnes -lpetscts -L$(X11_DIR) -lX11 -L$(PETSC_DIR)/externalpackages/fblaslapack/$(PETSC_ARCH) -lflapack -lfblas -lm

#      ******************************************************************
#      *                                                                *
#      * Python Directories for includes                                *
#      *                                                                *
#      ******************************************************************

PYTHON_INCLUDE_DIR = /scinet/gpc/tools/Python/Python262/include/python2.6
NUMPY_ROOT_DIR = /scinet/gpc/tools/Python/Python262/lib/python2.6/site-packages/numpy-1.3.0-py2.6-linux-x86_64.egg

#      ******************************************************************
#      *                                                                *
#      * Archiver and archiver flags.                                   *
#      *                                                                *
#      ******************************************************************

AR       = ar
AR_FLAGS = -rvs

#      ******************************************************************
#      *                                                                *
#      * Linker and linker flags.                                       *
#      *                                                                *
#      ******************************************************************

LINKER       = $(FF90)
LINKER_FLAGS = $(FF90_OPTFLAGS) -nofor_main -lmkl
