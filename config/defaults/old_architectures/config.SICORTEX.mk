#      ******************************************************************
#      *                                                                *
#      * File:          config.SICORTEX.mk                              *
#      * Author:        Ruben E. Perez                                  *
#      * Starting date: 07-31-2009                                      *
#      * Last modified: 07-31-2009                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Defines the compiler settings and other commands  *
#      *              to have "make" function correctly. This file      *
#      *              defines the settings for a sequential executable  *
#      *              on the SICORTEX (MIPS64 processors) type machines *
#      *              using the pathscale compilers.                    *
#      *                                                                *
#      ******************************************************************

#      ==================================================================

#      ******************************************************************
#      *                                                                *
#      * Possibly overrule the make command to allow for parallel make. *
#      *                                                                *
#      ******************************************************************

MAKE = make -j 8

#      ******************************************************************
#      *                                                                *
#      * Suffix of the executable such that binaries of multiple        *
#      * platforms can be stored in the bin directory.                  *
#      *                                                                *
#      ******************************************************************

EXEC_SUFFIX = _sicortex

#      ******************************************************************
#      *                                                                *
#      * F90 and C compiler definitions.                                *
#      *                                                                *
#      ******************************************************************

FF90 = pathf95
CC   = pathcc

#      ******************************************************************
#      *                                                                *
#      * CGNS include and linker flags.                                 *
#      *                                                                *
#      ******************************************************************

#CGNS_INCLUDE_FLAGS = -DUSE_NO_CGNS
#CGNS_INCLUDE_FLAGS =
CGNS_INCLUDE_FLAGS = -I/usr/include
CGNS_LINKER_FLAGS  = -L/usr/lib64 -lcgns

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
#      * Compiler flags.                                                *
#      *                                                                *
#      ******************************************************************

COMMAND_SEARCH_PATH_MODULES = -I

FF90_GEN_FLAGS = -fPIC -fno-second-underscore -G6 -DHAS_ISNAN
CC_GEN_FLAGS   = -fPIC -fno-strict-aliasing

FF90_OPTFLAGS   = -O3 -OPT:Ofast -ffast-math
CC_OPTFLAGS     = -O3 -OPT:Ofast -ffast-math

#FF90_DEBUGFLAGS = -DDEBUG_MODE -g -Wall
#CC_DEBUGFLAGS   = -DDEBUG_MODE -g -Wall

FF90_FLAGS = $(FF90_GEN_FLAGS) $(FF90_OPTFLAGS) $(FF90_DEBUGFLAGS)
CC_FLAGS   = $(CC_GEN_FLAGS)   $(CC_OPTFLAGS)   $(CC_DEBUGFLAGS)

#      ******************************************************************
#      *                                                                *
#      * pV3 and pvm3 linker flags.                                     *
#      *                                                                *
#      ******************************************************************

#PV3_FLAGS         = -DUSE_PV3
#PV3_LINKER_FLAGS  = -L/usr/local/pV3/clients/LINUX -lpV3
#PVM3_LINKER_FLAGS = -L/usr/local/pvm3/lib/LINUX -lgpvm3 -lpvm3
#PV3_INT_SRC_DIR   = src/pv3Interface

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

#PETSC_DIR = /usr/share/petsc
#PETSC_ARCH = linux-mips-n64
#PETSC_INCLUDE_FLAGS = -DUSE_NO_PETSC
PETSC_INCLUDE_FLAGS = -I$(PETSC_DIR) -I$(PETSC_DIR)/bmake/$(PETSC_ARCH) -I$(PETSC_DIR)/include -I$(PETSC_DIR)/include/mpiuni
PETSC_LINKER_FLAGS  = -L/usr/lib64 -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc -lX11 -llapack -lf77blas -latlas -lscm -lm

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
LINKER_FLAGS = $(FF90_OPTFLAGS)
