#      ******************************************************************
#      *                                                                *
#      * File:          config.SGI_N32_MPI_ORIGIN.mk                    *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 12-03-2003                                      *
#      * Last modified: 02-23-2006                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Defines the compiler settings and other commands  *
#      *              to have "make" function correctly. This file      *
#      *              defines the settings for a parallel N32 bit SGI   *
#      *              executable in combination with SGI MPI, i.e. a    *
#      *              shared memory MPI version. This will typically    *
#      *              run on machines of the Origin class.              *
#      *                                                                *
#      ******************************************************************

#      ==================================================================

#      ******************************************************************
#      *                                                                *
#      * Possibly overrule the make command to allow for parallel make. *
#      *                                                                *
#      ******************************************************************

MAKE = gmake -j 8

#      ******************************************************************
#      *                                                                *
#      * Suffix of the executable such that binaries of multiple        *
#      * platforms can be stored in the bin directory.                  *
#      *                                                                *
#      ******************************************************************

EXEC_SUFFIX = _sgi_n32_origin_mpi

#      ******************************************************************
#      *                                                                *
#      * F90 and C compiler definitions.                                *
#      *                                                                *
#      ******************************************************************

FF90 = f90
CC   = cc

#      ******************************************************************
#      *                                                                *
#      * CGNS include and linker flags.                                 *
#      *                                                                *
#      ******************************************************************

CGNS_INCLUDE_FLAGS = -I/usr/local/CGNS2.3/include
CGNS_LINKER_FLAGS  = -L/usr/local/CGNS2.3/libN32 -lcgns

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
#      * The Ofast option is based on the 600 MHZ IP35 processor.       *
#      *                                                                *
#      ******************************************************************

COMMAND_SEARCH_PATH_MODULES = -I

FF90_32_64_FLAGS = -mips4 -n32
CC_32_64_FLAGS   = -mips4 -n32

FF90_GEN_FLAGS = $(FF90_32_64_FLAGS)
CC_GEN_FLAGS   = $(CC_32_64_FLAGS)

#FF90_OPTFLAGS   = -r14000 -Ofast=ip35 -IPA
#CC_OPTFLAGS     = -r14000 -Ofast=ip35 -IPA

#FF90_OPTFLAGS   = -O2
#CC_OPTFLAGS     = -O2

FF90_OPTFLAGS   = -O3
CC_OPTFLAGS     = -O

DEBUGFLAGS      = -DEBUG:conform_check=ON:div_check=3:subscript_check=ON:trap_uninitialized=ON:varargs_interface_check=ON:verbose_runtime=ON
#FF90_DEBUGFLAGS = -g -fullwarn -check_bounds $(DEBUGFLAGS) -DDEBUG_MODE
#FF90_DEBUGFLAGS = -g3 -O3 -DDEBUG_MODE
#FF90_DEBUGFLAGS = -g -DDEBUG_MODE
#CC_DEBUGFLAGS   = -g -fullwarn $(DEBUGFLAGS) -DDEBUG_MODE
#CC_DEBUGFLAGS   = -g -DDEBUG_MODE

FF90_FLAGS = $(FF90_GEN_FLAGS) $(FF90_OPTFLAGS) $(FF90_DEBUGFLAGS)
CC_FLAGS   = $(CC_GEN_FLAGS)   $(CC_OPTFLAGS)   $(CC_DEBUGFLAGS)

#      ******************************************************************
#      *                                                                *
#      * pV3 and pvm3 linker flags.                                     *
#      *                                                                *
#      ******************************************************************

#PV3_FLAGS          = -DUSE_PV3
#PV3_LINKER_FLAGS   = -L/usr/local/pV3/clients/SGI6 -lpV3
#PVM3_LINKER_FLAGS  = -L/usr/local/pvm3/lib/SGI6 -lgpvm3 -lpvm3
#PV3_INT_SRC_DIR    = src/pv3Interface

#      ******************************************************************
#      *                                                                *
#      * Archiver and archiver flags.                                   *
#      *                                                                *
#      ******************************************************************

AR       = ar
AR_FLAGS = -rv

#      ******************************************************************
#      *                                                                *
#      * Linker and linker flags.                                       *
#      *                                                                *
#      ******************************************************************

LINKER       = $(FF90)
LINKER_FLAGS = $(FF90_32_64_FLAGS) $(FF90_OPTFLAGS) -lmpi
