#      ******************************************************************
#      *                                                                *
#      * File:          config.ALTIX.mk                                 *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 04-25-2004                                      *
#      * Last modified: 02-23-2006                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Defines the compiler settings and other commands  *
#      *              to have "make" function correctly. This file      *
#      *              defines the settings for a sequential executable  *
#      *              on the SGI ALTIX (Itanium 2 processors) machine   *
#      *              using the efc and ecc compilers.                  *
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

EXEC_SUFFIX = _altix

#      ******************************************************************
#      *                                                                *
#      * F90 and C compiler definitions; efc and ecc is version 7.0,    *
#      * ifort and icc is version 8.0 of the intel compiler.            *
#      *                                                                *
#      ******************************************************************

#FF90 = ifort
#CC   = icc
FF90 = efc
CC   = ecc

#      ******************************************************************
#      *                                                                *
#      * CGNS include and linker flags.                                 *
#      *                                                                *
#      ******************************************************************

CGNS_INCLUDE_FLAGS = -I/usr/local/include
CGNS_LINKER_FLAGS  = -L/usr/local/lib -lcgns

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

FF90_GEN_FLAGS =
CC_GEN_FLAGS   =

FF90_OPTFLAGS   = -O2
CC_OPTFLAGS     = -O

#FF90_DEBUGFLAGS = -g -implicitnone -e90 -e95 -DDEBUG_MODE
#FF90_DEBUGFLAGS = -g -implicitnone -DDEBUG_MODE
#CC_DEBUGFLAGS   = -g -DDEBUG_MODE

FF90_FLAGS = $(FF90_GEN_FLAGS) $(FF90_OPTFLAGS) $(FF90_DEBUGFLAGS)
CC_FLAGS   = $(CC_GEN_FLAGS)   $(CC_OPTFLAGS)   $(CC_DEBUGFLAGS)

#      ******************************************************************
#      *                                                                *
#      * pV3 and pvm3 linker flags.                                     *
#      *                                                                *
#      ******************************************************************

#PV3_FLAGS          = -DUSE_PV3
#PV3_LINKER_FLAGS   = -L/opt/pV3/clients/LINUX64 -lpV3
#PVM3_LINKER_FLAGS  = #-L/opt/pvm3/lib/LINUX64 -lgpvm3 -lpvm3
#PV3_INT_SRC_DIR    = #src/pv3Interface

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

#      Intel compiler version 8.0

#LINKER       = $(FF90)
#LINKER_FLAGS = $(FF90_OPTFLAGS) -nofor_main

#      Intel compiler version 7.0

LINKER       = $(FF90)
LINKER_FLAGS = $(FF90_OPTFLAGS)
