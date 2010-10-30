#      ******************************************************************
#      *                                                                *
#      * File:          config.SGI_N32_MPICH.mk                         *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 01-18-2005                                      *
#      * Last modified: 02-21-2006                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Defines the compiler settings and other commands  *
#      *              to have "make" function correctly. This file      *
#      *              defines the settings for a parallel N32 bit SGI   *
#      *              executable in combination with MPICH.             *
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
#      * F90 and C compiler definitions.                                *
#      *                                                                *
#      ******************************************************************

FF90 = /usr/local/MPICH/N32/bin/mpif90
CC   = /usr/local/MPICH/N32/bin/mpicc

#      ******************************************************************
#      *                                                                *
#      * Compiler flags.                                                *
#      * The Ofast option is based on the 600 MHZ IP35 processor.       *
#      *                                                                *
#      ******************************************************************

COMMAND_SEARCH_PATH_MODULES = -I

FF90_32_64_FLAGS = -mips4 -n32
CC_32_64_FLAGS   = -mips4 -n32

FF90_GEN_FLAGS = $(FF90_32_64_FLAGS) -DUSE_MPI_INCLUDE_FILE
CC_GEN_FLAGS   = $(CC_32_64_FLAGS)

#FF90_OPTFLAGS   = -r14000 -Ofast=ip35 -IPA
#CC_OPTFLAGS     = -r14000 -Ofast=ip35 -IPA

#FF90_OPTFLAGS   = -O2
#CC_OPTFLAGS     = -O2

FF90_OPTFLAGS   = -O3
CC_OPTFLAGS     = -O3

DEBUGFLAGS      = -DEBUG:conform_check=ON:div_check=3:subscript_check=ON:trap_uninitialized=ON:varargs_interface_check=ON:verbose_runtime=ON
#FF90_DEBUGFLAGS = -g -fullwarn -check_bounds $(DEBUGFLAGS) -DDEBUG_MODE
#FF90_DEBUGFLAGS = -g3 -O3 -DDEBUG_MODE
#FF90_DEBUGFLAGS = -g -DDEBUG_MODE
#CC_DEBUGFLAGS   = -g -fullwarn $(DEBUGFLAGS) -DDEBUG_MODE

FF90_FLAGS = $(FF90_GEN_FLAGS) $(FF90_OPTFLAGS) $(FF90_DEBUGFLAGS)
CC_FLAGS   = $(CC_GEN_FLAGS)   $(CC_OPTFLAGS)   $(CC_DEBUGFLAGS)

#      ******************************************************************
#      *                                                                *
#      * Archiver and archiver flags.                                   *
#      *                                                                *
#      ******************************************************************

AR       = ar
AR_FLAGS = -rv
