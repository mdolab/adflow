#      ******************************************************************
#      *                                                                *
#      * File:          config.FLASH_PG.mk                              *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 01-08-2007                                      *
#      * Last modified: 01-08-200.                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Defines the compiler settings and other commands  *
#      *              to have "make" function correctly. This file      *
#      *              defines the settings for a parallel Linux         *
#      *              executable on the Flash machine in LANL. The      *
#      *              Fortran compiler is pgf90 and the C-compiler gcc. *
#      *                                                                *
#      ******************************************************************

#      ==================================================================

#      ******************************************************************
#      *                                                                *
#      * Possibly overrule the make command to allow for parallel make. *
#      *                                                                *
#      ******************************************************************

MAKE = make -j 2

#      ******************************************************************
#      *                                                                *
#      * F90 and C compiler definitions.                                *
#      *                                                                *
#      ******************************************************************

FF90 = pgf90
CC   = gcc

#      ******************************************************************
#      *                                                                *
#      * Compiler flags.                                                *
#      *                                                                *
#      ******************************************************************

COMMAND_SEARCH_PATH_MODULES = -I

FF90_GEN_FLAGS = -DUSE_MPI_INCLUDE_FILE $(MPI_COMPILE_FLAGS)
CC_GEN_FLAGS   = $(MPI_COMPILE_FLAGS)

FF90_OPTFLAGS   = -fast -tp k8-64
CC_OPTFLAGS     = -O3 -fexpensive-optimizations -frerun-cse-after-loop \
		  -fthread-jumps -funroll-loops -finline-functions

#FF90_DEBUGFLAGS = -g -Mbounds -Mdclchk -DDEBUG_MODE
#CC_DEBUGFLAGS   = -g -Wall -pedantic -DDEBUG_MODE

FF90_FLAGS = $(FF90_GEN_FLAGS) $(FF90_OPTFLAGS) $(FF90_DEBUGFLAGS)
CC_FLAGS   = $(CC_GEN_FLAGS)   $(CC_OPTFLAGS)   $(CC_DEBUGFLAGS)

#      ******************************************************************
#      *                                                                *
#      * Archiver and archiver flags.                                   *
#      *                                                                *
#      ******************************************************************

AR       = ar
AR_FLAGS = -rvs
