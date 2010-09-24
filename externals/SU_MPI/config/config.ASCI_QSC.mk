#      ******************************************************************
#      *                                                                *
#      * File:          config.ASCI_QSC.mk                              *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 01-18-2005                                      *
#      * Last modified: 02-21-2006                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Defines the compiler settings and other commands  *
#      *              to have "make" function correctly. This file      *
#      *              defines the settings for a parallel executable    *
#      *              on the ASCI qsc machine in Los Alamos.            *
#      *                                                                *
#      ******************************************************************

#      ==================================================================

#      ******************************************************************
#      *                                                                *
#      * Possibly overrule the make command to allow for parallel make. *
#      *                                                                *
#      ******************************************************************

MAKE = gmake

#      ******************************************************************
#      *                                                                *
#      * F90 and C compiler definitions.                                *
#      *                                                                *
#      ******************************************************************

FF90 = f90
CC   = cc

#      ******************************************************************
#      *                                                                *
#      * Compiler flags. The MPI_COMPILE_FLAGS are set when the modules *
#      * are loaded.                                                    *
#      *                                                                *
#      ******************************************************************

COMMAND_SEARCH_PATH_MODULES = -I

FF90_GEN_FLAGS = $(MPI_COMPILE_FLAGS) -DUSE_MPI_INCLUDE_FILE
CC_GEN_FLAGS   = $(MPI_COMPILE_FLAGS)

FF90_OPTFLAGS  = -fast -transform_loops
CC_OPTFLAGS    = -O3

DEBUGFLAGS      = -check bounds -check format -check output_conversion \
		  -check overflow -check underflow
#FF90_DEBUGFLAGS = -g $(DEBUGFLAGS)
#FF90_DEBUGFLAGS = -g3 -O3
#FF90_DEBUGFLAGS = -g
#CC_DEBUGFLAGS   = -g $(DEBUGFLAGS)

FF90_FLAGS = $(FF90_GEN_FLAGS) $(FF90_OPTFLAGS) $(FF90_DEBUGFLAGS)
CC_FLAGS   = $(CC_GEN_FLAGS)   $(CC_OPTFLAGS)   $(CC_DEBUGFLAGS)

#      ******************************************************************
#      *                                                                *
#      * Archiver and archiver flags.                                   *
#      *                                                                *
#      ******************************************************************

AR       = ar
AR_FLAGS = -rv
