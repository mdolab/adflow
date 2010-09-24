#      ******************************************************************
#      *                                                                *
#      * File:          config.ASCI_QSC.mk                              *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 01-18-2005                                      *
#      * Last modified: 01-18-2005                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Defines the compiler settings and other commands  *
#      *              to have "make" function correctly. This file      *
#      *              defines the settings for the ASCI qsc machine in  *
#      *              Los Alamos.                                       *
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
#      * F90 definition.                                                *
#      *                                                                *
#      ******************************************************************

FF90 = f90

#      ******************************************************************
#      *                                                                *
#      * Precision flags. When nothing is specified 4 byte integer and  *
#      * and 8 byte real types are used. The flag -DUSE_LONG_INT sets   *
#      * the 8 byte integer to the standard integer type. The flag      *
#      * -DUSE_SINGLE_PRECISION sets the 4 byte real to the standard    *
#      * real type.                                                     *
#      *                                                                *
#      ******************************************************************

INTEGER_PRECISION_FLAG =# -DUSE_LONG_INT
REAL_PRECISION_FLAG    =# -DUSE_SINGLE_PRECISION

#      ******************************************************************
#      *                                                                *
#      * Compiler flags.                                                *
#      *                                                                *
#      ******************************************************************

COMMAND_SEARCH_PATH_MODULES = -I

FF90_GEN_FLAGS =

FF90_OPTFLAGS  = -fast -transform_loops

DEBUGFLAGS      = -check bounds -check format -check output_conversion \
		  -check overflow -check underflow
#FF90_DEBUGFLAGS = -g $(DEBUGFLAGS)
#FF90_DEBUGFLAGS = -g3 -O3
#FF90_DEBUGFLAGS = -g

FF90_FLAGS = $(FF90_GEN_FLAGS) $(FF90_OPTFLAGS) $(FF90_DEBUGFLAGS)

#      ******************************************************************
#      *                                                                *
#      * Archiver and archiver flags.                                   *
#      *                                                                *
#      ******************************************************************

AR       = ar
AR_FLAGS = -rv
