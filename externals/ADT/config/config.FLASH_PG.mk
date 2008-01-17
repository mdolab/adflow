#      ******************************************************************
#      *                                                                *
#      * File:          config.FLASH_PG.mk                              *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 01-08-2007                                      *
#      * Last modified: 01-08-2007                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Defines the compiler settings and other commands  *
#      *              to have "make" function correctly. This file      *
#      *              defines the settings for an executable on the     *
#      *              LANL Flash system using the pgf90 compiler.       *
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
#      * F90 compiler definition.                                       *
#      *                                                                *
#      ******************************************************************

FF90 = pgf90

#      ******************************************************************
#      *                                                                *
#      * Precision flags. When nothing is specified 4 byte integer and  *
#      * and 8 byte real types are used. The flag -DUSE_LONG_INT sets   *
#      * the 8 byte integer to the standard integer type. The flag      *
#      * -DUSE_SINGLE_PRECISION sets the 4 byte real to the standard    *
#      * real type.                                                     *
#      *                                                                *
#      ******************************************************************

INTEGER_PRECISION_FLAG = #-DUSE_LONG_INT
REAL_PRECISION_FLAG    = #-DUSE_SINGLE_PRECISION

#      ******************************************************************
#      *                                                                *
#      * Compiler flags.                                                *
#      *                                                                *
#      ******************************************************************

COMMAND_SEARCH_PATH_MODULES = -I

FF90_GEN_FLAGS =

FF90_GEN_FLAGS =

#FF90_OPTFLAGS   = -fast -tp k8-64

FF90_DEBUGFLAGS = -g -Mbounds -Mdclchk -DDEBUG_MODE

FF90_FLAGS = $(FF90_GEN_FLAGS) $(FF90_OPTFLAGS) $(FF90_DEBUGFLAGS)

#      ******************************************************************
#      *                                                                *
#      * Archiver and archiver flags.                                   *
#      *                                                                *
#      ******************************************************************

AR       = ar
AR_FLAGS = -rvs
