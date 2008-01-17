#      ******************************************************************
#      *                                                                *
#      * File:          config.ALTIX.mk                                 *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 01-18-2005                                      *
#      * Last modified: 02-23-2006                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Defines the compiler settings and other commands  *
#      *              to have "make" function correctly. This file      *
#      *              defines the settings for the SGI ALTIX (Itanium 2 *
#      *              processors) machine using the efc/ifort compiler. *
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
#      * F90 definition; efc is version 7.0, ifort is version 8.1 of    *
#      * the intel compiler.                                            *
#      *                                                                *
#      ******************************************************************

FF90 = ifort
#FF90 = efc

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

FF90_GEN_FLAGS = -fpic

FF90_OPTFLAGS   = -O2

#FF90_DEBUGFLAGS = -g -implicitnone -e90 -e95 -DDEBUG_MODE
#FF90_DEBUGFLAGS = -g -implicitnone -DDEBUG_MODE

FF90_FLAGS = $(FF90_GEN_FLAGS) $(FF90_OPTFLAGS) $(FF90_DEBUGFLAGS)

#      ******************************************************************
#      *                                                                *
#      * Archiver and archiver flags.                                   *
#      *                                                                *
#      ******************************************************************

AR       = ar
AR_FLAGS = -rvs
