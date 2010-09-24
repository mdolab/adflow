#      ******************************************************************
#      *                                                                *
#      * File:          config.ALTIX.mk                                 *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 01-18-2005                                      *
#      * Last modified: 02-21-2006                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Defines the compiler settings and other commands  *
#      *              to have "make" function correctly. This file      *
#      *              defines the settings for a sequential executable  *
#      *              on the SGI ALTIX (Itanium 2 processors) machine   *
#      *              using the intel compilers.                        *
#      *                                                                *
#      ******************************************************************

#      ==================================================================

#      ******************************************************************
#      *                                                                *
#      * Possibly overrule the make command to allow for parallel make. *
#      *                                                                *
#      ******************************************************************

#MAKE = make -j 8

#      ******************************************************************
#      *                                                                *
#      * F90 and C compiler definitions; efc and ecc is version 7.0,    *
#      * ifort and icc is version 8.0 of the intel compiler.            *
#      *                                                                *
#      ******************************************************************

FF90 = ifort
CC   = icc
#FF90 = efc
#CC   = ecc

#      ******************************************************************
#      *                                                                *
#      * Compiler flags. Note the addition of -DSEQUENTIAL_MODE.        *
#      *                                                                *
#      ******************************************************************

COMMAND_SEARCH_PATH_MODULES = -I

FF90_GEN_FLAGS = -fpic -DSEQUENTIAL_MODE
CC_GEN_FLAGS   = -fpic -DSEQUENTIAL_MODE

FF90_OPTFLAGS   = -O2
CC_OPTFLAGS     = -O2

#FF90_DEBUGFLAGS = -g -implicitnone -e90 -e95 -DDEBUG_MODE
#FF90_DEBUGFLAGS = -g -implicitnone -DDEBUG_MODE
#CC_DEBUGFLAGS   = -g -DDEBUG_MODE

FF90_FLAGS = $(FF90_GEN_FLAGS) $(FF90_OPTFLAGS) $(FF90_DEBUGFLAGS)
CC_FLAGS   = $(CC_GEN_FLAGS)   $(CC_OPTFLAGS)   $(CC_DEBUGFLAGS)

#      ******************************************************************
#      *                                                                *
#      * Archiver and archiver flags.                                   *
#      *                                                                *
#      ******************************************************************

AR       = ar
AR_FLAGS = -rvs
