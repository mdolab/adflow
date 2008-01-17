#      ******************************************************************
#      *                                                                *
#      * File:          config.IBM_DATASTAR.mk                          *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 07-04-2005                                      *
#      * Last modified: 03-30-2006                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Defines the compiler settings and other commands  *
#      *              to have "make" function correctly. This file      *
#      *              defines the settings for a parallel executable    *
#      *              IBM xlf95 and cc compiler for the DataStar        *
#      *              machine of the San Diego Supercomputer Center.    *
#      *                                                                *
#      ******************************************************************

#      ==================================================================

#      ******************************************************************
#      *                                                                *
#      * Possibly overrule the make command to allow for parallel make. *
#      *                                                                *
#      ******************************************************************

#MAKE = gmake

#      ******************************************************************
#      *                                                                *
#      * F90 and C compiler definitions.                                *
#      *                                                                *
#      ******************************************************************

FF90 = mpxlf95_r
CC   = mpcc_r

#      ******************************************************************
#      *                                                                *
#      * Compiler flags.                                                *
#      *                                                                *
#      ******************************************************************

COMMAND_SEARCH_PATH_MODULES = -I

FF90_GEN_FLAGS = -qsuffix=f=f90 -qsuffix=cpp=F90 -qextname -q64
CC_GEN_FLAGS   = -q64

FF90_OPTFLAGS   = -qtune=pwr4 -qarch=pwr4 -O3
CC_OPTFLAGS     = -qtune=pwr4 -qarch=pwr4 -O3

#FF90_DEBUGFLAGS = -qsigtrap -g -C
#CC_DEBUGFLAGS   = -qsigtrap -g

FF90_FLAGS = $(FF90_GEN_FLAGS) $(FF90_OPTFLAGS) $(FF90_DEBUGFLAGS)
CC_FLAGS   = $(CC_GEN_FLAGS)   $(CC_OPTFLAGS)   $(CC_DEBUGFLAGS)

#      ******************************************************************
#      *                                                                *
#      * Archiver and archiver flags.                                   *
#      *                                                                *
#      ******************************************************************

AR       = ar -X64
AR_FLAGS = -rvs
