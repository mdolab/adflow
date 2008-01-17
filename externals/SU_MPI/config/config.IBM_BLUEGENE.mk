#      ******************************************************************
#      *                                                                *
#      * File:          config.IBM_BLUEGENE.mk                          *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 07-04-2005                                      *
#      * Last modified: 02-21-2006                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Defines the compiler settings and other commands  *
#      *              to have "make" function correctly. This file      *
#      *              defines the settings for a parallel executable    *
#      *              IBM xlf95 and cc compiler for the BlueGene        *
#      *              machine of the San Diego Supercomputer Center.    *
#      *                                                                *
#      ******************************************************************

#      ==================================================================

#      ******************************************************************
#      *                                                                *
#      * Possibly overrule the make command to allow for parallel make. *
#      *                                                                *
#      ******************************************************************

#MAKE = make

#      ******************************************************************
#      *                                                                *
#      * F90 and C compiler definitions.                                *
#      *                                                                *
#      ******************************************************************

FF90 = mpxlf90
CC   = mpcc

#      ******************************************************************
#      *                                                                *
#      * Compiler flags.                                                *
#      *                                                                *
#      ******************************************************************

COMMAND_SEARCH_PATH_MODULES = -I

FF90_GEN_FLAGS = -qsuffix=f=f90 -qsuffix=cpp=F90 -qarch=440 -WF,-DUSE_MPI_INCLUDE_FILE
CC_GEN_FLAGS   = -qarch=440

FF90_OPTFLAGS   = -O3
CC_OPTFLAGS     = -O3

#FF90_DEBUGFLAGS =
#CC_DEBUGFLAGS   =

FF90_FLAGS = $(FF90_GEN_FLAGS) $(FF90_OPTFLAGS) $(FF90_DEBUGFLAGS)
CC_FLAGS   = $(CC_GEN_FLAGS)   $(CC_OPTFLAGS)   $(CC_DEBUGFLAGS)

#      ******************************************************************
#      *                                                                *
#      * Archiver and archiver flags.                                   *
#      *                                                                *
#      ******************************************************************

AR       = ar
AR_FLAGS = -rvs
