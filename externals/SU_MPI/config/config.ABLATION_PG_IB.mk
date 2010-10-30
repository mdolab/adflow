#      ******************************************************************
#      *                                                                *
#      * File:          config.ABLATION_PG_IB.mk                        *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 08-23-2005                                      *
#      * Last modified: 02-21-2006                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Defines the compiler settings and other commands  *
#      *              to have "make" function correctly. This file      *
#      *              defines the settings for a parallel execcutable   *
#      *              for the Portland Group compilers on the Stanford  *
#      *              Linux cluster ablation.                           *
#      *                                                                *
#      ******************************************************************

#      ==================================================================

#      ******************************************************************
#      *                                                                *
#      * Possibly overrule the make command to allow for parallel make. *
#      *                                                                *
#      ******************************************************************

#MAKE = make -j 2

#      ******************************************************************
#      *                                                                *
#      * F90 and C compiler definitions.                                *
#      *                                                                *
#      ******************************************************************

FF90 = /usr/local/apps/mvapich/pgi/bin/mpif90
CC   = /usr/local/apps/mvapich/pgi/bin/mpicc

#      ******************************************************************
#      *                                                                *
#      * Compiler flags. It is assumed that mpif90 is based on the      *
#      * Portland Group compiler pgf90 and mpicc on the pgcc compiler.  *
#      *                                                                *
#      ******************************************************************

COMMAND_SEARCH_PATH_MODULES = -I

FF90_GEN_FLAGS = -DUSE_MPI_INCLUDE_FILE -fPIC
CC_GEN_FLAGS   = -fPIC

FF90_OPTFLAGS   = -fast
CC_OPTFLAGS     = -O3

#FF90_DEBUGFLAGS = -g -Mbounds -Mdclchk -DDEBUG_MODE
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
