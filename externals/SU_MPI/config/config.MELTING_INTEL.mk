#      ******************************************************************
#      *                                                                *
#      * File:          config.MELTING_INTEL.mk                         *
#      * Author:        Juan J. Alonso, Edwin van der Weide, M. Herrmann*
#      * Starting date: 10-03-2006                                      *
#      * Last modified: 10-03-2006                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Defines the compiler settings and other commands  *
#      *              to have "make" function correctly. This file      *
#      *              defines the settings for a parallel executable    *
#      *              for the intel compilers on the Stanford Linux     *
#                     cluster melting.                                  *
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

FF90 = /share/apps/mvapich/intel/bin/mpif90
CC   = /share/apps/mvapich/intel/bin/mpicc

#      ******************************************************************
#      *                                                                *
#      * Compiler flags. It is assumed that mpif90 is based on the      *
#      * intel compiler ifort and mpicc on the intel compiler icc.      *
#      *                                                                *
#      ******************************************************************

COMMAND_SEARCH_PATH_MODULES = -I

FF90_GEN_FLAGS = -DUSE_MPI_INCLUDE_FILE -fPIC
CC_GEN_FLAGS   =

#FF90_OPTFLAGS   = -O3 -ipo -ipo_obj
FF90_OPTFLAGS   = -O2 -tpp7 -xW -unroll -ip
CC_OPTFLAGS     = -O3

#F90_DEBUGFLAGS = -g -CA -CB -CS -CU -implicitnone -DDEBUG_MODE
#FF90_DEBUGFLAGS = -g -implicitnone -e90 -e95 -DDEBUG_MODE
#CC_DEBUGFLAGS   = -g -Wall -Wcheck -DDEBUG_MODE

FF90_FLAGS = $(FF90_GEN_FLAGS) $(FF90_OPTFLAGS) $(FF90_DEBUGFLAGS)
CC_FLAGS   = $(CC_GEN_FLAGS)   $(CC_OPTFLAGS)   $(CC_DEBUGFLAGS)

#      ******************************************************************
#      *                                                                *
#      * Archiver and archiver flags.                                   *
#      *                                                                *
#      ******************************************************************

AR       = ar
AR_FLAGS = -rvs
