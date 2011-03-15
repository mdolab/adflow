#      ******************************************************************
#      *                                                                *
#      * File:          config.LINUX_INTEL.mk                           *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 01-18-2005                                      *
#      * Last modified: 02-21-2006                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Defines the compiler settings and other commands  *
#      *              to have "make" function correctly. This file      *
#      *              defines the settings for a sequential Linux       *
#      *              executable using the intel and the gcc compilers. *
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

FF90 = ifort
#FF90 = ifc
CC   = gcc

#      ******************************************************************
#      *                                                                *
#      * Compiler flags. Note the addition of -DSEQUENTIAL_MODE.        *
#      *                                                                *
#      ******************************************************************

COMMAND_SEARCH_PATH_MODULES = -I

FF90_GEN_FLAGS = -DSEQUENTIAL_MODE
CC_GEN_FLAGS   = -DSEQUENTIAL_MODE

#FF90_OPTFLAGS   = -O2 -tpp7 -xW -unroll -ipo -ipo_obj
#FF90_OPTFLAGS   = -O2 -tpp7 -xW -unroll -ip
#FF90_OPTFLAGS   = -O2
FF90_OPTFLAGS   = -O2 -tpp7 -axW -ip
#CC_OPTFLAGS     = -O3 -fexpensive-optimizations -frerun-cse-after-loop \
		  -fthread-jumps -funroll-loops -finline-functions

#FF90_DEBUGFLAGS = -g -CA -CB -CS -CU -implicitnone -e90 -e95 -DDEBUG_MODE
#FF90_DEBUGFLAGS = -g -C -implicitnone -e90 -e95 -DDEBUG_MODE
#FF90_DEBUGFLAGS = -g -implicitnone -DDEBUG_MODE
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
