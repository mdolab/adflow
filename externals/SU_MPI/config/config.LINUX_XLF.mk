#      ******************************************************************
#      *                                                                *
#      * File:          config.LINUX_XLF.mk                             *
#      * Author:        Patrick LeGresley                               *
#      * Starting date: 07-28-2006                                      *
#      * Last modified: 07-28-2006                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Defines the compiler settings and other commands  *
#      *              to have "make" function correctly. This file      *
#      *              defines the settings for a sequential executable  *
#      *              using the IBM xlf95 and the gcc compilers.        *
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

FF90 = xlf95
CC   = gcc

#      ******************************************************************
#      *                                                                *
#      * Compiler flags.                                                *
#      *                                                                *
#      ******************************************************************

COMMAND_SEARCH_PATH_MODULES = -I

FF90_GEN_FLAGS = -WF,-DSEQUENTIAL_MODE \
                 -qsuffix=f=f90 -qsuffix=cpp=F90 -qextname
CC_GEN_FLAGS   = -fno-common -DSEQUENTIAL_MODE

FF90_OPTFLAGS   = -O3
CC_OPTFLAGS     = -O3 -fexpensive-optimizations -frerun-cse-after-loop \
                  -fthread-jumps -funroll-loops -finline-functions

#FF90_DEBUGFLAGS = -g -C -qfullpath -qfloat=nans -WF,-DDEBUG_MODE
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
