#      ******************************************************************
#      *                                                                *
#      * File:          config.APPLE_MAC_NAG.mk                         *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 01-18-2005                                      *
#      * Last modified: 02-21-2006                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Defines the compiler settings and other commands  *
#      *              to have "make" function correctly. This file      *
#      *              defines the settings for a sequential Mac         *
#      *              executable using the NAG f95 and the gcc          *
#      *              compilers.                                        *
#      *                                                                *
#      ******************************************************************

#      ==================================================================

#      ******************************************************************
#      *                                                                *
#      * Possibly overrule the make command to allow for parallel make. *
#      *                                                                *
#      ******************************************************************

#MAKE       = make -j 2

#      ******************************************************************
#      *                                                                *
#      * F90 and C compiler definitions.                                *
#      *                                                                *
#      ******************************************************************

FF90 = f95
CC   = gcc

#      ******************************************************************
#      *                                                                *
#      * Compiler flags. Note the addition of -DSEQUENTIAL_MODE.        *
#      *                                                                *
#      ******************************************************************

COMMAND_SEARCH_PATH_MODULES = -I

FF90_GEN_FLAGS = -fpp -kind=byte -mismatch_all -Wc,-fno-common \
	         -DSEQUENTIAL_MODE
CC_GEN_FLAGS   = -fno-common -DSEQUENTIAL_MODE

FF90_OPTFLAGS   = -O4
CC_OPTFLAGS     = -O3 -fexpensive-optimizations -frerun-cse-after-loop \
		  -fthread-jumps -funroll-loops -finline-functions

#FF90_DEBUGFLAGS = -g -ieee=stop -nan -C -gline -DDEBUG_MODE
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
