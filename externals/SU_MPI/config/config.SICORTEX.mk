#      ******************************************************************
#      *                                                                *
#      * File:          config.SICORTEX.mk                              *
#      * Author:        Ruben E. Perez                                  *
#      * Starting date: 07-31-2009                                      *
#      * Last modified: 07-31-2009                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Defines the compiler settings and other commands  *
#      *              to have "make" function correctly. This file      *
#      *              defines the settings for a sequential executable  *
#      *              on the SICORTEX (MIPS64 processors) type machines *
#      *              using the pathscale compilers.                    *
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
#      * F90 and C compiler definitions.                                *
#      *                                                                *
#      ******************************************************************

FF90 = pathf95
CC   = pathcc

#      ******************************************************************
#      *                                                                *
#      * Compiler flags. Note the addition of -DSEQUENTIAL_MODE.        *
#      *                                                                *
#      ******************************************************************

COMMAND_SEARCH_PATH_MODULES = -I

FF90_GEN_FLAGS = -fPIC -DSEQUENTIAL_MODE
CC_GEN_FLAGS   = -fPIC -DSEQUENTIAL_MODE

#FF90_OPTFLAGS   = -O3 -IPA -OPT:Ofast
#CC_OPTFLAGS     = -O3 -IPA -OPT:Ofast
FF90_OPTFLAGS   = -O3 -OPT:Ofast
CC_OPTFLAGS     = -O3 -OPT:Ofast

#FF90_DEBUGFLAGS = -g -Wall -DDEBUG_MODE
#CC_DEBUGFLAGS   = -g -Wall -DDEBUG_MODE

FF90_FLAGS = $(FF90_GEN_FLAGS) $(FF90_OPTFLAGS) $(FF90_DEBUGFLAGS)
CC_FLAGS   = $(CC_GEN_FLAGS)   $(CC_OPTFLAGS)   $(CC_DEBUGFLAGS)

#      ******************************************************************
#      *                                                                *
#      * Archiver and archiver flags.                                   *
#      *                                                                *
#      ******************************************************************

AR       = ar
AR_FLAGS = -rvs
