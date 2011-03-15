#      ******************************************************************
#      *                                                                *
#      * File:          config.REDSTORM.mk                              *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 10-12-2007                                      *
#      * Last modified: 10-12-2007                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Defines the compiler settings and other commands  *
#      *              to have "make" function correctly. This file      *
#      *              defines the settings for an executable on         *
#      *              RedStorm.                                         *
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

FF90 = ftn
CC   = cc

#      ******************************************************************
#      *                                                                *
#      * Compiler flags.                                                *
#      *                                                                *
#      ******************************************************************

COMMAND_SEARCH_PATH_MODULES = -I

FF90_GEN_FLAGS =
CC_GEN_FLAGS   =

FF90_OPTFLAGS   = -fastsse -Mipa=fast
CC_OPTFLAGS     = -fast

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
