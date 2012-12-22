#      ******************************************************************
#      *                                                                *
#      * File:          ADT_Common.mk                                   *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 11-26-2004                                      *
#      * Last modified: 02-23-2006                                      *
#      *                                                                *
#      ******************************************************************

ADT_MODDIR = $(ADT_DIR)/mod
ADT_OBJDIR = $(ADT_DIR)/obj

# We use the sumb config file which already has the necessary
# information
ADT_COMPILERS = $(ADT_DIR)/../../config.mk

ifneq ($(MAKECMDGOALS),clean)
include ${ADT_COMPILERS}
endif


#      ******************************************************************
#      *                                                                *
#      * Redefine .SUFFIXES to be sure all the desired ones are         *
#      * included.                                                      *
#      *                                                                *
#      ******************************************************************

.SUFFIXES: .o .f .F .f90 .F90 .c .C .cc .cpp .h .hh .H

#      ******************************************************************
#      *                                                                *
#      * Arguments of make clean.                                       *
#      *                                                                *
#      ******************************************************************

MAKE_CLEAN_ARGUMENTS = *~ *.o *.mod *.il *.stb

#      ******************************************************************
#      *                                                                *
#      * Compiler flags to compile the sources.                         *
#      *                                                                *
#      ******************************************************************

FF90_ALL_FLAGS   = $(FF90_FLAGS) $(FF90_PRECISION_FLAGS)
