#      ******************************************************************
#      *                                                                *
#      * File:          SUmb_Common.mk                                  *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 12-07-2002                                      *
#      * Last modified: 03-02-2006                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Defines the root directory and includes the files *
#      *              describing the architecture and the compiler      *
#      *              settings.                                         *
#      *                                                                *
#      ******************************************************************

#      ==================================================================

#      ******************************************************************
#      *                                                                *
#      * Set the root directory for SU_MPI and ADT and the modules,     *
#      * object, lib and bin directories for SUmb. Note that SUMB_DIR   *
#      * must be set in the Makefiles that include this common file.    *
#      *                                                                *
#      ******************************************************************

SU_MPI_DIR = $(SUMB_DIR)/externals/SU_MPI
ADT_DIR    = $(SUMB_DIR)/externals/ADT

SUMB_MODDIR = $(SUMB_DIR)/mod
SUMB_OBJDIR = $(SUMB_DIR)/obj
SUMB_BINDIR = $(SUMB_DIR)/bin
SUMB_LIBDIR = $(SUMB_DIR)/lib

#      ******************************************************************
#      *                                                                *
#      * Include the file describing the compiler settings.             *
#      *                                                                *
#      ******************************************************************

MAKE = gmake

SUMB_COMPILERS = $(SUMB_DIR)/config.mk
ifneq ($(MAKECMDGOALS),clean)
include ${SUMB_COMPILERS}
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
#      * The macro's ADDITIONAL_FF90_FLAGS and ADDITIONAL_CC_FLAGS      *
#      * make it possible that every subdirectory adds its own specific *
#      * compiler flags, if necessary.                                  *
#      *                                                                *
#      ******************************************************************

FF90_LOCAL_FLAGS = $(COMMAND_SEARCH_PATH_MODULES)$(SUMB_MODDIR) \
		   $(COMMAND_SEARCH_PATH_MODULES)$(SU_MPI_DIR)/mod \
		   $(COMMAND_SEARCH_PATH_MODULES)$(ADT_DIR)/mod
FF90_ALL_FLAGS   = $(FF90_LOCAL_FLAGS) $(CGNS_INCLUDE_FLAGS) \
		   $(FF90_FLAGS) $(ADDITIONAL_FF90_FLAGS) $(PV3_FLAGS) \
		   $(FF90_PRECISION_FLAGS)
CC_LOCAL_FLAGS   = -I$(SUMB_DIR)/src/c_defines -I.
CC_ALL_FLAGS     = $(CC_LOCAL_FLAGS) $(CC_FLAGS) $(ADDITIONAL_CC_FLAGS) \
		   $(CC_PRECISION_FLAGS)
