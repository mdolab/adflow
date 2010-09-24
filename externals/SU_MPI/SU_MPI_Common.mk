#      ******************************************************************
#      *                                                                *
#      * File:          SU_MPI_Common.mk                                *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 11-26-2004                                      *
#      * Last modified: 02-21-2006                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Defines the root directory, includes the files    *
#      *              describing the architecture and the compiler      *
#      *              settings.                                         *
#      *                                                                *
#      ******************************************************************

#      ==================================================================

#      ******************************************************************
#      *                                                                *
#      * Set the modules and object directories for SU_MPI. Note that   *
#      * SU_MPI_DIR must be set in the Makefiles that include this      *
#      * common file.                                                   *
#      *                                                                *
#      ******************************************************************

SU_MPI_MODDIR = $(SU_MPI_DIR)/mod
SU_MPI_OBJDIR = $(SU_MPI_DIR)/obj

#      ******************************************************************
#      *                                                                *
#      * Include the file describing the compiler settings.             *
#      *                                                                *
#      ******************************************************************

MAKE = gmake

SU_MPI_COMPILERS = $(SU_MPI_DIR)/config.mk
ifneq ($(MAKECMDGOALS),clean)
include ${SU_MPI_COMPILERS}
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

FF90_LOCAL_FLAGS = $(COMMAND_SEARCH_PATH_MODULES)$(SU_MPI_MODDIR)
FF90_ALL_FLAGS   = $(FF90_LOCAL_FLAGS) $(FF90_FLAGS)
CC_LOCAL_FLAGS   =
CC_ALL_FLAGS     = $(CC_LOCAL_FLAGS) $(CC_FLAGS)

#      ******************************************************************
#      *                                                                *
#      * Rules to make objects.                                         *
#      *                                                                *
#      ******************************************************************

.F90.o:	Makefile
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $(SU_MPI_OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.F90 successfully ---"
	@echo

.f90.o:	Makefile
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $(SU_MPI_OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.f90 successfully ---"
	@echo

.c.o:   Makefile
	$(CC) $(CC_ALL_FLAGS) -c $< -o $(SU_MPI_OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.c successfully ---"
	@echo
