#      ******************************************************************
#      *                                                                *
#      * File:          ADT_Common.mk                                   *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 11-26-2004                                      *
#      * Last modified: 02-23-2006                                      *
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
#      * Set the root directory for SU_MPI and the modules and object   *
#      * directories for ADT. Note that ADT_DIR must be set in the      *
#      * Makefiles that include this common file.                       *
#      *                                                                *
#      ******************************************************************

SU_MPI_DIR = $(ADT_DIR)/../SU_MPI

ADT_MODDIR = $(ADT_DIR)/mod
ADT_OBJDIR = $(ADT_DIR)/obj

#      ******************************************************************
#      *                                                                *
#      * Include the file describing the compiler settings.             *
#      *                                                                *
#      ******************************************************************

MAKE = gmake

ADT_COMPILERS = $(ADT_DIR)/config.mk
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

FF90_LOCAL_FLAGS = $(COMMAND_SEARCH_PATH_MODULES)$(SU_MPI_DIR)/mod
FF90_ALL_FLAGS   = $(FF90_LOCAL_FLAGS) $(FF90_FLAGS) \
		   $(INTEGER_PRECISION_FLAG) $(REAL_PRECISION_FLAG)

#      ******************************************************************
#      *                                                                *
#      * Rules to make objects.  Note that Cygwin requires a different  *
#      * format.                                                        *
#      *                                                                *
#      ******************************************************************

ifneq ($(origin CYGWIN), undefined)

.F90.o: Makefile
	$(FF90) $(FF90_ALL_FLAGS) -c $< -object:$(ADT_OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.F90 successfully ---"
	@echo

.f90.o: Makefile
	$(FF90) $(FF90_ALL_FLAGS) -c $< -object:$(ADT_OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.f90 successfully ---"
	@echo

else

.F90.o:	Makefile
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $(ADT_OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.F90 successfully ---"
	@echo

.f90.o:	Makefile
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $(ADT_OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.f90 successfully ---"
	@echo

endif
