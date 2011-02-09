#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 12-10-2002                                      *
#      * Last modified: 10-27-2005                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Makefile to create the executable.                *
#      *                                                                *
#      ******************************************************************

#      ==================================================================

#      ******************************************************************
#      *                                                                *
#      * Include the settings for the entire code.                      *
#      *                                                                *
#      ******************************************************************

SUMB_DIR = ../..
SUMB_COMMON_FILE = $(SUMB_DIR)/SUmb_Common.mk
include ${SUMB_COMMON_FILE}

#      ******************************************************************
#      *                                                                *
#      * Define the rules to make objects. Do not use the rules defined *
#      * in rulesSources.mk, because the objects files from this        *
#      * directory should not be stored in -lsumb.                      *
#      *                                                                *
#      ******************************************************************

.F90.o: Makefile
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $*.o
	@echo
	@echo "        --- Compiled $*.F90 successfully ---"
	@echo

.f90.o: Makefile $(MODULE_SRC)
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $*.o
	@echo
	@echo "        --- Compiled $*.f90 successfully ---"
	@echo

.c.o:   Makefile
	$(CC) $(CC_ALL_FLAGS) -c $< -o $*.o
	@echo
	@echo "        --- Compiled $*.c successfully ---"
	@echo

#      ******************************************************************
#      *                                                                *
#      * Linker flags.                                                  *
#      *                                                                *
#      ******************************************************************

LINKER_ALL_FLAGS = -L$(SU_MPI_DIR)/lib -L$(ADT_DIR)/lib \
		   -L$(SUMB_LIBDIR) -lsumb $(PV3_LINKER_FLAGS) -ladt \
		   -lsu_mpi $(CGNS_LINKER_FLAGS) $(PVM3_LINKER_FLAGS) \
		   $(LINKER_FLAGS)

#      ******************************************************************
#      *                                                                *
#      * Names of the C and Fortran object files to be created.         *
#      *                                                                *
#      ******************************************************************

CC_OBJECTS   = SUmb.o
FF90_OBJECTS = SUmb_fortran.o

OBJECTS      = $(CC_OBJECTS) $(FF90_OBJECTS)

#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all:
	$(MAKE) src
	gmake SUmb

clean:
	@echo "        Making clean in exec... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

#      ******************************************************************
#      *                                                                *
#      * Rule to make the src part. This could be incorporated directly *
#      * in "make SUmb", but is kept separately for consistency.        *
#      *                                                                *
#      ******************************************************************

src:	$(OBJECTS)

#      ******************************************************************
#      *                                                                *
#      * Rule to make the executable.                                   *
#      *                                                                *
#      ******************************************************************

SUmb:	$(OBJECTS)
	@echo "     creating SUmb$(EXEC_SUFFIX) ..."
	@echo
	$(LINKER) $(OBJECTS) $(LINKER_ALL_FLAGS) \
	-o $(SUMB_BINDIR)/SUmb$(EXEC_SUFFIX)
	@echo
	@echo "     SUmb$(EXEC_SUFFIX) created."
	@echo
