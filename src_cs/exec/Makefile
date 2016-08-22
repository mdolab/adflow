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

LINKER_ALL_FLAGS =  -L$(ADT_DIR)/lib \
		   -L$(SUMB_LIBDIR) -lsumb_cs $(CGNS_LINKER_FLAGS) \
		    $(PETSC_LINKER_FLAGS) $(LINKER_FLAGS)

#      ******************************************************************
#      *                                                                *
#      * Names of the C and Fortran object files to be created.         *
#      *                                                                *
#      ******************************************************************

CC_OBJECTS   = SUmb.o
FF90_OBJECTS = SUmb_fortran.o

FILES_TO_COMPLEXIFY = 

FILES_TO_COPY = SUmb.c\
		SUmb_fortran.f90\


OBJECTS      = $(CC_OBJECTS) $(FF90_OBJECTS)

#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all:
	$(PMAKE) src
	make SUmb

clean:
	@echo "        Making clean in exec... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

complexify:
	@for file in $(FILES_TO_COMPLEXIFY); do \
		echo Checking $$file; \
		if [ $$file -nt ../../src_cs/exec/$$file -o ! -e ../../src_cs/exec/$$file ]; then \
			python ../complexify.py $$file; \
			mv c_$$file ../../src_cs/exec/$$file; \
		fi;\
		done
	@for file in $(FILES_TO_COPY); do \
		if [ $$file -nt ../../src_cs/exec/$$file -o ! -e ../../src_cs/exec/$$file ]; then \
			echo copying file $$file; \
			cp $$file ../../src_cs/exec/$$file; \
		fi;\
		done

	cp Makefile ../../src_cs/exec


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
