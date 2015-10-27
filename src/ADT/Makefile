#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 11-26-2004                                      *
#      * Last modified: 04-26-2006                                      *
#      *                                                                *
#      ******************************************************************
#      ******************************************************************
#      *                                                                *
#      * Include the common file.                                       *
#      *                                                                *
#      ******************************************************************

SUMB_DIR = ../..
SUMB_COMMON_FILE = $(SUMB_DIR)/SUmb_Common.mk
SUMB_RULES_FILE  = $(SUMB_DIR)/rulesSources.mk
include ${SUMB_COMMON_FILE}
include ${SUMB_RULES_FILE}

#      ******************************************************************
#      *                                                                *
#      * Names of the Fortran object files to be created.               *
#      * These are split in sets, such that no dependencies are present *
#      * within a set. This guarantees that parallel make works         *
#      * correctly.                                                     *
#      *                                                                *
#      ******************************************************************

vpath %.o $(SUMB_OBJDIR)

FF90_OBJECTS_2 = gradD2Hexa.o\
		 hessD2Hexa.o\
		 minD2Hexa.o\
		 newtonStep.o
FF90_OBJECTS_3 = adtUtils.o
FF90_OBJECTS_4 = adtBuild.o\
		 adtLocalSearch.o
FF90_OBJECTS_5 = adtSearch.o
FF90_OBJECTS_6 = adtAPI.o

FILES_TO_COMPLEXIFY =	gradD2Hexa.f90\
			hessD2Hexa.f90\
			minD2Hexa.f90\
			newtonStep.f90\
			adtUtils.F90\
			adtBuild.f90\
			adtLocalSearch.f90\
			adtSearch.f90\
			adtAPI.f90

FILES_TO_COPY = 	adtData.F90\



#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all: 
	$(PMAKE) src_2
	$(PMAKE) src_3
	$(PMAKE) src_4
	$(PMAKE) src_5
	$(PMAKE) src_6
	cp -f *.mod $(SUMB_MODDIR)

clean:
	@echo "        Making clean in ADT ... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)


complexify:
	@for file in $(FILES_TO_COMPLEXIFY); do \
		echo Checking $$file; \
		if [ $$file -nt ../../src_cs/ADT/$$file -o ! -e ../../src_cs/ADT/$$file ]; then \
			python ../complexify.py $$file; \
			mv c_$$file ../../src_cs/ADT/$$file; \
		fi;\
		done
	@for file in $(FILES_TO_COPY); do \
		if [ $$file -nt ../../src_cs/ADT/$$file -o ! -e ../../src_cs/ADT/$$file ]; then \
			echo copying file $$file; \
			cp $$file ../../src_cs/ADT/$$file; \
		fi;\
		done
	cp Makefile ../../src_cs/ADT/

#      ******************************************************************
#      *                                                                *
#      * Rules to make the individual parts; these are introduced to    *
#      * have parallel make function correctly. The objects of src_1    *
#      * are independent; the ones of src_2 depend on src_1, etc.       *
#      *                                                                *
#      ******************************************************************


src_2:	$(FF90_OBJECTS_2)
src_3:	$(FF90_OBJECTS_3)
src_4:	$(FF90_OBJECTS_4)
src_5:	$(FF90_OBJECTS_5)
src_6:	$(FF90_OBJECTS_6)
