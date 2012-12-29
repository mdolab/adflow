#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 02-27-2003                                      *
#      * Last modified: 02-24-2006                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Makefile to create the object files of this       *
#      * directory.                                                     *
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
FF90_OBJECTS_1 = viscSurface.o

FF90_OBJECTS_2 = allocExtraWallDistanceMemory.o\
		 computeNormalSpacing.o\
		 determineDistance.o\
		 determineDistanceQuickly.o\
		 initWallDistance.o\
		 localViscousSurfaceMesh.o\
		 updateWallDistanceAllLevels.o\
		 viscousSurfaceMesh.o\
		 wallDistance.o

FILES_TO_COMPLEXIFY = viscSurface.f90\
		 allocExtraWallDistanceMemory.f90\
		 computeNormalSpacing.f90\
		 determineDistance.f90\
		 determineDistanceQuickly.f90\
		 initWallDistance.f90\
		 localViscousSurfaceMesh.f90\
		 updateWallDistanceAllLevels.f90\
		 viscousSurfaceMesh.f90\
		 wallDistance.f90

#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all:
	$(PMAKE) src_1
	$(PMAKE) src_2

clean:
	@echo "        Making clean in wallDistance... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

complexify:
	@for file in $(FILES_TO_COMPLEXIFY); do \
		echo Checking $$file; \
		if [ $$file -nt ../../src_cs/wallDistance/$$file -o ! -e ../../src_cs/wallDistance/$$file ]; then \
			python ../complexify.py $$file; \
			mv c_$$file ../../src_cs/wallDistance/$$file; \
		fi;\
		done
	@for file in $(FILES_TO_COPY); do \
		if [ $$file -nt ../../src_cs/wallDistance/$$file -o ! -e ../../src_cs/wallDistance/$$file ]; then \
			echo copying file $$file; \
			cp $$file ../../src_cs/wallDistance/$$file; \
		fi;\
		done

	cp Makefile ../../src_cs/wallDistance


#      ******************************************************************
#      *                                                                *
#      * Rules to make the individual parts; these are introduced to    *
#      * have parallel make function correctly. The objects of src_1    *
#      * are independent; the ones of src_2 depend on src_1.            *
#      *                                                                *
#      ******************************************************************

src_1:  $(FF90_OBJECTS_1)
src_2:  $(FF90_OBJECTS_2)
