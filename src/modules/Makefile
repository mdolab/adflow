#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Edwin van der Weide, Seonghyeon Hahn            *
#      * Starting date: 12-10-2002                                      *
#      * Last modified: 03-02-2006                                      *
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
FF90_OBJECTS_1 = precision.o\
		 su_cgns.o

FF90_OBJECTS_2 = bleedFlows.o\
		 constants.o\
		 coolingModelLevel0.o\
		 killSignals.o

FF90_OBJECTS_3 = block.o\
		 IOModule.o

FF90_OBJECTS_4 = BCTypes.o\
		 blockPointers.o\
		 cgnsGrid.o\
		 cgnsNames.o\
		 commMixing.o\
		 commSliding.o\
		 communication.o\
		 couplerParam.o\
		 CpCurveFits.o\
		 extraOutput.o\
		 flowVarRefState.o\
		 inputParam.o\
		 interfaceGroups.o\
		 iteration.o\
		 monitor.o\
		 paramTurb.o\
		 section.o

#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all:
	$(MAKE) src_1
	$(MAKE) src_2
	$(MAKE) src_3
	$(MAKE) src_4
	cp -f *.mod $(SUMB_MODDIR)

clean:
	@echo "        Making clean in modules... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

#      ******************************************************************
#      *                                                                *
#      * Rules to make the individual parts; these are introduced to    *
#      * have parallel make function correctly. The objects of src_1    *
#      * are independent; the ones of src_2 depend on src_1, etc.       *
#      *                                                                *
#      ******************************************************************

src_1:	$(FF90_OBJECTS_1)
src_2:	$(FF90_OBJECTS_2)
src_3:	$(FF90_OBJECTS_3)
src_4:	$(FF90_OBJECTS_4)
