#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Edwin van der Weide, Seonghyeon Hahn            *
#      * Starting date: 09-26-2004                                      *
#      * Last modified: 09-13-2007                                      *
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
FF90_OBJECTS_1 = BCDataMod.o

FF90_OBJECTS_2 = allocMemBCData.o\
		 BCDataDomainInterface.o\
		 BCDataIsothermalWall.o\
		 BCDataMassBleedInflow.o\
		 BCDataMassBleedOutflow.o\
		 BCDataSubsonicInflow.o\
		 BCDataSubsonicOutflow.o\
		 BCDataSupersonicInflow.o\
		 cart1D_InterpolSubface.o\
		 computeHtot.o\
		 extractFromDataSet.o\
		 initBCData.o\
		 initBCDataDomainInterfaces.o\
		 nondimBoundData.o\
		 radialInterpolSubface.o\
		 setBCDataCoarseGrid.o\
		 setBCDataFineGrid.o\
		 setBCTurb.o\
		 setInletFreestreamTurb.o\
		 setSupersonicInletFreeStream.o\
		 storeFreestreamSubface.o\
		 unitVectorsCylSystem.o

#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all:
	$(MAKE) src_1
	$(MAKE) src_2

clean:
	@echo "        Making clean in bcdata... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

#      ******************************************************************
#      *                                                                *
#      * Rules to make the individual parts; these are introduced to    *
#      * have parallel make function correctly. The objects of src_1    *
#      * are independent; the ones of src_2 depend on src_1.            *
#      *                                                                *
#      ******************************************************************

src_1:  $(FF90_OBJECTS_1)
src_2:  $(FF90_OBJECTS_2)
