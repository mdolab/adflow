#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 06-12-2002                                      *
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
FF90_OBJECTS_1 = turbMod.o

FF90_OBJECTS_2 = applyAllTurbBC.o\
		 bcEddyNoWall.o\
		 bcEddyWall.o\
		 bcTurbFarfield.o\
		 bcTurbInflow.o\
		 bcTurbInterface.o\
		 bcTurbOutflow.o\
		 bcTurbSymm.o\
		 bcTurbTreatment.o\
		 bcTurbWall.o\
		 bl.o\
		 computeEddyViscosity.o\
		 curveFit.o\
		 f1SST.o\
		 initCurveFitDataKtau.o\
		 initCurveFitDataKw.o\
		 initCurveFitDataKwMod.o\
		 initCurveFitDataSa.o\
		 initCurveFitDataSae.o\
		 initCurveFitDataSST.o\
		 initCurveFitDataVf.o\
		 initKOmega.o\
		 keSolve.o\
		 kt.o\
		 ktCDterm.o\
		 ktSolve.o\
		 kw.o\
		 kwCDterm.o\
		 kwSolve.o\
		 prodKatoLaunder.o\
		 prodSmag2.o\
		 prodWmag2.o\
		 sa.o\
		 saNuFromEddyRatio.o\
		 saSolve.o\
		 SST.o\
		 SSTSolve.o\
		 tdia3.o\
		 turb2ndHalo.o\
		 turbAdvection.o\
		 turbBCNSWall.o\
		 turbResidual.o\
		 turbSolveSegregated.o\
		 unsteadyTurbSpectral.o\
		 unsteadyTurbTerm.o\
		 vf.o\
		 vfScale.o\
		 vfSolve.o\
		 writeLoglaw.o

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
	@echo "        Making clean in turbulence... "
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
