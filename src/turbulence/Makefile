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

FILES_TO_COMPLEXIFY = turbMod.F90\
		 applyAllTurbBC.F90\
		 bcEddyNoWall.f90\
		 bcEddyWall.f90\
		 bcTurbFarfield.f90\
		 bcTurbInflow.f90\
		 bcTurbInterface.f90\
		 bcTurbOutflow.f90\
		 bcTurbSymm.f90\
		 bcTurbTreatment.F90\
		 bcTurbWall.F90\
		 bl.f90\
		 computeEddyViscosity.F90\
		 curveFit.f90\
		 f1SST.f90\
		 initCurveFitDataKtau.f90\
		 initCurveFitDataKw.f90\
		 initCurveFitDataKwMod.f90\
		 initCurveFitDataSa.f90\
		 initCurveFitDataSae.f90\
		 initCurveFitDataSST.f90\
		 initCurveFitDataVf.f90\
		 initKOmega.f90\
		 keSolve.f90\
		 kt.f90\
		 ktCDterm.f90\
		 ktSolve.f90\
		 kw.f90\
		 kwCDterm.f90\
		 kwSolve.f90\
		 prodKatoLaunder.f90\
		 prodSmag2.f90\
		 prodWmag2.f90\
		 sa.F90\
		 saNuFromEddyRatio.f90\
		 saSolve.F90\
		 SST.f90\
		 SSTSolve.f90\
		 tdia3.f90\
		 turb2ndHalo.f90\
		 turbAdvection.f90\
		 turbBCNSWall.f90\
		 turbResidual.f90\
		 turbSolveSegregated.f90\
		 unsteadyTurbSpectral.f90\
		 unsteadyTurbTerm.f90\
		 vf.f90\
		 vfScale.f90\
		 vfSolve.f90\
		 writeLoglaw.f90\

FILES_TO_COPY = 

ALL_FILES = $(FILES_TO_COMPLEXIFY) $(FILES_TO_COPY)
OFILES := $(ALL_FILES:%.f90=%.o) $(ALL_FILES:%.F90=%.o) $(ALL_FILES:%.c=%.o)

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
	@echo "        Making clean in turbulence... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

complexify:
	@for file in $(FILES_TO_COMPLEXIFY); do \
		echo Checking $$file; \
		if [ $$file -nt ../../src_cs/turbulence/$$file -o ! -e ../../src_cs/turbulence/$$file ]; then \
			python ../complexify.py $$file; \
			mv c_$$file ../../src_cs/turbulence/$$file; \
		fi;\
		done
	@for file in $(FILES_TO_COPY); do \
		if [ $$file -nt ../../src_cs/turbulence/$$file -o ! -e ../../src_cs/turbulence/$$file ]; then \
			echo copying file $$file; \
			cp $$file ../../src_cs/turbulence/$$file; \
		fi;\
		done

	cp Makefile ../../src_cs/turbulence


#      ******************************************************************
#      *                                                                *
#      * Rules to make the individual parts; these are introduced to    *
#      * have parallel make function correctly. The objects of src_1    *
#      * are independent; the ones of src_2 depend on src_1.            *
#      *                                                                *
#      ******************************************************************

src_1:  $(FF90_OBJECTS_1)
src_2:  $(OFILES)
