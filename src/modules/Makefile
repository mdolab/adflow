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

FF90_OBJECTS_0 = complexify.o\

FF90_OBJECTS_1 = precision.o\
		 su_cgns.o

FF90_OBJECTS_2 = bleedFlows.o\
		 constants.o\
		 coolingModelLevel0.o\
		 killSignals.o\
		 communication.o\

FF90_OBJECTS_3 = NKsolverVars.o\
		 block.o\
		 IOModule.o\
		 costFunctions.o\
		 liftDistributionData.o\

FF90_OBJECTS_4 = ADjointVars.o\
		 ADjointPETSc.o\
		 BCTypes.o\
		 blockPointers.o\
		 cgnsGrid.o\
		 cgnsNames.o\
		 commMixing.o\
		 commSliding.o\
		 couplerParam.o\
		 CpCurveFits.o\
		 extraOutput.o\
		 flowVarRefState.o\
		 inputParam.o\
		 interfaceGroups.o\
		 iteration.o\
		 monitor.o\
		 paramTurb.o\
		 section.o\
		 stencils.o\
	         diffSizes.o\
		 wallDistanceData.o\
		 myPushPopLib.o\

FILES_TO_COMPLEXIFY =  precision.F90\
		 bleedFlows.f90\
		 constants.F90\
		 coolingModelLevel0.f90\
		 liftDistributionData.F90\
		 killSignals.f90\
		 communication.f90\
		 NKsolverVars.F90\
		 block.F90\
		 IOModule.f90\
		 costFunctions.F90\
		 ADjointVars.F90\
		 ADjointPETSc.F90\
		 BCTypes.f90\
		 blockPointers.F90\
		 cgnsGrid.f90\
		 cgnsNames.f90\
		 commMixing.f90\
		 commSliding.f90\
		 couplerParam.f90\
		 CpCurveFits.f90\
		 extraOutput.f90\
		 flowVarRefState.F90\
		 inputParam.F90\
		 interfaceGroups.f90\
		 iteration.f90\
		 monitor.f90\
		 paramTurb.F90\
		 section.f90\
		 stencils.f90\
	         blockpointers_d.f90\
	         blockpointers_t.f90\
	         blockpointers_b.f90\
	         diffSizes.f90\
		 wallDistanceData.F90\
		 myPushPopLib.F90\

FILES_TO_COPY = su_cgns.F90 \
		complexify.f90 



#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all:
	$(PMAKE) src_0
	$(PMAKE) src_1
	$(PMAKE) src_2
	$(PMAKE) src_3
	$(PMAKE) src_4
	cp -f *.mod $(SUMB_MODDIR)

clean:
	@echo "        Making clean in modules... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

complexify:
	@for file in $(FILES_TO_COMPLEXIFY); do \
		if [ $$file -nt ../../src_cs/modules/$$file -o ! -e ../../src_cs/modules/$$file ]; then \
			python ../complexify.py $$file; \
			mv c_$$file ../../src_cs/modules/$$file; \
		fi;\
		done
	@for file in $(FILES_TO_COPY); do \
		if [ $$file -nt ../../src_cs/modules/$$file -o ! -e ../../src_cs/modules/$$file ]; then \
			echo copying file $$file; \
			cp $$file ../../src_cs/modules/$$file; \
		fi;\
		done

	cp Makefile ../../src_cs/modules

#      ******************************************************************
#      *                                                                *
#      * Rules to make the individual parts; these are introduced to    *
#      * have parallel make function correctly. The objects of src_1    *
#      * are independent; the ones of src_2 depend on src_1, etc.       *
#      *                                                                *
#      ******************************************************************

src_0:	$(FF90_OBJECTS_0)
src_1:	$(FF90_OBJECTS_1)
src_2:	$(FF90_OBJECTS_2)
src_3:	$(FF90_OBJECTS_3)
src_4:	$(FF90_OBJECTS_4)
