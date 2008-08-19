#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        C.A.(Sandy) Mader                               *
#      * Starting date: 01-15-2007                                      *
#      * Last modified: 05-08-2008                                      *
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

SUMB_DIR = ../../..
SUMB_COMMON_FILE = $(SUMB_DIR)/SUmb_Common.mk
SUMB_RULES_FILE  = $(SUMB_DIR)/rulesSources.mk
include ${SUMB_COMMON_FILE}
include ${SUMB_RULES_FILE}

#      ******************************************************************
#      *                                                                *
#      * Names of the Fortran and C object files to be created.         *
#      *                                                                *
#      ******************************************************************

vpath %.o $(SUMB_OBJDIR)

FF90_OBJECTS = adjustinflowangleforcecouplingadj_b.o\
	       applyallbcforcecouplingadj_b.o \
               bceulerwallforcecouplingadj_b.o \
               bcfarfieldforcecouplingadj_b.o \
	       bcsymmforcecouplingadj_b.o \
	       checkinputparamforcecouplingadj_b.o\
	       computeforcecouplingadj_b.o\
	       computeforcecouplingpressureadj_b.o\
	       eintarrayforcecouplingadj_b.o\
	       etotarrayforcecouplingadj_b.o\
	       extractbcstatesforcecouplingadj_b.o \
               extrapolate2ndhaloforcecouplingadj_b.o\
	       forcescouplingadj_b.o\
	       getdirvectorforcecoupling_b.o\
               getsurfacenormalscouplingadj_b.o\
	       referencestateforcecouplingadj_b.o \
	       replacebcstatesforcecouplingadj_b.o \
               setflowinfinitystateforcecouplingadj_b.o\
	       vectorrotationforcecoupling_b.o

#CC_OBJECTS   = connect_signals.o

#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all:
	$(MAKE) src

clean:
	@echo "        Making clean in adjoint... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

#      ******************************************************************
#      *                                                                *
#      * Rule to make the src part.                                     *
#      *                                                                *
#      ******************************************************************

src:	$(FF90_OBJECTS) $(CC_OBJECTS)
