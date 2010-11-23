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

FF90_OBJECTS = 	adjustinflowangleforcesadj_b.o\
		bceulerwallforcesadj_b.o\
		computeforceandmomentadj_b.o\
		computeforcesadj_b.o\
		computepressureforcesadj_b.o\
		forcesadj_b.o\
		getdirvectorforces_b.o\
		getsurfacenormalsforcesadj_b.o\
		vectorrotationforces_b.o\


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
