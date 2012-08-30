#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Andre C. Marta                                  *
#      * Starting date: 07-11-2006                                      *
#      * Last modified: 02-28-2007                                      *
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

# Running Tapenade in the reverse mode, requires linking the
# differentiated program to a separate library, written in Fortran and C,
# that defines all the PUSH... and POP... routines. These sources can be
# downloaded from the Tapenade home page:
# http://www-sop.inria.fr/tropics/tapenade.html.

ADKIT_OBJECTS = adBuffer.o \
                adStack.o  \
                dpStack.o

#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all:
	$(PMAKE) adkit

clean:
	@echo "        Making clean in adjoint/ADFirstAidKit... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)
	rm -f $(ADKIT_OBJECTS)

#      ******************************************************************
#      *                                                                *
#      * Rule to make the src part.                                     *
#      *                                                                *
#      * It is necessary to include an extra rule for the *.f files     *
#      * created by TAPENADE, since the general rulesSources.mk only    *
#      * account for f90 source files.                                  *
#      *                                                                *
#      ******************************************************************

adkit:	$(ADKIT_OBJECTS) 

.f.o: Makefile
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $(SUMB_OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.f successfully ---"
	@echo
