#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Edwin van der Weide, Seonghyeon Hahn            *
#      * Starting date: 03-10-2003                                      *
#      * Last modified: 11-21-2007                                      *
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
#      * Names of the Fortran and C object files to be created.         *
#      *                                                                *
#      ******************************************************************

vpath %.o $(SUMB_OBJDIR)
FF90_OBJECTS = applyPC.o\
	       computeResidualNK.o\
	       computeResidualNK2.o\
	       destroyNKSolver.o\
	       formFunction.o\
	       formJacobian.o\
	       formJacobian_petsc.o\
	       getCurrentResidual.o\
	       getFreeStreamResidual.o\
	       NKSolver.o\
	       NKSolver_petsc.o\
	       NKUtilities.o\
	       setupNKSolver.o\
	       setupNKSolver_petsc.o\
	       snesMonitor.o\
	       LSCheck.o\
	       temp.o\

#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all:
	$(MAKE) src

clean:
	@echo "        Making clean in solver... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

#      ******************************************************************
#      *                                                                *
#      * Rule to make the src part.                                     *
#      *                                                                *
#      ******************************************************************

src:   $(FF90_OBJECTS)
