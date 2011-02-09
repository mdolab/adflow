#      ******************************************************************
#      *                                                                *
#      * File:          SUmb_MB_Common.mk                               *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 03-02-2006                                      *
#      * Last modified: 03-02-2006                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: General include file for the routines to make     *
#      *              SUmb working in a multi-disciplinary environment. *
#      *                                                                *
#      ******************************************************************

#      ==================================================================

#      ******************************************************************
#      *                                                                *
#      * Include the settings for the entire code.                      *
#      *                                                                *
#      ******************************************************************

SUMB_DIR = $(SUMB_MD_DIR)/../../..
SUMB_COMMON_FILE = $(SUMB_DIR)/SUmb_Common.mk
include ${SUMB_COMMON_FILE}

#      ******************************************************************
#      *                                                                *
#      * Set the object directory for the routines of these directories.*
#      *                                                                *
#      ******************************************************************

SUMB_MD_OBJDIR = $(SUMB_MD_DIR)/obj

#      ******************************************************************
#      *                                                                *
#      * Define the rules to make objects. Do not use the rules defined *
#      * in rulesSources.mk, because the objects files from this        *
#      * directory should be stored in a different directory.           *
#      *                                                                *
#      ******************************************************************

.F90.o:	Makefile
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $(SUMB_MD_OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.F90 successfully ---"
	@echo

.f90.o:	Makefile
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $(SUMB_MD_OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.f90 successfully ---"
	@echo

.c.o:   Makefile
	$(CC) $(CC_ALL_FLAGS) -c $< -o $(SUMB_MD_OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.c successfully ---"
	@echo
