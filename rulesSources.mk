#      ******************************************************************
#      *                                                                *
#      * File:          rulesSources.mk                                 *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 12-07-2002                                      *
#      * Last modified: 02-24-2006                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Rules to make the objects. These are the general  *
#      *              rules. If in a subdirectory different rules must  *
#      *              be used, this file should not be included.        *
#      *                                                                *
#      ******************************************************************

.F90.o:	Makefile
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $(SUMB_OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.F90 successfully ---"
	@echo

.f90.o:	Makefile
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $(SUMB_OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.f90 successfully ---"
	@echo

.c.o:   Makefile
	$(CC) $(CC_ALL_FLAGS) -c $< -o $(SUMB_OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.c successfully ---"
	@echo
