# adjointRules: A special set of rules for the adjoint:

.F90.o:	Makefile
	$(FF90) $(FF90_ADJOINT_FLAGS) -c $< -o $(SUMB_OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.F90 successfully ---"
	@echo

.f90.o:	Makefile
	$(FF90) $(FF90_ADJOINT_FLAGS) -c $< -o $(SUMB_OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.f90 successfully ---"
	@echo

.c.o:   Makefile
	$(CC) $(CC_ALL_FLAGS) -c $< -o $(SUMB_OBJDIR)/$(@F)
	@echo
	@echo "        --- Compiled $*.c successfully ---"
	@echo
