#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 02-24-2006                                      *
#      * Last modified: 03-02-2006                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Driver Makefile for the routines needed for       *
#      *              multi-disciplinary applications.                  *
#      *                                                                *
#      ******************************************************************

#      ==================================================================

#      ******************************************************************
#      *                                                                *
#      * Include the settings for the directory containing the routines *
#      * for running SUmb for multi-disciplinary applications.          *
#      *                                                                *
#      ******************************************************************

SUMB_MD_DIR = .
SUMB_MD_COMMON_FILE = $(SUMB_MD_DIR)/SUmb_MD_Common.mk
include ${SUMB_MD_COMMON_FILE}

#      ******************************************************************
#      *                                                                *
#      * The subdirectories where the sources are located.              *
#      *                                                                *
#      ******************************************************************

SUBDIR_SRC_MD = aeroElastic\
		couplerAPI\
		hqp\
		suggar

#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all:
	@for subdir in $(SUBDIR_SRC_MD) ; \
		do \
			echo "making $@ in $$subdir"; \
			echo; \
			(cd $$subdir && gmake) || exit 1; \
		done
	gmake libpyfort

clean:
	@echo "        Making clean in python/fortran... "
	rm -f *~ obj/*
	@for subdir in $(SUBDIR_SRC_MD) ; \
		do \
			echo; \
			echo "making $@ in $$subdir"; \
			echo; \
			(cd $$subdir && gmake $@) || exit 1; \
		done

#      ******************************************************************
#      *                                                                *
#      * Rule to make the library.                                      *
#      *                                                                *
#      ******************************************************************

libpyfort:
	@echo "        creating library ..."
	@echo
	$(AR) $(AR_FLAGS) $(SUMB_LIBDIR)/libsumb_pyfort.a \
	$(SUMB_MD_OBJDIR)/*.o
	@echo
	@echo "        library libsumb_pyfort.a created."
	@echo
