#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 02-12-2003                                      *
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
#      * Names of the C object files to be created.                     *
#      *                                                                *
#      ******************************************************************

vpath %.o $(SUMB_OBJDIR)

FILES_TO_COMPLEXIFY = 

FILES_TO_COPY = balance.c\
	     bucketsort.c\
	     ccgraph.c\
	     coarsen.c\
	     compress.c\
	     debug.c\
	     estmem.c\
	     fm.c\
	     fortran.c\
	     frename.c\
	     graph.c\
	     initpart.c\
	     kmetis.c\
	     kvmetis.c\
	     kwayfm.c\
	     kwayrefine.c\
	     kwayvolfm.c\
	     kwayvolrefine.c\
	     match.c\
	     mbalance.c\
	     mbalance2.c\
	     mcoarsen.c\
	     memory.c\
	     mesh.c\
	     meshpart.c\
	     mfm.c\
	     mfm2.c\
	     mincover.c\
	     minitpart.c\
	     minitpart2.c\
	     mkmetis.c\
	     mkwayfmh.c\
	     mkwayrefine.c\
	     mmatch.c\
	     mmd.c\
	     mpmetis.c\
	     mrefine.c\
	     mrefine2.c\
	     mutil.c\
	     myqsort.c\
	     ometis.c\
	     parmetis.c\
	     pmetis.c\
	     pqueue.c\
	     refine.c\
	     separator.c\
	     sfm.c\
	     srefine.c\
	     stat.c\
	     subdomains.c\
	     timing.c\
	     util.c\
	     metis.h\
	     defs.h\
	     struct.h\
	     macros.h\
	     rename.h\
	     proto.h\

ALL_FILES = $(FILES_TO_COMPLEXIFY) $(FILES_TO_COPY)
OFILES := $(ALL_FILES:%.f90=%.o) $(ALL_FILES:%.F90=%.o) $(ALL_FILES:%.c=%.o)

#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all:
	$(PMAKE) src

clean:
	@echo "        Making clean in metis... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

complexify:
	@for file in $(FILES_TO_COMPLEXIFY); do \
		echo Checking $$file; \
		if [ $$file -nt ../../src_cs/metis-4.0/$$file -o ! -e ../../src_cs/metis-4.0/$$file ]; then \
			python complexify.py $$file; \
			mv c_$$file ../../src_cs/metis-4.0/$$file; \
		fi;\
		done
	@for file in $(FILES_TO_COPY); do \
		if [ $$file -nt ../../src_cs/metis-4.0/$$file -o ! -e ../../src_cs/metis-4.0/$$file ]; then \
			echo copying file $$file; \
			cp $$file ../../src_cs/metis-4.0/$$file; \
		fi;\
		done
	cp Makefile ../../src_cs/metis-4.0

#      ******************************************************************
#      *                                                                *
#      * Rule to make the src part.                                     *
#      *                                                                *
#      ******************************************************************

src:	$(OFILES)
