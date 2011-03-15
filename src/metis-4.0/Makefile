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
CC_OBJECTS = balance.o\
	     bucketsort.o\
	     ccgraph.o\
	     coarsen.o\
	     compress.o\
	     debug.o\
	     estmem.o\
	     fm.o\
	     fortran.o\
	     frename.o\
	     graph.o\
	     initpart.o\
	     kmetis.o\
	     kvmetis.o\
	     kwayfm.o\
	     kwayrefine.o\
	     kwayvolfm.o\
	     kwayvolrefine.o\
	     match.o\
	     mbalance.o\
	     mbalance2.o\
	     mcoarsen.o\
	     memory.o\
	     mesh.o\
	     meshpart.o\
	     mfm.o\
	     mfm2.o\
	     mincover.o\
	     minitpart.o\
	     minitpart2.o\
	     mkmetis.o\
	     mkwayfmh.o\
	     mkwayrefine.o\
	     mmatch.o\
	     mmd.o\
	     mpmetis.o\
	     mrefine.o\
	     mrefine2.o\
	     mutil.o\
	     myqsort.o\
	     ometis.o\
	     parmetis.o\
	     pmetis.o\
	     pqueue.o\
	     refine.o\
	     separator.o\
	     sfm.o\
	     srefine.o\
	     stat.o\
	     subdomains.o\
	     timing.o\
	     util.o

#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all:
	$(MAKE) src

clean:
	@echo "        Making clean in metis... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

#      ******************************************************************
#      *                                                                *
#      * Rule to make the src part.                                     *
#      *                                                                *
#      ******************************************************************

src:	$(CC_OBJECTS)
