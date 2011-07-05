#! /usr/bin/env python

import sys
import pymdobuildutil

# Make sure the version of pymdobuilutil we got is new enough.
if not(int(pymdobuildutil.__version__[0]) >= 1 and
       int(pymdobuildutil.__version__[2]) >= 1):
    print "Error: You need to update pymdobuildutil to version 1.1.0 or greater"
    print "       to use this setup.py file.\n"
    sys.exit(1)

from pymdobuildutil import build

build(signature="src/python/f2py/sumb.PYF",
      f90modules=["$(MODULE_SRCDIR)",
                  "$(SU_MPI_ROOTDIR)/src",
                  "$(ADT_ROOTDIR)/src",
                  "src/python/fortran"],
      f2pydir="src/python/f2py",
      libs=["-L$(SU_MPI_LIBDIR) -L$(ADT_LIBDIR) -L$(SUMB_LIBDIR) ",
            "-linputParam -lpartitioning -linitFlow -lsolver ",
            "-lpreprocessing -lturbulence -lbcdata ",
            "-lslidingComm -lwallDistance -loutput -loverset ",
            "-lparallelIO $(PV3_LINKER_FLAGS) -lutils -lmodules -lpyFort ",
            " -lmetis -ladt -lsu_mpi $(CGNS_LINKER_FLAGS) $(PVM3_LINKER_FLAGS)"],
      builddir="python",
      makefiles=[".","src/python/fortran"],
      configfiles="SUmb_Common.mk",
      buildf2cmap=True,
      dryrun=False,
      cleanup=True)
