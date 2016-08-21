.. _sumb_install:

Install
--------

All the core computations in ``SUmb`` are coded in Fortran.  It
is therefore necessary to build this library before using ``SUmb``.

To see a list of architectures that ``SUmb`` has been known to
compile on run::
   
   make

from the root directory. 

The easiest approach to to try the closest one to your system and
attempt a build using (for example)::

   make LINUX_INTEL_OPENMPI

SUmb has been successfully compiled on LINUX, and OS X with either
ifort or gfortran.

If everything was successful, the following lines will be printed to
the screen (near the end)::

   Testing if module sumb can be imported...
   Module sumb was successfully imported.

If you don't see this, it will be necessary to configure the build
manually. To configure manually, first copy a default configuration
file from the defaults folder like this (run this in the root
directory)::
  
   cp config/defaults/config.LINUX_INTEL_OPENMPI.mk config

Now open ``config/config.LINUX_INTEL_OPENMPI.mk`` which should look like::

  # ----------------------------------------------------------------------
  # Config file for Intel ifort  with OpenMPI
  # ----------------------------------------------------------------------
  
  # ------- Define a possible parallel make (use PMAKE = make otherwise)--
  PMAKE = make -j 4
  
  # ------- Define the MPI Compilers--------------------------------------
  FF90 = mpif90
  CC   = mpicc
  
  #-------- Define Exec Suffix for executable ----------------------------
  EXEC_SUFFIX = _linux_intel_openmpi
  
  # ------- Define Precision Flags ---------------------------------------
  # Options for Integer precision flags: -DUSE_LONG_INT
  # Options for Real precision flags: -DUSE_SINGLE_PRECISION, -DUSE_QUADRUPLE_PRECISION
  # Default (nothing specified) is 4 byte integers and 8 byte reals
  
  FF90_INTEGER_PRECISION_FLAG =
  FF90_REAL_PRECISION_FLAG    = 
  CC_INTEGER_PRECISION_FLAG   =
  CC_REAL_PRECISION_FLAG      = 

  # ------- Define CGNS Inlcude and linker flags -------------------------
  CGNS_INCLUDE_FLAGS = -I/usr/local/include
  CGNS_LINKER_FLAGS = -Wl,-rpath,/usr/local/lib -lcgns
  
  # ------- Define Compiler Flags ----------------------------------------
  
  FF90_GEN_FLAGS = -DHAS_ISNAN 
  CC_GEN_FLAGS   = -DHAS_ISNAN  
  
  FF90_OPT_FLAGS   =  -fPIC -r8 -g -O2
  CC_OPT_FLAGS     = -O -fPIC
  
  FF90_DEBUG_FLAGS = #-check bounds -check all
  CC_DEBUG_FLAGS   = #-g -Wall -pedantic -DDEBUG_MODE
  
  # ------- Define Archiver  and Flags -----------------------------------
  AR       = ar
  AR_FLAGS = -rvs

  # ------- Define Linker Flags ------------------------------------------
  LINKER       = $(FF90)
  LINKER_FLAGS = -nofor_main
  
  # ------- Define Petsc Info --- Should not need to modify this -----
  include ${PETSC_DIR}/conf/variables
  PETSC_INCLUDE_FLAGS=${PETSC_CC_INCLUDES} -I$(PETSC_DIR)
  PETSC_LINKER_FLAGS=${PETSC_LIB}
  
  # Combine flags from above -- don't modify here
  FF90_FLAGS = $(FF90_GEN_FLAGS) $(FF90_OPT_FLAGS) $(FF90_DEBUG_FLAGS)
  CC_FLAGS   = $(CC_GEN_FLAGS)   $(CC_OPT_FLAGS)   $(CC_DEBUG_FLAGS)
  FF90_PRECISION_FLAGS = $(FF90_INTEGER_PRECISION_FLAG)$(FF90_REAL_PRECISION_FLAG)
  CC_PRECISION_FLAGS   = $(CC_INTEGER_PRECISION_FLAG) $(CC_REAL_PRECISION_FLAG)
  
  # Define potentially different python, python-config and f2py executables:
  PYTHON = python
  PYTHON-CONFIG = python-config
  F2PY = f2py

It is most likely that you need to modify the ``CGNS_INCLUDE_FLAGS``
and the ``CGNS_LINKER_FLAGS`` variables. It is also necessary to
``PETSc`` already compiled including support for
``SuperLU_dist``. After changes to the configuration file, run ``make
clean`` before attempting a new build. 

