#      ******************************************************************
#      *                                                                *
#      * File:          config.LINUX_INTEL_64_MPICH.mk                  *
#      * Author:        Edwin van der Weide, Andre C. Marta             *
#      * Starting date: 12-07-2002                                      *
#      * Last modified: 02-05-2007                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Defines the compiler settings and other commands  *
#      *              to have "make" function correctly. This file      *
#      *              defines the settings for a parallel Linux         *
#      *              executable in combination with MPICH. Assumed is  *
#      *              that mpif90 and mpicc are based on the intel and  *
#      *              gcc compilers respectively.                       *
#      *              Tuned for the ADL 64-bit workstations.            *
#      *                                                                *
#      ******************************************************************

#      ==================================================================

#      ******************************************************************
#      *                                                                *
#      * Possibly overrule the make command to allow for parallel make. *
#      *                                                                *
#      ******************************************************************

MAKE = make -j 8

#      ******************************************************************
#      *                                                                *
#      * Suffix of the executable such that binaries of multiple        *
#      * platforms can be stored in the bin directory.                  *
#      *                                                                *
#      ******************************************************************

EXEC_SUFFIX = _altix_mpich2

#      ******************************************************************
#      *                                                                *
#      * F90 and C compiler definitions.                                *
#      *                                                                *
#      ******************************************************************

FF90 = /opt/mpich2-icc10/bin/mpif90
CC   = /opt/mpich2-icc10/bin/mpicc

#      ******************************************************************
#      *                                                                *
#      * CGNS include and linker flags.                                 *
#      *                                                                *
#      ******************************************************************

#CGNS_INCLUDE_FLAGS = -DUSE_NO_CGNS
#CGNS_INCLUDE_FLAGS = -I/usr/local/include
#CGNS_LINKER_FLAGS  = -L/usr/local/lib64 -lcgns.intel
CGNS_INCLUDE_FLAGS = -I$(HOME)/SRC/cgnslib_2.4
CGNS_LINKER_FLAGS  = -L$(HOME)/SRC/cgnslib_2.4/LINUX -lcgns


#      ******************************************************************
#      *                                                                *
#      * Precision flags. When nothing is specified 4 byte integer and  *
#      * and 8 byte real types are used. The flag -DUSE_LONG_INT sets   *
#      * the 8 byte integer to the standard integer type. The flag      *
#      * -DUSE_SINGLE_PRECISION sets the 4 byte real to the standard    *
#      * real type and -DUSE_QUADRUPLE_PRECISION will make 16 byte      *
#      * reals the default type. The latter option may not work on all  *
#      * platforms.                                                     *
#      *                                                                *
#      ******************************************************************

#FF90_INTEGER_PRECISION_FLAG = -DUSE_LONG_INT
#CC_INTEGER_PRECISION_FLAG   = -DUSE_LONG_INT

#FF90_REAL_PRECISION_FLAG = -DUSE_SINGLE_PRECISION
#FF90_REAL_PRECISION_FLAG = -DUSE_QUADRUPLE_PRECISION
#CC_REAL_PRECISION_FLAG   = -DUSE_SINGLE_PRECISION
#CC_REAL_PRECISION_FLAG   = -DUSE_QUADRUPLE_PRECISION

FF90_PRECISION_FLAGS = $(FF90_INTEGER_PRECISION_FLAG) \
		       $(FF90_REAL_PRECISION_FLAG)
CC_PRECISION_FLAGS   = $(CC_INTEGER_PRECISION_FLAG) \
		       $(CC_REAL_PRECISION_FLAG)

#      ******************************************************************
#      *                                                                *
#      * Compiler flags. It is assumed that mpif90 is based on the      *
#      * intel compiler ifort and mpicc on the gcc gnu compiler.        *
#      *                                                                *
#      ******************************************************************

# Intel® Fortran Compiler for Linux version 9.x
#
# Type 'ifort -help' at the command prompt for information about the
# available compiler options.
#
# Performance
# -----------
# -O2        Enable optimizations (DEFAULT)
# -axW       Generate code specialized for Intel Pentium 4 and compatible 
#            Intel processors.
# -ip        Enable single-file IP optimizations (within files)
# -unroll[n] Set maximum number of times to unroll loops.  Omit n to use
#            default heuristics.  Use n=0 to disable loop unroller.
#
# Debug
# -----------
# -g               Produce symbolic debug information in object file
#                  (implies -O0 when another optimization option is 
#                   not explicitly set)
# -debug [keyword] Enable debug information and control output of 
#                  enhanced debug information.
#                  keywords:all, full, minimal, none (same as -nodebug),
#                           inline-debug-info, variable-locations,
#                           semantic-stepping, extended
# -ftrapuv         Trap uninitialized variables
# -print-multi-lib Print information about libraries being used
#
# Language
# --------
# -implicitnone    Set IMPLICIT NONE by default
# -CB              Runtime checks for out-of-bounds array 
#                  subscript/substring refs (same as -check bounds)
# -check <keyword> Check run-time conditions
#                  keywords: all, none (same as -nocheck),
#                   [no]arg_temp_created, [no]bounds,[no]format, [no]overflow,
#                   [no]output_conversion, [no]power, [no]args
#
# Compiler Diagnostics
# --------------------
# -w90, -w95 suppress messages about use of non-standard Fortran
# -warn <keyword>  Specifies the level of warning messages issued.
#                   keywords: all, none (same as -nowarn),
#                     [no]alignments, [no]declarations, [no]errors,
#                     [no]general, [no]ignore_loc, [no]interfaces,
#                     [no]stderrors, [no]truncated_source, [no]uncalled,
#                     [no]uninitialized, [no]unused, [no]usage
#
# Miscellaneous
# -------------
#
#-f[no-]pic, -f[no-]PIC
#              generate position independent code (OFF by default)

COMMAND_SEARCH_PATH_MODULES = -I

FF90_GEN_FLAGS = -DHAS_ISNAN
CC_GEN_FLAGS   =

FF90_OPTFLAGS   = -O2 -axW -unroll -ip -fPIC
CC_OPTFLAGS     = -O -fPIC

#FF90_DEBUGFLAGS = -g -CB -implicitnone -DDEBUG_MODE
#CC_DEBUGFLAGS   = -g -Wall -pedantic -DDEBUG_MODE

FF90_FLAGS = $(FF90_GEN_FLAGS) $(FF90_OPTFLAGS) $(FF90_DEBUGFLAGS)
CC_FLAGS   = $(CC_GEN_FLAGS)   $(CC_OPTFLAGS)   $(CC_DEBUGFLAGS)

#      ******************************************************************
#      *                                                                *
#      * pV3 and pvm3 linker flags.                                     *
#      *                                                                *
#      ******************************************************************

#PV3_FLAGS          = -DUSE_PV3
#PV3_LINKER_FLAGS   = -L/usr/local/pV3/clients/LINUX-INTEL -lpV3
#PVM3_LINKER_FLAGS  = -L/usr/local/pvm3/lib/LINUX -lgpvm3 -lpvm3
#PV3_INT_SRC_DIR    = src/pv3Interface

#      ******************************************************************
#      *                                                                *
#      * PETSc include and linker flags. The flag -DUSE_NO_PETSC        *
#      * indicates that the PETSc libraries should not be included in   *
#      * the built process. These are only needed when running ADjoint  *
#      * problems, that use PETSc to solve the system of equations.     *
#      * If -DUSE_NO_PETSC is uncommented then the PETSc include and    *
#      * linker flags should be commented out.                          *
#      *                                                                *
#      ******************************************************************
#      ******************************************************************
#      *                                                                *
#      * PETSc include and linker flags. The flag -DUSE_NO_PETSC        *
#      * indicates that the PETSc libraries should not be included in   *
#      * the built process. These are only needed when running ADjoint  *
#      * problems, that use PETSc to solve the system of equations.     *
#      * If -DUSE_NO_PETSC is uncommented then the PETSc include and    *
#      * linker flags should be commented out.                          *
#      *                                                                *
#      ******************************************************************


X11_DIR    = /usr/X11R6/lib
#PETSC_DIR  = /opt/petsc-2.3.1-p16
#PETSC_ARCH = linux-gnu-c-real-debug
#PETSC_DIR  = /opt/petsc-2.3.1-p16-opt
PETSC_DIR  = /opt/petsc-2.3.3-p7-mpich2-icc10
PETSC_ARCH = linux-gnu-c-opt#linux-gnu-c-real-opt
#PETSC_INCLUDE_FLAGS = -DUSE_NO_PETSC
PETSC_INCLUDE_FLAGS = -I$(PETSC_DIR) -I$(PETSC_DIR)/bmake/$(PETSC_ARCH) -I$(PETSC_DIR)/include -I$(PETSC_DIR)/include/mpiuni
PETSC_LINKER_FLAGS  = -L$(PETSC_DIR)/lib/$(PETSC_ARCH) -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc -L$(X11_DIR) -lX11 -L$(PETSC_DIR)/externalpackages/fblaslapack/$(PETSC_ARCH) -lflapack -lfblas -lm

#      ******************************************************************
#      *                                                                *
#      * Archiver and archiver flags.                                   *
#      *                                                                *
#      ******************************************************************

AR       = ar
AR_FLAGS = -rvs

#      ******************************************************************
#      *                                                                *
#      * Linker and linker flags.                                       *
#      *                                                                *
#      ******************************************************************

LINKER       = $(FF90)
LINKER_FLAGS = $(FF90_OPTFLAGS) -nofor_main -llapack -lg2c
