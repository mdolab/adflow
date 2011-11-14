#      ******************************************************************
#      *                                                                *
#      * File:          config.ASCI_QSC.mk                              *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 04-12-2004                                      *
#      * Last modified: 02-23-2006                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Defines the compiler settings and other commands  *
#      *              to have "make" function correctly. This file      *
#      *              defines the settings for a parallel executable    *
#      *              on the ASCI qsc machine in Los Alamos.            *
#      *                                                                *
#      ******************************************************************

#      ==================================================================

#      ******************************************************************
#      *                                                                *
#      * Possibly overrule the make command to allow for parallel make. *
#      *                                                                *
#      ******************************************************************

MAKE = gmake

#      ******************************************************************
#      *                                                                *
#      * Suffix of the executable such that binaries of multiple        *
#      * platforms can be stored in the bin directory.                  *
#      *                                                                *
#      ******************************************************************

EXEC_SUFFIX = _asci_qsc

#      ******************************************************************
#      *                                                                *
#      * F90 and C compiler definitions.                                *
#      *                                                                *
#      ******************************************************************

FF90 = f90
CC   = cc

#      ******************************************************************
#      *                                                                *
#      * CGNS include and linker flags.                                 *
#      *                                                                *
#      ******************************************************************

CGNS_INCLUDE_FLAGS = -I$(HOME)/CGNS2.3/include
CGNS_LINKER_FLAGS  = -L$(HOME)/CGNS2.3/lib -lcgns.ALPHA

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
#      * Compiler flags.                                                *
#      *                                                                *
#      ******************************************************************

COMMAND_SEARCH_PATH_MODULES = -I

FF90_GEN_FLAGS =
CC_GEN_FLAGS   =

FF90_OPTFLAGS  = -fast -transform_loops
CC_OPTFLAGS    = -O

DEBUGFLAGS      = -check bounds -check format -check output_conversion \
		  -check overflow -check underflow
#FF90_DEBUGFLAGS = -g $(DEBUGFLAGS)
#FF90_DEBUGFLAGS = -g3 -O3
#FF90_DEBUGFLAGS = -g
#CC_DEBUGFLAGS   = -g $(DEBUGFLAGS)

FF90_FLAGS = $(FF90_GEN_FLAGS) $(FF90_OPTFLAGS) $(FF90_DEBUGFLAGS)
CC_FLAGS   = $(CC_GEN_FLAGS)   $(CC_OPTFLAGS)   $(CC_DEBUGFLAGS)

#      ******************************************************************
#      *                                                                *
#      * pV3 and pvm3 linker flags.                                     *
#      *                                                                *
#      ******************************************************************

#PV3_FLAGS          =
#PV3_LINKER_FLAGS   =
#PVM3_LINKER_FLAGS  =
#PV3_INT_SRC_DIR    =

#      ******************************************************************
#      *                                                                *
#      * Archiver and archiver flags.                                   *
#      *                                                                *
#      ******************************************************************

AR       = ar
AR_FLAGS = -rv

#      ******************************************************************
#      *                                                                *
#      * Linker and linker flags.                                       *
#      *                                                                *
#      ******************************************************************

LINKER       = $(FF90)
LINKER_FLAGS = $(MPI_LD_FLAGS) $(FF90_OPTFLAGS) \
		-lmpi -lfmpi -lmpio -lelan -nofor_main
