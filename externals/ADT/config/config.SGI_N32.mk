#      ******************************************************************
#      *                                                                *
#      * File:          config_N32.SGI.mk                               *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 01-18-2005                                      *
#      * Last modified: 02-23-2006                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Defines the compiler settings and other commands  *
#      *              to have "make" function correctly. This file      *
#      *              defines the settings for a 32 bit SGI version.    *
#      *                                                                *
#      ******************************************************************

#      ==================================================================

#      ******************************************************************
#      *                                                                *
#      * Possibly overrule the make command to allow for parallel make. *
#      *                                                                *
#      ******************************************************************

MAKE = gmake -j 8

#      ******************************************************************
#      *                                                                *
#      * F90 definition.                                                *
#      *                                                                *
#      ******************************************************************

FF90 = f90

#      ******************************************************************
#      *                                                                *
#      * Precision flags. When nothing is specified 4 byte integer and  *
#      * and 8 byte real types are used. The flag -DUSE_LONG_INT sets   *
#      * the 8 byte integer to the standard integer type. The flag      *
#      * -DUSE_SINGLE_PRECISION sets the 4 byte real to the standard    *
#      * real type.                                                     *
#      *                                                                *
#      ******************************************************************

INTEGER_PRECISION_FLAG =# -DUSE_LONG_INT
REAL_PRECISION_FLAG    =# -DUSE_SINGLE_PRECISION

#      ******************************************************************
#      *                                                                *
#      * Compiler flags.                                                *
#      * The Ofast option is based on the 600 MHZ IP35 processor.       *
#      *                                                                *
#      ******************************************************************

COMMAND_SEARCH_PATH_MODULES = -I

FF90_32_64_FLAGS = -mips4 -n32

FF90_GEN_FLAGS = $(FF90_32_64_FLAGS)

#FF90_OPTFLAGS   = -r14000 -Ofast=ip35 -IPA

FF90_OPTFLAGS   = -O3

DEBUGFLAGS      = -DEBUG:conform_check=ON:div_check=3:subscript_check=ON:trap_uninitialized=ON:varargs_interface_check=ON:verbose_runtime=ON
#FF90_DEBUGFLAGS = -g -fullwarn -check_bounds $(DEBUGFLAGS) -DDEBUG_MODE
#FF90_DEBUGFLAGS = -g

FF90_FLAGS = $(FF90_GEN_FLAGS) $(FF90_OPTFLAGS) $(FF90_DEBUGFLAGS)

#      ******************************************************************
#      *                                                                *
#      * Archiver and archiver flags.                                   *
#      *                                                                *
#      ******************************************************************

AR       = ar
AR_FLAGS = -rv
