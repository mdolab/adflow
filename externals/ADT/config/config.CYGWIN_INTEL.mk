#      ******************************************************************
#      *                                                                *
#      * File:          config.CYGWIN_INTEL.mk                          *
#      * Author:        Erich Elsen                                     *
#      * Starting date: 01-11-2007                                      *
#      * Last modified: 01-11-2007                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Defines the compiler settings and other commands  *
#      *              to have "make" function correctly. This file      *
#      *              defines the settings Cygwin using the intel       *
#      *              ifc/ifort compiler.                               *
#      *                                                                *
#      ******************************************************************

#      ==================================================================

#      ******************************************************************
#      *                                                                *
#      * Possibly overrule the make command to allow for parallel make. *
#      *                                                                *
#      ******************************************************************

#MAKE = make -j 2

#      ******************************************************************
#      *                                                                *
#      * F90 compiler definition; ifort is for version 8; ifc for 7.    *
#      *                                                                *
#      ******************************************************************

FF90 = ifort
#FF90 = ifc

#      ******************************************************************
#      *                                                                *
#      * Precision flags. When nothing is specified 4 byte integer and  *
#      * and 8 byte real types are used. The flag -DUSE_LONG_INT sets   *
#      * the 8 byte integer to the standard integer type. The flag      *
#      * -DUSE_SINGLE_PRECISION sets the 4 byte real to the standard    *
#      * real type.                                                     *
#      *                                                                *
#      ******************************************************************

INTEGER_PRECISION_FLAG = #-DUSE_LONG_INT
REAL_PRECISION_FLAG    = -DUSE_SINGLE_PRECISION

#      ******************************************************************
#      *                                                                *
#      * Compiler flags.                                                *
#      *                                                                *
#      ******************************************************************

COMMAND_SEARCH_PATH_MODULES = -I

FF90_GEN_FLAGS = -nologo

#FF90_OPTFLAGS   = -O2 -tpp7 -xW -unroll -ipo -ipo_obj
#FF90_OPTFLAGS   = -O2 -tpp7 -xW -unroll -ip
#FF90_OPTFLAGS   = -O2
FF90_OPTFLAGS   = -O2 -Qfpp
#FF90_OPTFLAGS    = -Od -Zi -Qfpp -CU -CB

#FF90_DEBUGFLAGS = -g -CA -CB -CS -CU -implicitnone -e90 -e95 -DDEBUG_MODE
#FF90_DEBUGFLAGS = -g -C -implicitnone -e90 -e95 -DDEBUG_MODE
#FF90_DEBUGFLAGS = -g -implicitnone -DDEBUG_MODE

FF90_FLAGS = $(FF90_GEN_FLAGS) $(FF90_OPTFLAGS) $(FF90_DEBUGFLAGS)

#      ******************************************************************
#      *                                                                *
#      * Archiver and archiver flags.                                   *
#      *                                                                *
#      ******************************************************************

AR       = ar
AR_FLAGS = -rvs
