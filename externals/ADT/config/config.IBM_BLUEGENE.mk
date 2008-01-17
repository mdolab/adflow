#      ******************************************************************
#      *                                                                *
#      * File:          config.IBM_BLUEGENE.mk                          *
#      * Author:        Edwin van der Weide                             *
#      * Starting date: 07-07-2005                                      *
#      * Last modified: 02-23-2006                                      *
#      *                                                                *
#      ******************************************************************

#      ******************************************************************
#      *                                                                *
#      * Description: Defines the compiler settings and other commands  *
#      *              to have "make" function correctly. This file      *
#      *              defines the settings for the IBM xlf90 and cc     *
#      *              compiler for the BlueGene machine of the          *
#      *              San Diego Supercomputer Center.                   *
#      *                                                                *
#      ******************************************************************

#      ==================================================================

#      ******************************************************************
#      *                                                                *
#      * Possibly overrule the make command to allow for parallel make. *
#      *                                                                *
#      ******************************************************************

MAKE = make

#      ******************************************************************
#      *                                                                *
#      * F90 definition.                                                *
#      *                                                                *
#      ******************************************************************

FF90 = mpxlf90

#      ******************************************************************
#      *                                                                *
#      * Precision flags. When nothing is specified 4 byte integer and  *
#      * and 8 byte real types are used. The flag -DUSE_LONG_INT sets   *
#      * the 8 byte integer to the standard integer type. The flag      *
#      * -DUSE_SINGLE_PRECISION sets the 4 byte real to the standard    *
#      * real type.                                                     *
#      *                                                                *
#      ******************************************************************

#INTEGER_PRECISION_FLAG = -WF,-DUSE_LONG_INT
#REAL_PRECISION_FLAG    = -WF,-DUSE_SINGLE_PRECISION
#REAL_PRECISION_FLAG    = -WF,-DUSE_QUADRUPOLE_PRECISION

#      ******************************************************************
#      *                                                                *
#      * Compiler flags.                                                *
#      *                                                                *
#      ******************************************************************

COMMAND_SEARCH_PATH_MODULES = -I

FF90_GEN_FLAGS = -qsuffix=f=f90 -qsuffix=cpp=F90 -qarch=440

FF90_OPTFLAGS  = -O3

#FF90_DEBUGFLAGS = -g -C -qfullpath -qfloat=nans -qsigtrap -WF,-DDEBUG_MODE

FF90_FLAGS = $(FF90_GEN_FLAGS) $(FF90_OPTFLAGS) $(FF90_DEBUGFLAGS)

#      ******************************************************************
#      *                                                                *
#      * Archiver and archiver flags.                                   *
#      *                                                                *
#      ******************************************************************

AR       = ar
AR_FLAGS = -rvs
