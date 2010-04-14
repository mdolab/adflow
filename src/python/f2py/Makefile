#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Edwin van der Weide, Seonghyeon Hahn            *
#      * Starting date: 03-10-2003                                      *
#      * Last modified: 11-21-2007                                      *
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

SUMB_DIR = ../../..
SUMB_COMMON_FILE = $(SUMB_DIR)/SUmb_Common.mk
SUMB_RULES_FILE  = $(SUMB_DIR)/rulesSources.mk
include ${SUMB_COMMON_FILE}
include ${SUMB_RULES_FILE}

#      ******************************************************************
#      *                                                                *
#      * Linker flags.                                                  *
#      *                                                                *
#      ******************************************************************

LINKER_ALL_FLAGS = -L$(SU_MPI_DIR)/lib -L$(ADT_DIR)/lib \
		   -L$(SUMB_LIBDIR) -lsumb -lsumb_pyfort $(PV3_LINKER_FLAGS) -ladt \
		   -lsu_mpi $(CGNS_LINKER_FLAGS) $(PVM3_LINKER_FLAGS) \
		    $(PETSC_LINKER_FLAGS) $(LINKER_FLAGS)

#      ******************************************************************
#      *                                                                *
#      * Names of the Fortran and C object files to be created.         *
#      *                                                                *
#      ******************************************************************

SUFFIX = _parallel

PYTHON_OBJECTS = fortranobject.o \
		 sumb$(SUFFIX)module.o \
		 sumb$(SUFFIX)-f2pywrappers.o \
		 sumb$(SUFFIX)-f2pywrappers2.o 

#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all:
#	cp sumb.PYF sumb.pyf
	python createpyfFile.py sumb.PYF $(SUFFIX)
#       need to have f2py_f2cmap in this folder....
	cp f2py_f2cmap.ref .f2py_f2cmap
	f2py sumb$(SUFFIX).pyf
	$(CC) $(CC_ALL_FLAGS) -I$(NUMPY_ROOT_DIR)/numpy/core/include -I$(NUMPY_ROOT_DIR)/numpy/f2py/src -I$(PYTHON_INCLUDE_DIR) -c sumb$(SUFFIX)module.c
	$(CC) $(CC_ALL_FLAGS) -I$(NUMPY_ROOT_DIR)/numpy/core/include -I$(NUMPY_ROOT_DIR)/numpy/f2py/src -I$(PYTHON_INCLUDE_DIR) -c $(NUMPY_ROOT_DIR)/numpy/f2py/src/fortranobject.c -o fortranobject.o
	$(FF90) $(FF90_ALL_FLAGS) -I$(SUMB_MODDIR) -c sumb$(SUFFIX)-f2pywrappers.f
	$(FF90) $(FF90_ALL_FLAGS) -I$(SUMB_MODDIR) -I$(SUMB_DIR)/src/python/fortran/suggar/ -I$(SUMB_DIR)/src/python/fortran/aeroElastic -c sumb$(SUFFIX)-f2pywrappers2.f90
	$(FF90) -shared $(PYTHON_OBJECTS) $(LINKER_ALL_FLAGS) -o sumb$(SUFFIX).so
	python importTest.py $(SUFFIX)

	mv sumb$(SUFFIX).so ../../../python
	-rm *.pyf
	-rm *.o
	-rm *.f90 *.f *.c
	-rm .f2py_f2cmap

clean:
	@echo "        Making clean in f2py... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

#      ******************************************************************
#      *                                                                *
#      * Rule to make the src part.                                     *
#      *                                                                *
#      ******************************************************************

src:   $(FF90_OBJECTS) $(CC_OBJECTS)
