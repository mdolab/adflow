#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author:        Gaetan Kenway                                   *
#      * Starting date: 03-10-2003                                      *
#      * Last modified: 28-05-2012                                      *
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

LINKER_ALL_FLAGS = -L$(SUMB_LIBDIR) -lsumb $(PV3_LINKER_FLAGS)\
		    $(CGNS_LINKER_FLAGS) $(PVM3_LINKER_FLAGS) \
		    $(PETSC_LINKER_FLAGS) $(LINKER_FLAGS)

#      ******************************************************************
#      *                                                                *
#      * Names of the Fortran and C object files to be created.         *
#      *                                                                *
#      ******************************************************************

PYTHON_OBJECTS = fortranobject.o \
		 sumbmodule.o \
		 sumb-f2pywrappers2.o 

#      ******************************************************************
#      *                                                                *
#      * Targets for make.                                              *
#      *                                                                *
#      ******************************************************************

default: all

all:

# Generate Python inlude directory
	 $(eval PYTHON_INCLUDES = $(shell $(PYTHON-CONFIG) --includes))
	 @echo "#------------------------------------------------------#"
	 @echo Python Inclue Flags $(PYTHON_INCLUDES)
	 @echo "#------------------------------------------------------#"

# Generate Numpy inlude directory
	$(eval NUMPY_INCLUDES = $(shell $(PYTHON) -c 'import numpy; print numpy.get_include()'))
	@echo "#------------------------------------------------------#"
	@echo Numpy Include Directory: $(NUMPY_INCLUDES)
	@echo "#------------------------------------------------------#"

# Generate f2py root directory
	$(eval F2PY_ROOT = $(shell $(PYTHON) get_f2py.py))
	@echo "#------------------------------------------------------#"
	@echo f2py root directory: $(F2PY_ROOT)
	@echo "#------------------------------------------------------#"

#       need to have f2py_f2cmap in this folder....
	cp f2py_f2cmap.ref .f2py_f2cmap

# Run f2py to get sumbmodule.c and sumb-f2pywrapper2.f90
	f2py sumb.pyf

# Compile c wrapper
	$(CC) $(CC_FLAGS) $(PYTHON_INCLUDES) -I$(NUMPY_INCLUDES) \
	-I$(F2PY_ROOT)/src -c sumbmodule.c

# Compile fortran object needed by all f2py modules
	$(CC) $(CC_FLAGS) $(PYTHON_INCLUDES) -I$(NUMPY_INCLUDES) \
	-c $(F2PY_ROOT)/src/fortranobject.c -o fortranobject.o

# Compiled f2py generated wrapper file
	$(FF90) $(FF90_ALL_FLAGS) -I$(SUMB_MODDIR) \
	-I$(SUMB_DIR)/src/python/fortran/suggar/ \
	-I$(SUMB_DIR)/src/python/fortran/aeroElastic -c sumb-f2pywrappers2.f90

# Final Link:
	$(FF90) -shared $(PYTHON_OBJECTS) $(LINKER_ALL_FLAGS) -o sumb.so
	python importTest.py

	mv sumb.so ../../../python
	-rm *.o
	-rm *.f90 *.f *.c
	-rm .f2py_f2cmap

clean:
	@echo "        Making clean in f2py... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)
