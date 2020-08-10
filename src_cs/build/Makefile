# Include the user supplied makefil
include ../../config/config.mk

# Group all the fortran, C and compiler flags together.
FF90_ALL_FLAGS   = $(FF90_FLAGS) $(CGNS_INCLUDE_FLAGS) -I. -DUSE_COMPLEX \
		   $(PETSC_CC_INCLUDES) $(FF90_PRECISION_FLAGS) \
		   $(CGNS_VERSION_FLAG)

CC_ALL_FLAGS     = $(C_FLAGS) -I../c_defines  -I../metis-4.0  $(PETSC_CC_INCLUDES) \
		   $(CC_PRECISION_FLAGS)

LINKER_ALL_FLAGS = $(CGNS_LINKER_FLAGS) $(PETSC_LIB) $(LINKER_FLAGS)

# Include the full filelist
include ../../src/build/fileList

# Include the directory list
include ../../src/build/directoryList

# Now, systematically remove the automatically generated tapenade AD
# code.

fortranFiles:=$(patsubst adjoint/outputForward/%,,$(fortranFiles))
fortranFiles:=$(patsubst adjoint/outputReverse/%,,$(fortranFiles))
fortranFiles:=$(patsubst adjoint/outputReverseFast/%,,$(fortranFiles))

dirs:=$(patsubst adjoint/outputForward/%,,$(dirs))
dirs:=$(patsubst adjoint/outputReverse/%,,$(dirs))
dirs:=$(patsubst adjoint/outputReverseFast/%,,$(dirs))

# Include the rule list:
include ../../src/build/rules

# Set the make VPATH variable to the "dirs" variable from
# directorylist. We must first append the '../'
dirs:=$(addprefix ../,$(dirs))
VPATH:=$(dirs)

# Add the special extra complexify module to the list of files:
fortranFiles:=$(fortranFiles) modules/complexify.f90

# Need to add the ../ to each of the file lists:
fortranFiles:=$(addprefix ../,$(fortranFiles))
cFiles:=$(addprefix ../,$(cFiles))
f77Files:=$(addprefix ../,$(f77Files))

# Now strip off the directory part. This is because fort_depend gives
# us targets that are like <filename.o>.
fortranFilesNoDir=$(notdir $(fortranFiles))
cFilesNoDir=$(notdir $(cFiles))
f77FilesNoDir=$(notdir $(f77Files))

# Generate two separate list of .F90 and .f90 files using the filter command
f90Files=$(filter %.f90,$(fortranFilesNoDir))
F90Files=$(filter %.F90,$(fortranFilesNoDir))

# Define the list of real files and where the complex ones will end up.
realFiles:=$(addprefix ../../src/,$(fortranFiles))
cmplxFiles:=$(addprefix ../,$(fortranFiles))

# Finally convert all source files to .o
OFILES=$(f90Files:%.f90=%.o) $(F90Files:%.F90=%.o) $(cFilesNoDir:%.c=%.o) $(f77FilesNoDir:%.f=%.o)

# Define the file name for the auto-generated dependency file.
DEP_FILE=adflow_project.dep

# Not real targets are "phony"
.PHONY: sources lib

default: all

# The 'source's makes all the .o files. It also depends on the dependence file "DEP_FILE".
sources: $(OFILES) $(DEP_FILE)

# The .dep file depends on the source files, so it automatically gets updated
# when you change your source
$(DEP_FILE): $(fortranFiles)
	@echo "Making dependencies!\n"
	@echo "fortran Files" $(fortranFiles)
	$(PYTHON) ../../src/build/fort_depend.py -w -o $(DEP_FILE) -f $(fortranFiles)
	@echo "VPATH" $(dirs)
# Include the automatically generated file.
include $(DEP_FILE)

# The 'lib' target builds the static archive. It depends on "sources".
lib: sources
	$(AR) $(AR_FLAGS) libadflow.a $(OFILES)

# These are the extra python objects for the python library
PYTHON_OBJECTS = fortranobject.o\
		 libadflow_csmodule.o\
		 libadflow_cs-f2pywrappers2.o

# For all (driven by default)...we depend on having the library 'lib'
# up to date as well changes to the pyf file. This is where we execute
# the "PMAKE" command so everything is built in parallel.
all:
	$(PMAKE) python

python: lib ../f2py/adflow.pyf
# Generate Python inlude directory
	 $(eval PYTHON_INCLUDES = $(shell $(PYTHON-CONFIG) --includes))
	 @echo "#------------------------------------------------------#"
	 @echo Python Inclue Flags $(PYTHON_INCLUDES)
	 @echo "#------------------------------------------------------#"

# Generate Numpy inlude directory
	$(eval NUMPY_INCLUDES = $(shell $(PYTHON) -c 'import numpy; print(numpy.get_include())'))
	@echo "#------------------------------------------------------#"
	@echo Numpy Include Directory: $(NUMPY_INCLUDES)
	@echo "#------------------------------------------------------#"

# Generate f2py root directory
	$(eval F2PY_ROOT = $(shell $(PYTHON) ../f2py/get_f2py.py))
	@echo "#------------------------------------------------------#"
	@echo f2py root directory: $(F2PY_ROOT)
	@echo "#------------------------------------------------------#"

# Check the consistency of the pyf file. All carachters (except comments and preprocessor tags) need to be lower case!!!
	python ../f2py/checkPyfForUpperCase.py

# Run the preprocessor on the pyf file:
	python ../f2py/pyf_preprocessor.py complex ../f2py/adflow.pyf

# Run f2py to get adflowmodule.c and adflow-f2pywrapper2.f90
	$(F2PY) ../f2py/adflow.pyf.autogen

# Compile c wrapper. Don't use CC_ALL_FLAGS...PETSc wil F-up this command.
	$(CC) $(CC_ALL_FLAGS) $(PYTHON_INCLUDES) -I$(NUMPY_INCLUDES) \
	-I$(F2PY_ROOT)/src -c libadflow_csmodule.c

# Compile fortranobject needed by all f2py modules
	$(CC) $(CC_ALL_FLAGS) $(PYTHON_INCLUDES) -I$(NUMPY_INCLUDES) \
	-c $(F2PY_ROOT)/src/fortranobject.c -o fortranobject.o

# Compiled f2py-generated wrapper file
	$(FF90) $(FF90_ALL_FLAGS) -I./ -c libadflow_cs-f2pywrappers2.f90

# Final Link:
	$(FF90) -shared $(PYTHON_OBJECTS) libadflow.a $(LINKER_ALL_FLAGS) -o libadflow_cs.so
	$(PYTHON) importTest.py
	mv libadflow_cs.so ../../adflow
