include config.mk

# Define project name. Will define the name of the xxxx.so output file and the Python importable module
PROJECT_NAME = test_f2py

# Group all the fortran, C and compiler flags together.
FF90_ALL_FLAGS   = $(FF90_FLAGS) $(FF90_PRECISION_FLAGS) $(CUFFLAGS)
CC_ALL_FLAGS     = $(C_FLAGS) $(CC_PRECISION_FLAGS)
LINKER_ALL_FLAGS = $(LINKER_FLAGS)

# Include the full list of FORTRAN source files
include fileList

# Include the full list of directories to scan which hold the source files in the fileList
include directoryList

#The rules file includes instructions for how to compile source files into binary .o files including flags
include rules

# Set the make VPATH variable to the "dirs" variable from
# directorylist.
# dirs:=$(addprefix ./,$(dirs))
VPATH:=$(dirs)

# Now strip off the directory part. This is because fort_depend gives
# us targets that are like <filename.o>.
fortranFilesNoDir=$(notdir $(fortranFiles))
cFilesNoDir=$(notdir $(cFiles))
f77FilesNoDir=$(notdir $(f77Files))

# Generate two separate list of .F90, .f90, and f95 files using the filter command
f90Files=$(filter %.f90,$(fortranFilesNoDir))
F90Files=$(filter %.F90,$(fortranFilesNoDir))
f95Files=$(filter %.f95,$(fortranFilesNoDir))

# Finally convert all source files to .o
# This produces a list of all source files with .o endings
OFILES=$(f90Files:%.f90=%.o) $(F90Files:%.F90=%.o) $(f95Files:%.f95=%.o) $(cFilesNoDir:%.c=%.o) $(f77FilesNoDir:%.f=%.o)

# Define the file name for the auto-generated dependency file.
# The .dep file is automatically generated and contains data on the order in which the source files need to be compiled
DEP_FILE=$(PROJECT_NAME).dep

# The .pyf file tells f2py how to generate the Python interface into your Fortran source
PYF_FILE=./$(PROJECT_NAME).pyf

# ==============================================================================
#                         Recipes
# ==============================================================================
default: all

# Optionally, create the pyf file without compiling everything
pyf : $(PYF_FILE)

# The pyf signature file tells f2py how to generate a Python interface
# This command will generate a .pyf automatically
# Often, you will want to customize the pyf file to reduce the scope of your Python interface.
# In that case, this build step will leave your file as is.
# $(PYF_FILE): $(fortranFiles)
# 	@echo "Generating new pyf signature file\n"
# 	$(F2PY) $(fortranFiles) -m $(PROJECT_NAME) -h $(PYF_FILE) --overwrite-signature

# The 'source's makes all the .o files. It also depends on the dependence file "DEP_FILE".
# This command, while nondescript, does the bulk of the work in the build process
# It compiles all your Fortran and C source files into .o binary files which get linked together later
# The commands used to compile your sources are defined in the `rules` file that was included earlier
sources: $(OFILES) $(DEP_FILE)

# The .dep file depends on the source files, so it automatically gets updated
# when you change your sources
$(DEP_FILE): $(fortranFiles)
	@echo "Making dependencies file!\n"
	$(PYTHON) fort_depend.py -w -o $(DEP_FILE) -f $(fortranFiles)

# Include the automatically generated file.
include $(DEP_FILE)

# The 'lib' target builds the static archive. It depends on "sources".
lib: sources
	$(AR) $(AR_FLAGS) $(PROJECT_NAME).a $(OFILES)

# These binaries are compiled from auto-generated f2py wrapper source files and linked later
PYTHON_OBJECTS = fortranobject.o\
		 $(PROJECT_NAME)module.o\
		 $(PROJECT_NAME)-f2pywrappers2.o

# For all (driven by default)...we depend on having the library 'lib'
# up to date as well changes to the pyf file. This is where we execute
# the "PMAKE" command so everything is built in parallel.
all:
	$(PMAKE) python

python: lib $(PYF_FILE)

# ************************************************
#     SETUP INCLUDE DIRECTORIES
# ************************************************
# Generate Python inlcude directory
	 $(eval PYTHON_INCLUDES = $(shell $(PYTHON-CONFIG) --includes))
	 @echo "#------------------------------------------------------#"
	 @echo Python Include Flags $(PYTHON_INCLUDES)
	 @echo "#------------------------------------------------------#"

# Generate Numpy include directory
	$(eval NUMPY_INCLUDES = $(shell $(PYTHON) -c 'from __future__ import print_function; import numpy; print(numpy.get_include())'))
	@echo "#------------------------------------------------------#"
	@echo Numpy Include Directory: $(NUMPY_INCLUDES)
	@echo "#------------------------------------------------------#"

# Generate f2py root directory
	$(eval F2PY_ROOT = $(shell $(PYTHON) ./get_f2py.py))
	@echo "#------------------------------------------------------#"
	@echo f2py root directory: $(F2PY_ROOT)
	@echo "#------------------------------------------------------#"

# Check the consistency of the pyf file. All characters (except comments and preprocessor tags) need to be lower case!!!
	$(PYTHON) ./checkPyfForUpperCase.py

# Run the preprocessor on the pyf file:
	$(PYTHON) ./pyf_preprocessor.py real $(PYF_FILE)

# ************************************************
#     GENERATE CUSTOM PYTHON C API WRAPPER
# ************************************************
# Run f2py to get $(PROJECT_NAME)module.c and $(PROJECT_NAME)-f2pywrapper2.f90
# These automatically-generated source files constitute the Python-Fortran interface
	$(F2PY) ./$(PROJECT_NAME).pyf.autogen

# Compile the Python C interface which was automatically generated by f2py in the previous step
	$(CC) $(CC_ALL_FLAGS) $(PYTHON_INCLUDES) -I$(NUMPY_INCLUDES) \
	-I$(F2PY_ROOT)/src -c $(PROJECT_NAME)module.c

# Compile the Fortran wrapper which was automatically generated by f2py
	$(FF90) $(FF90_ALL_FLAGS) -I./ -c $(PROJECT_NAME)-f2pywrappers2.f90

# This "fortranobject" is an f2py source file which is needed by the C and Fortran wrappers
# It is the same for all f2py interfaces (not unique to this project)
	$(CC) $(CC_ALL_FLAGS) $(PYTHON_INCLUDES) -I$(NUMPY_INCLUDES) \
	-c $(F2PY_ROOT)/src/fortranobject.c -o fortranobject.o

# ************************************************
#     LINK IT ALL TOGETHER
# ************************************************
# This step links your custom Fortran code (the .a static library) with the Python API into a shared object file
# This is what Python can import and use
	$(FF90) $(CUFFLAGS) -shared $(PYTHON_OBJECTS) $(PROJECT_NAME).a $(LINKER_ALL_FLAGS) -o $(PROJECT_NAME).so

# Verify that the .so file can be imported successfully
	$(PYTHON) importTest.py -n $(PROJECT_NAME)

clean:
	rm -fr *.o *.mod
	rm -fr *.so
	rm -fr *.a
	rm -fr *.dep
	rm -fr *.c
	rm -fr *.autogen
	rm -fr *f2pywrappers*
	rm -fr *.exe
