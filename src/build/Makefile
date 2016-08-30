# Include the user supplied makefil
include ../../config/config.mk

# Group all the fortran, C and compiler flags together. 
FF90_ALL_FLAGS   = $(FF90_FLAGS) $(CGNS_INCLUDE_FLAGS) -I. \
		   $(PETSC_CC_INCLUDES) $(FF90_PRECISION_FLAGS)

CC_ALL_FLAGS     = $(C_FLAGS) -I../c_defines  -I../metis-4.0  $(PETSC_CC_INCLUDES) \
		   $(CC_PRECISION_FLAGS)

LINKER_ALL_FLAGS = $(CGNS_LINKER_FLAGS) $(PETSC_LIB) $(LINKER_FLAGS)

# Set the make VPATH variable. This tells make where to find source 
VPATH:=../adjoint/:../bcdata/:../inputParam/:../NKSolver/:../adjoint/:../partitioning/:../solver/:../utils/:../ADT/:../metis-4.0/:../output/:../wallDistance/:../initFlow/:../modules/:../preprocessing/:../slidingComm/:../turbulence/:../warping/:../adjoint/residualInput/:../adjoint/outputForward/:../adjoint/outputReverse/:../adjoint/outputReverseFast/:../metis-4.0/:../adjoint/ADFirstAidKit/:../overset


# These are the rules to make each of the 4 different types of files. 

%.o : %.f90
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $*.o
	@echo
	@echo "        --- Compiled $*.f90 successfully ---"
	@echo

%.o : %.F90
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $*.o
	@echo
	@echo "        --- Compiled $*.F90 successfully ---"
	@echo

%.o : %.c
	$(CC) $(CC_ALL_FLAGS) -c $< -o $*.o
	@echo
	@echo "        --- Compiled $*.c successfully ---"
	@echo


%.o: %.f
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $*.o
	@echo
	@echo "        --- Compiled $*.f successfully ---"
	@echo

# List out all the fortran files: (grouped by directory)

fortranFiles=	ADT/adtUtils.F90\
		ADT/adtBuild.F90\
		ADT/adtLocalSearch.F90\
		ADT/adtSearch.F90\
		ADT/adtAPI.F90\
		\
		adjoint/outputForward/utils_d.f90\
		adjoint/outputForward/flowutils_d.f90\
		adjoint/outputForward/bcroutines_d.f90\
		adjoint/outputForward/turbutils_d.f90\
		adjoint/outputForward/solverutils_d.f90\
		adjoint/outputForward/aleutils_d.f90\
		adjoint/outputForward/flowutils_d.f90\
		adjoint/outputForward/bcroutines_d.f90\
		adjoint/outputForward/turbutils_d.f90\
		adjoint/outputForward/solverutils_d.f90\
		adjoint/outputForward/fluxes_d.f90\
		adjoint/outputForward/turbbcroutines_d.f90\
		adjoint/outputForward/surfaceintegrations_d.f90\
		adjoint/outputForward/residuals_d.f90\
		adjoint/outputForward/sa_d.f90\
		adjoint/outputForward/walldistance_d.f90\
		adjoint/outputForward/initializeflow_d.f90\
		adjoint/outputForward/xhalo_block_d.f90\
		adjoint/outputForward/metric_block_d.f90\
		adjoint/outputForward/getcostfunction2_d.f90\
		adjoint/outputForward/boundarynormals_d.f90\
		adjoint/outputForward/volume_block_d.f90\
		adjoint/outputForward/block_res_d.f90  \
		\
		adjoint/outputReverse/getcostfunction2_b.f90\
		adjoint/outputReverse/metric_block_b.f90\
		adjoint/outputReverse/xhalo_block_b.f90\
		adjoint/outputReverse/block_res_b.f90  \
		adjoint/outputReverse/getcostfunction2_b.f90\
		adjoint/outputReverse/boundarynormals_b.f90\
		adjoint/outputReverse/volume_block_b.f90\
		adjoint/outputReverseFast/BCExtra.F90\
		adjoint/outputReverse/utils_b.f90\
		adjoint/outputReverse/flowutils_b.f90\
		adjoint/outputReverse/bcroutines_b.f90\
		adjoint/outputReverse/turbutils_b.f90\
		adjoint/outputReverse/solverutils_b.f90\
		adjoint/outputReverse/aleutils_b.f90\
		adjoint/outputReverse/flowutils_b.f90\
		adjoint/outputReverse/bcroutines_b.f90\
		adjoint/outputReverse/turbutils_b.f90\
		adjoint/outputReverse/solverutils_b.f90\
		adjoint/outputReverse/turbbcroutines_b.f90\
		adjoint/outputReverse/turbbcroutines_b.f90\
		adjoint/outputReverse/residuals_b.f90\
		adjoint/outputReverse/surfaceintegrations_b.f90\
		adjoint/outputReverse/sa_b.f90\
		adjoint/outputReverse/fluxes_b.f90\
		adjoint/outputReverse/walldistance_b.f90\
		adjoint/outputReverse/initializeflow_b.f90\
		\
		adjoint/outputReverseFast/utils_fast_b.f90\
		adjoint/outputReverseFast/flowutils_fast_b.f90\
		adjoint/outputReverseFast/bcroutines_fast_b.f90\
		adjoint/outputReverseFast/turbutils_fast_b.f90\
		adjoint/outputReverseFast/solverutils_fast_b.f90\
		adjoint/outputReverseFast/aleutils_fast_b.f90\
		adjoint/outputReverseFast/flowutils_fast_b.f90\
		adjoint/outputReverseFast/bcroutines_fast_b.f90\
		adjoint/outputReverseFast/turbutils_fast_b.f90\
		adjoint/outputReverseFast/solverutils_fast_b.f90\
		adjoint/outputReverseFast/fluxes_fast_b.f90\
		adjoint/outputReverseFast/turbbcroutines_fast_b.f90\
		adjoint/outputReverseFast/sa_fast_b.f90\
		adjoint/outputReverseFast/BCExtra.F90\
		\
		adjoint/residualInput/block_res.F90\
		adjoint/residualInput/metric_block.F90\
		adjoint/residualInput/xhalo_block.F90\
		adjoint/residualInput/costFunction2.F90\
		\
		adjoint/adjointPythonInterfaceFunctions.F90\
		adjoint/computeAeroCoef.f90\
		adjoint/createPETScVars.F90\
		adjoint/destroyPETScVars.F90\
		adjoint/preallocation.F90\
		adjoint/getSolution.f90\
		adjoint/initializePETSc.F90\
		adjoint/finalizePETSc.F90\
		adjoint/preprocessingADjoint.F90\
		adjoint/setupPETScKsp.F90\
		adjoint/updateFlow.f90\
		adjoint/updateReferencePoint.f90\
		adjoint/updateRotationRate.f90\
		adjoint/matrixFreeRoutines.F90\
		adjoint/alloc_derivative_values.F90\
		adjoint/dealloc_derivative_values.F90\
		adjoint/allocPCMem.F90\
		adjoint/determineStencil.F90\
		adjoint/setupAllResidualMatricesFwd.F90\
		adjoint/setupColoring.F90\
		adjoint/setPointers_d.f90\
		adjoint/setDiffSizes.F90\
		adjoint/zeroADSeeds.F90\
		adjoint/setFDReference.F90\
		adjoint/referenceShockSensor.F90\
		adjoint/PCMatMult.F90\
		adjoint/solveADjointTransposePETSc.F90\
		adjoint/saveAdjoint.F90\
		adjoint/setupStateResidualMatrix.F90\
		adjoint/setupPCMatrix.F90\
		\
		bcdata/BCData.F90\
		\
		initFlow/variableReading.F90\
		initFlow/initializeFlow.F90\
		\
		inputParam/inputParamRoutines.F90\
		\
		modules/complexify.f90\
		modules/precision.F90\
		modules/su_cgns.F90\
		modules/constants.F90\
		modules/killSignals.f90\
		modules/communication.f90\
		modules/block.F90\
		modules/blockPointers.F90\
		modules/kd_tree.f90\
		modules/IOModule.f90\
		modules/costFunctions.F90\
		modules/surfaceFamilies.F90\
		modules/adtData.F90\
		modules/BCPointers.F90\
		modules/BCPointers_b.F90\
		modules/BCPointers_d.F90\
		modules/BCPointers_fast_b.F90\
		modules/wallDistanceData.F90\
		modules/ADjointVars.F90\
		modules/ADjointPETSc.F90\
		modules/cgnsGrid.f90\
		modules/cgnsNames.f90\
		modules/commMixing.f90\
		modules/commSliding.f90\
		modules/couplerParam.f90\
		modules/CpCurveFits.f90\
		modules/extraOutput.f90\
		modules/flowVarRefState.F90\
		modules/inputParam.F90\
		modules/interfaceGroups.f90\
		modules/iteration.f90\
		modules/monitor.f90\
		modules/paramTurb.F90\
		modules/section.f90\
		modules/stencils.f90\
		modules/diffSizes.f90\
		modules/overset.F90\
		\
		NKSolver/NKSolvers.F90\
		\
		output/outputMod.F90\
		output/writeCGNSSurface.F90\
		output/tecplotIO.F90\
		output/writeCGNSGrid.F90\
		output/writeCGNSVolume.F90\
		output/writeSol.F90\
		\
		overset/stringOps.F90\
		overset/cartMesh.F90\
		overset/computeInterpolationParallel.F90\
		overset/buildClusterWalls.F90\
		overset/createZipperMesh.F90\
		overset/bowtieAndIsolationElimination.F90\
		overset/determineClusters.F90 \
		overset/determineClusterAreas.F90 \
		overset/determineClusterMarchDist.F90\
		overset/checkOverset.F90 \
		overset/computeCellWallPoint.F90\
		overset/exchangeSurfaceDelta.F90\
		overset/floodInteriorCells.F90\
		overset/oversetLoadBalance.F90\
		overset/qsortFringeType.F90\
		overset/qsortEdgeType.F90\
		overset/qsortPocketEdgeType.F90\
		overset/binSearchPocketEdgeType.F90\
		overset/binSearchNodes.F90\
		overset/initializeOBlock.F90\
		overset/initializeOWall.F90\
		overset/initializeFringes.F90\
		overset/initializeOFringes.F90\
		overset/determineDonors.F90\
		overset/packingRoutines.F90\
		overset/fringeSearch.F90\
		overset/oversetUtilities.F90\
		overset/irregularCellCorrection.F90\
		overset/fringeReduction.F90\
		overset/finalOversetCommStructures.F90\
		overset/setIBlankArray.F90\
		overset/exchangeFringes.f90\
		overset/surfaceCorrection.F90\
		overset/emptyOversetComm.F90\
		overset/wallSearch.F90\
		overset/initBCDataiBlank.F90\
		overset/exchangeSurfaceIBlanks.F90\
		overset/determineZipperWallAssociation.F90\
		overset/makeBoundaryStrings.F90\
		overset/writeDualGrid.F90\
		overset/exchangeStatus.F90\
		overset/exchangeStatusTranspose.F90\
		overset/surfaceDeviation.F90\
		\
		partitioning/partitionMod.F90\
		partitioning/loadBalance.F90\
		partitioning/gridChecking.F90\
		partitioning/readCGNSGrid.F90\
		partitioning/partitioning.F90\
		\
		preprocessing/bcHalo.f90\
		preprocessing/coarse1to1Subface.f90\
		preprocessing/haloList.F90\
		preprocessing/indirectHalo.f90\
		preprocessing/periodicInfo.f90\
		preprocessing/coarse1to1Subface.f90\
		preprocessing/allocMemHaloList.f90\
		preprocessing/cellRangeSubface.f90\
		preprocessing/checkCoarse1to1.f90\
		preprocessing/checkSymmetry.f90\
		preprocessing/closestDirectHalos.f90\
		preprocessing/coarseDonorInfo.f90\
		preprocessing/createCoarseBlocks.f90\
		preprocessing/determineCommPattern.f90\
		preprocessing/determineFaceHalos.f90\
		preprocessing/determineIndirectHalos.f90\
		preprocessing/determineNumberOfHalos.f90\
		preprocessing/determineNcellGlobal.f90\
		preprocessing/determinePeriodicFaces.f90\
		preprocessing/exchangeCoor.f90\
		preprocessing/exchangeCoor_d.f90\
		preprocessing/exchangeCoor_b.f90\
		preprocessing/faceRotationMatrices.f90\
		preprocessing/finalCommStructures.f90\
		preprocessing/indirectHalosPerLevel.f90\
		preprocessing/init2ndLevelCellHalos.f90\
		preprocessing/setSurfaceFamilyInfo.F90\
		preprocessing/mdUpdateRoutines.f90\
		preprocessing/metric.f90\
		preprocessing/preprocessing.f90\
		preprocessing/qsortHaloListType.f90\
		preprocessing/qsortIndHaloType.f90\
		preprocessing/qsortPeriodicSubfacesHaloType.f90\
		preprocessing/setFamilyInfoFaces.f90\
		preprocessing/setGlobalCellsAndNodes.f90\
		preprocessing/setPorosities.f90\
		preprocessing/shiftCoorAndVolumes.f90\
		preprocessing/unitVectorsInAxialPlane.f90\
		preprocessing/update1to1Coarse.f90\
		preprocessing/updateCoorFineMesh.f90\
		preprocessing/viscSubfaceInfo.f90\
		preprocessing/xhalo.f90\
		preprocessing/setReferenceVolume.F90\
		preprocessing/bsearchCGNSPeriodicType.f90\
		preprocessing/determinePeriodicData.F90\
		\
		slidingComm/cumulativeNSendReceives.f90\
		slidingComm/cylCoorNodesOnBlockFace.f90\
		slidingComm/determineNSlices.f90\
		slidingComm/determineRotationMatrices.f90\
		slidingComm/determineSections.f90\
		slidingComm/haloType.f90\
		slidingComm/initMemMixingPlane.f90\
		slidingComm/initMemSliding.f90\
		slidingComm/initThisSlide.f90\
		slidingComm/interpolateSlide.f90\
		slidingComm/localNodesAndQuads.f90\
		slidingComm/localSurfaceGrids.f90\
		slidingComm/mixingDonorInterpol.f90\
		slidingComm/mixingHaloInterpol.f90\
		slidingComm/mixingIntervals.f90\
		slidingComm/mixingPlaneComm.f90\
		slidingComm/mixingPlaneInterpol.f90\
		slidingComm/mySubfacesSlide.f90\
		slidingComm/nLocalNodesAndQuads.f90\
		slidingComm/procSliding.f90\
		slidingComm/qsortMixingIntervalType.f90\
		slidingComm/qsortSortedDonorCommType.f90\
		slidingComm/releaseMemSliding.f90\
		slidingComm/slidingComm.f90\
		slidingComm/slidingIntersections.f90\
		slidingComm/slidingMesh.f90\
		slidingComm/slidingMeshGroups.f90\
		slidingComm/updateExternalCommSlide.f90\
		slidingComm/updateInterpolSendBuf.f90\
		slidingComm/updateLocalCommSlide.f90\
		slidingComm/updateRotationInfo.f90\
		slidingComm/updateSlidingAllLevels.f90\
		slidingComm/updateSlidingCommPattern.f90\
		slidingComm/volumeInterpol.f90\
		slidingComm/localSubfaces.f90\
		slidingComm/mixingData.F90\
		slidingComm/thisSlide.f90\
		slidingComm/tmpSliding.f90\
		slidingComm/updateComm.F90\
		slidingComm/bsearchSortedDonorCommType.f90\
		\
		solver/BCRoutines.F90\
		solver/solverUtils.F90\
		solver/fluxes.F90\
		solver/surfaceIntegrations.F90\
		solver/ALEUtils.F90\
		solver/residuals.F90\
		solver/smoothers.F90\
		solver/multiGrid.F90\
		solver/solvers.F90\
		\
		turbulence/turbMod.F90\
		turbulence/turbBCRoutines.F90\
		turbulence/turbUtils.F90\
		turbulence/turbCurveFits.F90\
		turbulence/turbUtils.F90\
		turbulence/sa.F90\
		turbulence/kw.F90\
		turbulence/kt.F90\
		turbulence/SST.F90\
		turbulence/vf.F90\
		turbulence/turbAPI.F90\
		\
		utils/flowUtils.F90\
		utils/utils.F90\
		utils/haloExchange.F90\
		utils/surfaceUtils.F90\
		utils/sorting.F90\
		utils/signals.F90\
		\
		wallDistance/wallDistance.F90\
		\
		warping/warping.F90\
		warping/getForces.F90\
		warping/getAreas.f90\

cFiles = 	metis-4.0/metis.c\
		utils/myIsNaNC.c\
		utils/connect_signals.c\
		utils/signals.c\
		adjoint/initPETScWrap.c\
		partitioning/metisInterface.c\
		adjoint/ADFirstAidKit/adStack.c\
		adjoint/ADFirstAidKit/dpStack.c\

f77Files = 	adjoint/ADFirstAidKit/adBuffer.f\
		adjoint/ADFirstAidKit/debugAD.f\


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

# Finally convert all source files to .o 
OFILES=$(f90Files:%.f90=%.o) $(F90Files:%.F90=%.o) $(cFilesNoDir:%.c=%.o) $(f77FilesNoDir:%.f=%.o)

# Define the file name for the auto-generated dependency file.
DEP_FILE=sumb_project.dep

# Not real targets are "phony"
.PHONY: sources lib

default: all

# The 'source's makes all the .o files. It also depends on the dependence file "DEP_FILE". 
sources: $(OFILES) $(DEP_FILE) 

# The .dep file depends on the source files, so it automatically gets updated
# when you change your source
$(DEP_FILE): $(fortranFiles)
	@echo "Making dependencies!\n"
	python fort_depend.py -w -o $(DEP_FILE) -f $(fortranFiles)

# Include the automatically generated file.
include $(DEP_FILE)

# The 'lib' target builds the static archive. It depends on "sources".
lib: sources
	$(AR) $(AR_FLAGS) libsumb.a $(OFILES)

# These are the extra python objects for the python library
PYTHON_OBJECTS = fortranobject.o\
		 libsumbmodule.o\
		 libsumb-f2pywrappers2.o

# For all (driven by default)...we depend on having the library 'lib'
# up to date as well changes to the pyf file. This is where we execute
# the "PMAKE" command so everything is built in parallel.
all: 
	$(PMAKE) python

python: lib ../f2py/sumb.pyf
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
	$(eval F2PY_ROOT = $(shell $(PYTHON) ../f2py/get_f2py.py))
	@echo "#------------------------------------------------------#"
	@echo f2py root directory: $(F2PY_ROOT)
	@echo "#------------------------------------------------------#"

# Check the consistency of the pyf file. All carachters (except comments and preprocessor tags) need to be lower case!!!
	python ../f2py/checkPyfForUpperCase.py

# Run the preprocessor on the pyf file:
	python ../f2py/pyf_preprocessor.py real ../f2py/sumb.pyf 

# Run f2py to get sumbmodule.c and sumb-f2pywrapper2.f90
	f2py ../f2py/sumb.pyf.autogen

# Compile c wrapper. Don't use CC_ALL_FLAGS...PETSc wil F-up this command.
	$(CC) $(CC_ALL_FLAGS) $(PYTHON_INCLUDES) -I$(NUMPY_INCLUDES) \
	-I$(F2PY_ROOT)/src -c libsumbmodule.c

# Compile fortranobject needed by all f2py modules
	$(CC) $(CC_ALL_FLAGS) $(PYTHON_INCLUDES) -I$(NUMPY_INCLUDES) \
	-c $(F2PY_ROOT)/src/fortranobject.c -o fortranobject.o

# Compiled f2py-generated wrapper file
	$(FF90) $(FF90_ALL_FLAGS) -I./ -c libsumb-f2pywrappers2.f90

# Final Link:
	$(FF90) -shared $(PYTHON_OBJECTS) libsumb.a $(LINKER_ALL_FLAGS) -o libsumb.so
	python ../f2py/importTest.py
	mv libsumb.so ../../python



