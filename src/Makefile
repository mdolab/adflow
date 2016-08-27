# This is the main makefile for *all* of sumb code. 

%.o: %.f90
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $*.o
	@echo
	@echo "        --- Compiled $*.f90 successfully ---"
	@echo

%.o: %.f
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $*.o
	@echo
	@echo "        --- Compiled $*.f90 successfully ---"
	@echo

%.o: %.F90
	$(FF90) $(FF90_ALL_FLAGS) -c $< -o $*.o
	@echo
	@echo "        --- Compiled $*.F90 successfully ---"
	@echo

%.o: %.c
	$(CC) $(CC_ALL_FLAGS) -c $< -o $*.o
	@echo
	@echo "        --- Compiled $*.c successfully ---"
	@echo


# Note that PETSC_CC_INCLUDES and PETSC_LIB get defined by the petsc
# import.

include ../config.mk
FF90_ALL_FLAGS   = $(FF90_FLAGS) $(CGNS_INCLUDE_FLAGS) -I. \
		   $(PETSC_CC_INCLUDES) $(FF90_PRECISION_FLAGS)

CC_ALL_FLAGS     = $(C_FLAGS) -I./c_defines  -I./metis-4.0  $(PETSC_CC_INCLUDES) \
		   $(CC_PRECISION_FLAGS)

LINKER_ALL_FLAGS = -L../lib -lsumb $(CGNS_LINKER_FLAGS) $(PETSC_LIB) $(LINKER_FLAGS)

# List out all the objects in order:

FF90_OBJECTS_0 =	modules/complexify.o\
			modules/precision.o\
			modules/su_cgns.o\
			modules/constants.o\
			modules/killSignals.o\
			modules/communication.o\
			modules/block.o\
			modules/blockPointers.o\
			modules/kd_tree.o\

FF90_OBJECTS_1 =	modules/IOModule.o\
			modules/costFunctions.o\
			modules/surfaceFamilies.o\
			modules/adtData.o\
			modules/kd_tree.o\
			modules/BCPointers.o\
			modules/BCPointers_b.o\
			modules/BCPointers_d.o\
			modules/BCPointers_fast_b.o\
			modules/wallDistanceData.o\
			modules/ADjointVars.o\
			modules/ADjointPETSc.o\
			modules/cgnsGrid.o\
			modules/cgnsNames.o\
			modules/commMixing.o\
			modules/commSliding.o\
			modules/couplerParam.o\
			modules/CpCurveFits.o\
			modules/extraOutput.o\
			modules/flowVarRefState.o\
			modules/inputParam.o\
			modules/interfaceGroups.o\
			modules/iteration.o\
			modules/monitor.o\
			modules/paramTurb.o\
			modules/section.o\
			modules/stencils.o\
			modules/diffSizes.o\
			modules/overset.o\
			modules/liftDistributionData.o\


FF90_OBJECTS_2= 	ADT/adtUtils.o\
			utils/utils.o\
			adjoint/outputForward/utils_d.o\
			adjoint/outputReverse/utils_b.o\
			adjoint/outputReverseFast/utils_fast_b.o\
			utils/myIsNaNC.o\

FF90_OBJECTS_3 =utils/sorting.o\
		utils/flowUtils.o\
		adjoint/outputForward/flowutils_d.o\
		adjoint/outputReverse/flowutils_b.o\
		adjoint/outputReverseFast/flowutils_fast_b.o\
		ADT/adtBuild.o\
		ADT/adtLocalSearch.o\
		bcdata/BCData.o\
		partitioning/partitionMod.o\
		partitioning/metisInterface.o\
		warping/warping.o\
		solver/BCRoutines.o\
		adjoint/outputForward/bcroutines_d.o\
		adjoint/outputReverse/bcroutines_b.o\
		adjoint/outputReverseFast/bcroutines_fast_b.o\

FF90_OBJECTS_4 =utils/haloExchange.o\
		utils/surfaceUtils.o\
		ADT/adtSearch.o\
		initFlow/variableReading.o\
		partitioning/loadBalance.o\
		partitioning/gridChecking.o\
		partitioning/readCGNSGrid.o\
		adjoint/outputReverseFast/BCExtra.o\

FF90_OBJECTS_5 =inputParam/inputParamRoutines.o\
		ADT/adtAPI.o\
		NKSolver/NKSolvers.o\
		partitioning/partitioning.o\
		preprocessing/bcHalo.o\
		preprocessing/coarse1to1Subface.o\
		preprocessing/haloList.o\
		preprocessing/indirectHalo.o\
		preprocessing/periodicInfo.o\
		output/outputMod.o\
		turbulence/turbMod.o\
		turbulence/saModule.o\
		adjoint/outputForward/samodule_d.o\
		adjoint/outputReverse/samodule_b.o\
		adjoint/outputReverseFast/samodule_fast_b.o\
		overset/stringOps.o\
		overset/cartMesh.o\

FF90_OBJECTS_6 =wallDistance/wallDistance.o\
		adjoint/outputForward/walldistance_d.o\
		adjoint/outputReverse/walldistance_b.o\
		initFlow/initializeFlow.o\
		adjoint/outputForward/initializeflow_d.o\
		adjoint/outputReverse/initializeflow_b.o\
		slidingComm/localSubfaces.o\
		slidingComm/mixingData.o\
		slidingComm/thisSlide.o\
		slidingComm/tmpSliding.o\
		slidingComm/updateComm.o\

# These are the solver ones
FF90_OBJECTS_7 =slidingComm/bsearchSortedDonorCommType.o\
		metis-4.0/metis.o\
		solver/applyAllBC.o\
		solver/convergenceHeader.o\
		solver/adjustInflowAngle.o\
		solver/convergenceInfo.o\
		solver/computedwD3ADI.o\
		solver/computeLamViscosity.o\
		solver/computeUtau.o\
		solver/coorUtilsALE.o\
		solver/DADISmoother.o\
		solver/derivativeRotMatrixRigid.o\
		solver/executeMGCycle.o\
		solver/farFieldDrag.o\
		solver/farFieldInducedDrag.o\
		solver/forcesAndMoments.o\
		solver/forcesAndMomentsZipper.o\
		solver/farFieldDrag.o\
		solver/geomUtilsALE.o\
		solver/getDirAngle.o\
		solver/getDirVector.o\
		solver/gridVelocitiesCoarseLevels.o\
		solver/gridVelocitiesFineLevel.o\
		solver/gridVelocitiesFineLevel_ALE.o\
		solver/initres.o\
		solver/inviscidCentralFlux.o\
		solver/inviscidDissFluxMatrix.o\
		solver/inviscidDissFluxMatrixCoarse.o\
		solver/inviscidDissFluxScalar.o\
		solver/inviscidDissFluxScalarCoarse.o\
		solver/inviscidUpwindFlux.o\
		solver/killFunctions.o\
		solver/connect_signals.o\
		solver/maxEddyv.o\
		solver/maxHdiffMach.o\
		solver/normalVelocities.o\
		solver/residual.o\
		solver/residual2.o\
		solver/residualAveraging.o\
		solver/rungeKuttaSmoother.o\
		solver/setCoefTimeIntegrator.o\
		solver/setCornerRowHalos.o\
		solver/setCorrectionsCoarseHalos.o\
		solver/setCycleStrategy.o\
		solver/shiftSolution.o\
		solver/slipVelocities.o\
		solver/slipVelocities_ALE.o\
		solver/solveState.o\
		solver/solver.o\
		solver/solverSteady.o\
		solver/solverUnsteadyExplicitRK.o\
		solver/solverUnsteadyBDF.o\
		solver/solverUnsteady_ALE.o\
		solver/solverUnsteadyWrap.o\
		solver/timeStep.o\
		solver/transferToCoarseGrid.o\
		solver/transferToFineGrid.o\
		solver/unsteadyHeader.o\
		solver/utauWF.o\
		solver/vectorRotation.o\
		solver/viscousFlux.o\
		preprocessing/coarse1to1Subface.o\
		preprocessing/allocMemHaloList.o\
		preprocessing/cellRangeSubface.o\
		preprocessing/checkCoarse1to1.o\
		preprocessing/checkSymmetry.o\
		preprocessing/closestDirectHalos.o\
		preprocessing/coarseDonorInfo.o\
		preprocessing/createCoarseBlocks.o\
		preprocessing/determineCommPattern.o\
		preprocessing/determineFaceHalos.o\
		preprocessing/determineIndirectHalos.o\
		preprocessing/determineNumberOfHalos.o\
		preprocessing/determineNcellGlobal.o\
		preprocessing/determinePeriodicFaces.o\
		preprocessing/exchangeCoor.o\
		preprocessing/exchangeCoor_d.o\
		preprocessing/exchangeCoor_b.o\
		preprocessing/faceRotationMatrices.o\
		preprocessing/finalCommStructures.o\
		preprocessing/indirectHalosPerLevel.o\
		preprocessing/init2ndLevelCellHalos.o\
		preprocessing/setSurfaceFamilyInfo.o\
		preprocessing/mdUpdateRoutines.o\
		preprocessing/metric.o\
		preprocessing/preprocessing.o\
		preprocessing/qsortHaloListType.o\
		preprocessing/qsortIndHaloType.o\
		preprocessing/qsortPeriodicSubfacesHaloType.o\
		preprocessing/setFamilyInfoFaces.o\
		preprocessing/setGlobalCellsAndNodes.o\
		preprocessing/setPorosities.o\
		preprocessing/shiftCoorAndVolumes.o\
		preprocessing/unitVectorsInAxialPlane.o\
		preprocessing/update1to1Coarse.o\
		preprocessing/updateCoorFineMesh.o\
		preprocessing/viscSubfaceInfo.o\
		preprocessing/xhalo.o\
		preprocessing/setReferenceVolume.o\
		preprocessing/bsearchCGNSPeriodicType.o\
		preprocessing/determinePeriodicData.o\
		output/copyDataBuf.o\
		output/computeIsoVariable.o\
		output/describeScheme.o\
		output/gridFileNamesWrite.o\
		output/isoSurfNames.o\
		output/liftDistribution.o\
		output/marchingCubesData.o\
		output/numberOfSurfSolVariables.o\
		output/numberOfVolSolVariables.o\
		output/numberOfIsoVariables.o\
		output/openCGNSVolumeSol.o\
		output/releaseHelpVariablesWriting.o\
		output/setHelpVariablesWriting.o\
		output/storeOldSolInBuffer.o\
		output/storeSolInBuffer.o\
		output/storeSurfsolInBuffer.o\
		output/surfSolFileNamesWrite.o\
		output/surfSolNames.o\
		output/updateParamfile.o\
		output/volSolFileNamesWrite.o\
		output/volSolNames.o\
		output/writeCGNSConvInfo.o\
		output/writeCGNSGridFile.o\
		output/writeCGNSGridFrame.o\
		output/writeCGNSHeader.o\
		output/writeCGNSKomegaModifiedInfo.o\
		output/writeCGNSKomegaWilcoxInfo.o\
		output/writeCGNSKtauInfo.o\
		output/writeCGNSMenterSSTInfo.o\
		output/writeCGNSReferenceState.o\
		output/writeCGNSSaInfo.o\
		output/writeCGNSSaeInfo.o\
		output/writeCGNSSurfaceSol.o\
		output/writeCGNSTimeHistory.o\
		output/writeCGNSV2fInfo.o\
		output/writeCGNSVolumeSol.o\
		output/writeCoorCGNSZone.o\
		output/writeInterpolFromSpectral.o\
		output/writeIsoSurface.o\
		output/writeSol.o\
		output/writeSolCGNSZone.o\
		output/writeSurfsolCGNSZone.o\
		output/writeUnsteadyFromSpectral.o\
		slidingComm/cumulativeNSendReceives.o\
		slidingComm/cylCoorNodesOnBlockFace.o\
		slidingComm/determineNSlices.o\
		slidingComm/determineRotationMatrices.o\
		slidingComm/determineSections.o\
		slidingComm/haloType.o\
		slidingComm/initMemMixingPlane.o\
		slidingComm/initMemSliding.o\
		slidingComm/initThisSlide.o\
		slidingComm/interpolateSlide.o\
		slidingComm/localNodesAndQuads.o\
		slidingComm/localSurfaceGrids.o\
		slidingComm/mixingDonorInterpol.o\
		slidingComm/mixingHaloInterpol.o\
		slidingComm/mixingIntervals.o\
		slidingComm/mixingPlaneComm.o\
		slidingComm/mixingPlaneInterpol.o\
		slidingComm/mySubfacesSlide.o\
		slidingComm/nLocalNodesAndQuads.o\
		slidingComm/procSliding.o\
		slidingComm/qsortMixingIntervalType.o\
		slidingComm/qsortSortedDonorCommType.o\
		slidingComm/releaseMemSliding.o\
		slidingComm/slidingComm.o\
		slidingComm/slidingIntersections.o\
		slidingComm/slidingMesh.o\
		slidingComm/slidingMeshGroups.o\
		slidingComm/updateExternalCommSlide.o\
		slidingComm/updateInterpolSendBuf.o\
		slidingComm/updateLocalCommSlide.o\
		slidingComm/updateRotationInfo.o\
		slidingComm/updateSlidingAllLevels.o\
		slidingComm/updateSlidingCommPattern.o\
		slidingComm/volumeInterpol.o\
		turbulence/applyAllTurbBC.o\
		turbulence/bcEddyNoWall.o\
		turbulence/bcEddyWall.o\
		turbulence/bcTurbFarfield.o\
		turbulence/bcTurbInflow.o\
		turbulence/bcTurbInterface.o\
		turbulence/bcTurbOutflow.o\
		turbulence/bcTurbSymm.o\
		turbulence/bcTurbTreatment.o\
		turbulence/bcTurbWall.o\
		turbulence/computeEddyViscosity.o\
		turbulence/curveFit.o\
		turbulence/f1SST.o\
		turbulence/initCurveFitDataKtau.o\
		turbulence/initCurveFitDataKw.o\
		turbulence/initCurveFitDataKwMod.o\
		turbulence/initCurveFitDataSa.o\
		turbulence/initCurveFitDataSae.o\
		turbulence/initCurveFitDataSST.o\
		turbulence/initCurveFitDataVf.o\
		turbulence/initKOmega.o\
		turbulence/keSolve.o\
		turbulence/kt.o\
		turbulence/ktCDterm.o\
		turbulence/ktSolve.o\
		turbulence/kw.o\
		turbulence/kwCDterm.o\
		turbulence/kwSolve.o\
		turbulence/prodKatoLaunder.o\
		turbulence/prodSmag2.o\
		turbulence/prodWmag2.o\
		turbulence/saNuFromEddyRatio.o\
		turbulence/SST.o\
		turbulence/SSTSolve.o\
		turbulence/tdia3.o\
		turbulence/turb2ndHalo.o\
		turbulence/turbAdvection.o\
		turbulence/turbBCNSWall.o\
		turbulence/turbResidual.o\
		turbulence/turbSolveSegregated.o\
		turbulence/unsteadyTurbSpectral.o\
		turbulence/unsteadyTurbTerm.o\
		turbulence/vf.o\
		turbulence/vfScale.o\
		turbulence/vfSolve.o\
		turbulence/writeLoglaw.o\
		adjoint/outputForward/allnodalgradients_d.o\
		adjoint/outputForward/viscousflux_d.o \
		adjoint/outputForward/viscousfluxapprox_d.o\
		adjoint/outputForward/invisciddissfluxmatrix_d.o \
		adjoint/outputForward/computelamviscosity_d.o\
		adjoint/outputForward/computepressuresimple_d.o\
		adjoint/outputForward/timestep_block_d.o\
		adjoint/outputForward/invisciddissfluxscalar_d.o\
		adjoint/outputForward/invisciddissfluxscalarapprox_d.o\
		adjoint/outputForward/invisciddissfluxmatrixapprox_d.o\
		adjoint/outputForward/inviscidupwindflux_d.o  \
		adjoint/outputForward/block_res_d.o  \
		adjoint/outputForward/residual_block_d.o \
		adjoint/outputForward/inviscidcentralflux_d.o\
		adjoint/outputForward/saeddyviscosity_d.o\
		adjoint/outputForward/xhalo_block_d.o\
		adjoint/outputForward/metric_block_d.o\
		adjoint/outputForward/gridvelocitiesfinelevel_block_d.o\
		adjoint/outputForward/normalvelocities_block_d.o\
		adjoint/outputForward/adjustinflowangle_d.o\
		adjoint/outputForward/forcesandmoments_d.o\
		adjoint/outputForward/getdirangle_d.o\
		adjoint/outputForward/getdirvector_d.o\
		adjoint/outputForward/vectorrotation_d.o\
		adjoint/outputForward/computeeddyviscosity_d.o\
		adjoint/outputForward/sanuknowneddyratio_d.o\
		adjoint/outputForward/prodkatolaunder_d.o\
		adjoint/outputForward/prodwmag2_d.o\
		adjoint/outputForward/prodsmag2_d.o\
		adjoint/outputForward/turbadvection_d.o\
		adjoint/outputForward/unsteadyturbterm_d.o\
		adjoint/outputForward/slipvelocitiesfinelevel_block_d.o\
		adjoint/outputForward/getcostfunction2_d.o\
		adjoint/outputForward/bcturbwall_d.o\
		adjoint/outputForward/bcturbtreatment_d.o\
		adjoint/outputForward/applyallturbbcthisblock_d.o\
		adjoint/outputForward/bceddynowall_d.o\
		adjoint/outputForward/bceddywall_d.o\
		adjoint/outputForward/bcturbfarfield_d.o\
		adjoint/outputForward/bcturbinterface_d.o\
		adjoint/outputForward/turb2ndhalo_d.o\
		adjoint/outputForward/boundarynormals_d.o\
		adjoint/outputForward/volume_block_d.o\
		adjoint/outputForward/derivativerotmatrixrigid_d.o\
		adjoint/outputReverse/allnodalgradients_b.o\
		adjoint/outputReverse/adjustinflowangle_b.o\
		adjoint/outputReverse/getdirvector_b.o\
		adjoint/outputReverse/vectorrotation_b.o\
		adjoint/outputReverse/computelamviscosity_b.o\
		adjoint/outputReverse/computepressuresimple_b.o\
		adjoint/outputReverse/getcostfunction2_b.o\
		adjoint/outputReverse/timestep_block_b.o\
		adjoint/outputReverse/metric_block_b.o\
		adjoint/outputReverse/xhalo_block_b.o\
		adjoint/outputReverse/residual_block_b.o\
		adjoint/outputReverse/invisciddissfluxscalar_b.o\
		adjoint/outputReverse/block_res_b.o  \
		adjoint/outputReverse/viscousflux_b.o \
		adjoint/outputReverse/viscousfluxapprox_b.o\
		adjoint/outputReverse/inviscidcentralflux_b.o\
		adjoint/outputReverse/saeddyviscosity_b.o\
		adjoint/outputReverse/computeeddyviscosity_b.o\
		adjoint/outputReverse/unsteadyturbterm_b.o\
		adjoint/outputReverse/prodkatolaunder_b.o\
		adjoint/outputReverse/prodwmag2_b.o\
		adjoint/outputReverse/prodsmag2_b.o\
		adjoint/outputReverse/turbadvection_b.o\
		adjoint/outputReverse/forcesandmoments_b.o\
		adjoint/outputReverse/sanuknowneddyratio_b.o\
		adjoint/outputReverse/getcostfunction2_b.o\
		adjoint/outputReverse/invisciddissfluxscalarapprox_b.o\
		adjoint/outputReverse/invisciddissfluxmatrix_b.o\
		adjoint/outputReverse/inviscidupwindflux_b.o\
		adjoint/outputReverse/invisciddissfluxmatrixapprox_b.o\
		adjoint/outputReverse/applyallturbbcthisblock_b.o\
		adjoint/outputReverse/bcturbtreatment_b.o\
		adjoint/outputReverse/turb2ndhalo_b.o\
		adjoint/outputReverse/bceddywall_b.o\
		adjoint/outputReverse/bceddynowall_b.o\
		adjoint/outputReverse/bcturbwall_b.o\
		adjoint/outputReverse/bcturbinterface_b.o\
		adjoint/outputReverse/bcturbfarfield_b.o\
		adjoint/outputReverse/boundarynormals_b.o\
		adjoint/outputReverse/volume_block_b.o\
		adjoint/outputReverseFast/allnodalgradients_fast_b.o\
		adjoint/outputReverseFast/computelamviscosity_fast_b.o\
		adjoint/outputReverseFast/timestep_block_fast_b.o\
		adjoint/outputReverseFast/invisciddissfluxscalar_fast_b.o\
		adjoint/outputReverseFast/invisciddissfluxmatrix_fast_b.o\
		adjoint/outputReverseFast/viscousflux_fast_b.o\
		adjoint/outputReverseFast/inviscidcentralflux_fast_b.o\
		adjoint/outputReverseFast/computepressuresimple_fast_b.o\
		adjoint/outputReverseFast/BCExtra.o\
		adjoint/outputReverseFast/prodsmag2_fast_b.o\
		adjoint/outputReverseFast/prodwmag2_fast_b.o\
		adjoint/outputReverseFast/prodkatolaunder_fast_b.o\
		adjoint/outputReverseFast/turbadvection_fast_b.o\
		adjoint/outputReverseFast/saeddyviscosity_fast_b.o\
		adjoint/residualInput/block_res.o\
		adjoint/residualInput/metric_block.o\
		adjoint/residualInput/xhalo_block.o\
		adjoint/residualInput/viscousFluxApprox.o\
		adjoint/residualInput/inviscidDissFluxScalarApprox.o\
		adjoint/residualInput/inviscidDissFluxMatrixApprox.o\
		adjoint/residualInput/computePressureSimple.o\
		adjoint/residualInput/costFunction2.o\
		adjoint/adjointPythonInterfaceFunctions.o\
		adjoint/computeAeroCoef.o\
		adjoint/createPETScVars.o\
		adjoint/destroyPETScVars.o\
		adjoint/preallocation.o\
		adjoint/getSolution.o\
		adjoint/initializePETSc.o\
		adjoint/finalizePETSc.o\
		adjoint/preprocessingADjoint.o\
		adjoint/setupPETScKsp.o\
		adjoint/updateFlow.o\
		adjoint/updateReferencePoint.o\
		adjoint/updateRotationRate.o\
		adjoint/matrixFreeRoutines.o\
		adjoint/alloc_derivative_values.o\
		adjoint/dealloc_derivative_values.o\
		adjoint/allocPCMem.o\
		adjoint/determineStencil.o\
		adjoint/setupAllResidualMatricesFwd.o\
		adjoint/setupColoring.o\
		adjoint/setPointers_d.o\
		adjoint/setDiffSizes.o\
		adjoint/zeroADSeeds.o\
		adjoint/setFDReference.o\
		adjoint/referenceShockSensor.o\
		adjoint/PCMatMult.o\
		adjoint/solveADjointTransposePETSc.o\
		adjoint/initPETScWrap.o\
		adjoint/saveAdjoint.o\
		adjoint/setupStateResidualMatrix.o\
		adjoint/setupPCMatrix.o\
		adjoint/ADFirstAidKit/adBuffer.o\
		adjoint/ADFirstAidKit/adStack.o\
		adjoint/ADFirstAidKit/dpStack.o\
		adjoint/ADFirstAidKit/debugAD.o\
		warping/getForces.o\
		warping/getAreas.o\
		overset/computeInterpolationParallel.o\
		overset/buildClusterWalls.o\
		overset/createZipperMesh.o\
		overset/bowtieAndIsolationElimination.o\
		overset/determineClusters.o \
		overset/determineClusterAreas.o \
		overset/determineClusterMarchDist.o\
		overset/checkOverset.o \
		overset/computeCellWallPoint.o\
		overset/exchangeSurfaceDelta.o\
		overset/floodInteriorCells.o\
		overset/oversetLoadBalance.o\
		overset/qsortFringeType.o\
		overset/qsortEdgeType.o\
		overset/qsortPocketEdgeType.o\
		overset/binSearchPocketEdgeType.o\
		overset/binSearchNodes.o\
		overset/initializeOBlock.o\
		overset/initializeOWall.o\
		overset/initializeFringes.o\
		overset/initializeOFringes.o\
		overset/determineDonors.o\
		overset/packingRoutines.o\
		overset/fringeSearch.o\
		overset/oversetUtilities.o\
		overset/irregularCellCorrection.o\
		overset/fringeReduction.o\
		overset/finalOversetCommStructures.o\
		overset/setIBlankArray.o\
		overset/exchangeFringes.o\
		overset/surfaceCorrection.o\
		overset/emptyOversetComm.o\
		overset/wallSearch.o\
		overset/initBCDataiBlank.o\
		overset/exchangeSurfaceIBlanks.o\
		overset/determineZipperWallAssociation.o\
		overset/makeBoundaryStrings.o\
		overset/writeDualGrid.o\
		overset/exchangeStatus.o\
		overset/exchangeStatusTranspose.o\
		overset/surfaceDeviation.o\

PYTHON_OBJECTS = fortranobject.o\
		 libsumbmodule.o\
		 libsumb-f2pywrappers2.o

default: all

all:
	$(MAKE) src0
	$(PMAKE) src1
	$(PMAKE) src2
	$(PMAKE) src3
	$(PMAKE) src4
	$(PMAKE) src5
	$(PMAKE) src6
	$(PMAKE) src7
	$(AR) $(AR_FLAGS) ../lib/libsumb.a $(FF90_OBJECTS_0) $(FF90_OBJECTS_1) \
	$(FF90_OBJECTS_2) $(FF90_OBJECTS_3) $(FF90_OBJECTS_4) $(FF90_OBJECTS_5) \
	$(FF90_OBJECTS_6) $(FF90_OBJECTS_7)

# --------------------- f2py python module building ---------------

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
	$(eval F2PY_ROOT = $(shell $(PYTHON) f2py/get_f2py.py))
	@echo "#------------------------------------------------------#"
	@echo f2py root directory: $(F2PY_ROOT)
	@echo "#------------------------------------------------------#"

# Check the consistency of the pyf file. All carachters (except comments and preprocessor tags) need to be lower case!!!
	python f2py/checkPyfForUpperCase.py

# Run the preprocessor on the pyf file:
	python f2py/pyf_preprocessor.py real f2py/sumb.pyf 

# Run f2py to get sumbmodule.c and sumb-f2pywrapper2.f90
	f2py f2py/sumb.pyf.autogen

# Compile c wrapper. Don't use CC_ALL_FLAGS...PETSc wil F-up this command.

	$(CC) $(CC_ALL_FLAGS) $(PYTHON_INCLUDES) -I$(NUMPY_INCLUDES) \
	-I$(F2PY_ROOT)/src -c libsumbmodule.c

# Compile fortranobject needed by all f2py modules
	$(CC) $(CC_ALL_FLAGS) $(PYTHON_INCLUDES) -I$(NUMPY_INCLUDES) \
	-c $(F2PY_ROOT)/src/fortranobject.c -o fortranobject.o

# Compiled f2py-generated wrapper file
	$(FF90) $(FF90_ALL_FLAGS) -I./ -c libsumb-f2pywrappers2.f90

# Final Link:
	$(FF90) -shared $(PYTHON_OBJECTS) $(LINKER_ALL_FLAGS) -o libsumb.so
	python f2py/importTest.py
	mv libsumb.so ../../../python
	-rm *.o
	-rm *.f90
	-rm *.c

clean: 
	rm -f *.o *.mod

src0: $(FF90_OBJECTS_0)
src1: $(FF90_OBJECTS_1)
src2: $(FF90_OBJECTS_2)
src3: $(FF90_OBJECTS_3)
src4: $(FF90_OBJECTS_4)
src5: $(FF90_OBJECTS_5)
src6: $(FF90_OBJECTS_6)
src7: $(FF90_OBJECTS_7)
