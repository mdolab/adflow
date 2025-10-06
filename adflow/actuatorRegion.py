from abc import ABC, abstractmethod
from typing import Tuple
from types import ModuleType
import os

from mpi4py import MPI

import numpy as np
import numpy.typing as npt
import scipy.interpolate as interpolate


from baseclasses import AeroProblem


class AbstractActuatorRegion(ABC):
    def __init__(self, centerPoint: npt.NDArray, thrustVector: npt.NDArray, thrust: float, heat: float):
        if centerPoint.shape != (3,):
            raise ValueError('"centerPoint" must have shape "(3,)"')
        if thrustVector.shape != (3,):
            raise ValueError('"thrustVector" must have shape "(3,)"')
        if np.linalg.norm(thrustVector) == 0:
            raise ValueError('"trustVector" can not have a length of "0"')
        if thrust < 0:
            raise ValueError('"thrust" must not be smaller than 0.')
        if heat < 0:
            raise ValueError('"heat" must not be smaller than 0.')

        self.centerPoint = centerPoint
        self.thrustVector = thrustVector / np.linalg.norm(thrustVector)
        self.iRegion = -1
        self.nLocalCells = -1
        self.familyName = ""

        self._heatBcName = "Heat"
        self._thrustBcName = "Thrust"

        self._boundaryConditions = {
            self._thrustBcName: thrust,
            self._heatBcName: heat,
        }

    def setBCValue(self, bcName, bcValue):
        if bcName not in self._boundaryConditions.keys():
            raise ValueError(
                f'"{bcName}" is not a valid Boundary Condition. Valid Conditions are: {list(self._boundaryConditions.keys())}'
            )

        self._boundaryConditions[bcName] = bcValue

    def _computeScalingConstant(self, cellValues, totalValue, comm):
        computedTotalvalue = comm.allreduce(np.sum(cellValues), op=MPI.SUM)

        if np.isclose(computedTotalvalue, 0):
            return 0.0

        return totalValue / computedTotalvalue

    @abstractmethod
    def tagActiveCells(
        self, distance2axis: npt.NDArray, distance2plane: npt.NDArray, tangent: npt.NDArray
    ) -> npt.NDArray:
        flags = np.zeros_like(distance2axis)
        return flags

    @abstractmethod
    def computeCellForceVector(
        self,
        distance2axis: npt.NDArray,
        distance2plane: npt.NDArray,
        tangent: npt.NDArray,
        cellVolume: npt.NDArray,
        totalVolume: float,
        comm: MPI.Intracomm,
    ) -> npt.NDArray:
        force = np.zeros((3, len(distance2axis)))
        return force

    @abstractmethod
    def computeCellHeatVector(
        self,
        distance2axis: npt.NDArray,
        distance2plane: npt.NDArray,
        tangent: npt.NDArray,
        cellVolume: npt.NDArray,
        totalVolume: float,
        comm: MPI.Intracomm,
    ) -> npt.NDArray:
        heat = np.zeros_like(distance2axis)
        return heat


class ActuatorRegionHandler:
    def __init__(self, adflow_fortran: ModuleType, comm):
        self._comm = comm
        self._adflow_fortran = adflow_fortran

        self._actuatorRegions = list()

    def addActuatorRegion(
        self,
        actuatorRegion: AbstractActuatorRegion,
        familyName: str,
        familyID: int,
        relaxStart: float,
        relaxEnd: float,
    ):
        # compute the distance of each cell to the AR plane and axis
        ncells = self._adflow_fortran.adjointvars.ncellslocal[0]

        distance2plane, distance2axis, tangent = self._adflow_fortran.actuatorregion.computeinitialspatialmetrics(
            actuatorRegion.centerPoint, actuatorRegion.thrustVector, ncells
        )
        tangent = tangent.T

        # tag the active cells
        flag = actuatorRegion.tagActiveCells(distance2axis, distance2plane, tangent)
        iRegion, nLocalCells = self._adflow_fortran.actuatorregion.addactuatorregion(
            flag,
            familyName,
            familyID,
            relaxStart,
            relaxEnd,
        )
        actuatorRegion.iRegion = iRegion
        actuatorRegion.nLocalCells = nLocalCells
        actuatorRegion.familyName = familyName

        # book keep the new region
        self._actuatorRegions.append(actuatorRegion)

    def updateActuatorRegionsBC(self, aeroproblem: AeroProblem):
        for actuatorRegion in self._actuatorRegions:

            # set BC values
            for (bcName, familyName), bcValue in aeroproblem.bcVarData.items():
                if familyName != actuatorRegion.familyName:
                    continue
                actuatorRegion.setBCValue(bcName, bcValue)

            # compute local metrics
            distance2plane, distance2axis, tangent, volume = self._adflow_fortran.actuatorregion.computespatialmetrics(
                actuatorRegion.iRegion,
                actuatorRegion.nLocalCells,
                actuatorRegion.centerPoint,
                actuatorRegion.thrustVector,
            )
            tangent = tangent.T

            totalVolume = self._comm.allreduce(np.sum(volume), op=MPI.SUM)

            # compute the source terms
            force = actuatorRegion.computeCellForceVector(
                distance2axis, distance2plane, tangent, volume, totalVolume, self._comm
            )
            heat = actuatorRegion.computeCellHeatVector(
                distance2axis, distance2plane, tangent, volume, totalVolume, self._comm
            )

            # skip this proc if no cells are active
            if actuatorRegion.nLocalCells == 0:
                continue

            # apply the source terms in Fortran
            self._adflow_fortran.actuatorregion.populatebcvalues(
                actuatorRegion.iRegion,
                actuatorRegion.nLocalCells,
                force.T,
                heat.T,
            )

    def writeActuatorRegions(self, fileName: str, outputDir: str):
        # Join to get the actual filename root
        fileName = os.path.join(outputDir, fileName)

        # Ensure extension is .plt even if the user didn't specify
        fileName, ext = os.path.splitext(fileName)
        fileName += ".plt"

        # just call the underlying fortran routine
        self._adflow_fortran.actuatorregion.writeactuatorregions(fileName)


class CircularActuatorRegion(AbstractActuatorRegion):
    def __init__(
        self,
        centerPoint: npt.NDArray,
        thrustVector: npt.NDArray,
        innerDiameter: float,
        outerDiameter: float,
        regionDepth: float,
        thrust: float = 0,
        heat: float = 0,
    ):
        super().__init__(centerPoint, thrustVector, thrust, heat)

        if regionDepth <= 0:
            raise ValueError('"depth" must be greater than 0.')
        if innerDiameter < 0:
            raise ValueError('"innerDiameter" must not be smaller than 0.')
        if outerDiameter <= 0:
            raise ValueError('"outerDiameter" must be greater than 0.')
        if outerDiameter <= innerDiameter:
            raise ValueError('"outerDiameter" must be greater than "innerDiameter".')

        self._innerDiameter = innerDiameter
        self._outerDiameter = outerDiameter
        self._regionDepth = regionDepth

    def tagActiveCells(self, distance2axis, distance2plane, tangent):
        flags = np.zeros_like(distance2axis)

        indices = np.logical_and(
            np.logical_and(
                distance2axis >= self._innerDiameter / 2,
                distance2axis <= self._outerDiameter / 2,
            ),
            np.logical_and(distance2plane >= 0, distance2plane <= self._regionDepth),
        )

        flags[indices] = 1

        return flags


class UniformActuatorRegion(CircularActuatorRegion):
    def computeCellForceVector(self, distance2axis, distance2plane, tangent, cellVolume, totalVolume, comm):
        thrustBC = self._boundaryConditions[self._thrustBcName]

        cellThrusts = thrustBC * cellVolume / totalVolume
        scaling_constant = self._computeScalingConstant(cellThrusts, thrustBC, comm)
        force = np.outer(scaling_constant * cellThrusts, -self.thrustVector)

        return force

    def computeCellHeatVector(self, distance2axis, distance2plane, tangent, cellVolume, totalVolume, comm):
        heatBC = self._boundaryConditions[self._heatBcName]

        cellHeats = heatBC * cellVolume / totalVolume
        scaling_constant = self._computeScalingConstant(cellHeats, heatBC, comm)

        return cellHeats * scaling_constant


class BSplineActuatorRegion(CircularActuatorRegion):
    def __init__(
        self,
        centerPoint: npt.NDArray,
        thrustVector: npt.NDArray,
        innerDiameter: float,
        outerDiameter: float,
        regionDepth: float,
        thrustDistribution: Tuple[npt.NDArray, npt.NDArray, float] = (
            np.array([0.0, 0.0, 1.0, 1.0]),
            np.array([0.0, 0.0, 0.0, 0.0]),
            1,
        ),
        tangentDistribution: Tuple[npt.NDArray, npt.NDArray, float] = (
            np.array([0.0, 0.0, 1.0, 1.0]),
            np.array([0.0, 0.0, 0.0, 0.0]),
            1,
        ),
        radialDistribution: Tuple[npt.NDArray, npt.NDArray, float] = (
            np.array([0.0, 0.0, 1.0, 1.0]),
            np.array([0.0, 0.0, 0.0, 0.0]),
            1,
        ),
        heatDistribution: Tuple[npt.NDArray, npt.NDArray, float] = (
            np.array([0.0, 0.0, 1.0, 1.0]),
            np.array([0.0, 0.0, 0.0, 0.0]),
            1,
        ),
        thrust: float = 0,
        heat: float = 0,
    ):
        super().__init__(centerPoint, thrustVector, innerDiameter, outerDiameter, regionDepth, thrust, heat)

        self._thrustSpline = interpolate.BSpline(*thrustDistribution)
        self._tangentSpline = interpolate.BSpline(*tangentDistribution)
        self._radialSpline = interpolate.BSpline(*radialDistribution)
        self._heatSpline = interpolate.BSpline(*heatDistribution)

    def _formSpline(self, distribution: Tuple[npt.NDArray, npt.NDArray, float], scalingConstant: float):
        t = distribution[0]
        c = distribution[1]
        k = distribution[2]

        return interpolate.BSpline(t, c * scalingConstant, k)

    def _computeNormalizedRadius(self, distance2axis):
        normalizedRadius = distance2axis - self._innerDiameter / 2
        normalizedRadius /= self._outerDiameter / 2 - self._innerDiameter / 2
        return normalizedRadius

    def computeCellForceVector(self, distance2axis, distance2plane, tangent, cellVolume, totalVolume, comm):
        thrustBC = self._boundaryConditions[self._thrustBcName]

        normalizedRadius = self._computeNormalizedRadius(distance2axis)
        thrustFactor = thrustBC * cellVolume / totalVolume

        # add axial contribution (thrust)
        cellThrusts = self._thrustSpline(normalizedRadius) * thrustFactor

        # compute a scaling constant such that the summ of the thrust equals the prescribed thrust
        scaling_constant = self._computeScalingConstant(cellThrusts, thrustBC, comm)
        force = np.outer(scaling_constant * cellThrusts, -self.thrustVector)

        # add tangential contribution (swirl)
        force += np.multiply(
            tangent.T, np.array([scaling_constant * self._tangentSpline(normalizedRadius) * thrustFactor])
        ).T

        # add radial contribution
        radius = np.cross(self.thrustVector, tangent)
        force += np.multiply(
            radius.T, np.array([scaling_constant * self._radialSpline(normalizedRadius) * 2 * thrustFactor])
        ).T

        return force

    def computeCellHeatVector(self, distance2axis, distance2plane, tangent, cellVolume, totalVolume, comm):
        heatBC = self._boundaryConditions[self._heatBcName]

        normalizedRadius = self._computeNormalizedRadius(distance2axis)

        # compute a scaling constant such that the summ of the heat equals the prescribed heat

        cellHeats = self._heatSpline(normalizedRadius) * heatBC * cellVolume / totalVolume
        scaling_constant = self._computeScalingConstant(cellHeats, heatBC, comm)

        return cellHeats * scaling_constant
