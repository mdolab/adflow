from abc import ABC, abstractmethod
from typing import Tuple

import numpy as np
import numpy.typing as npt
import scipy.interpolate as interpolate


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
        self.familyName = None


        self._heatBcName = 'Heat'
        self._thrustBcName = 'Thrust'

        self._boundaryConditions = {
                self._thrustBcName: thrust,
                self._heatBcName: heat,
                }

    def setBCValue(self, bcName, bcValue):
        if bcName not in self._boundaryConditions.keys():
            raise ValueError(f'"{bcName}" is not a valid Boundary Condition. Valid Conditions are: {list(self._boundaryConditions.keys())}')
        
        self._boundaryConditions[bcName] = bcValue

    @abstractmethod
    def tagActiveCells(self, distance2axis: npt.NDArray, distance2plane: npt.NDArray, tangent: npt.NDArray) -> npt.NDArray:
        flags = np.zeros_like(distance2axis)
        return flags

    @abstractmethod
    def computeCellForceVector(self, distance2axis: npt.NDArray, distance2plane: npt.NDArray, tangent: npt.NDArray, cellVolume: npt.NDArray, totalVolume: float) -> npt.NDArray:
        force = np.zeros((3, len(distance2axis)))
        return force

    @abstractmethod
    def computeCellHeatVector(self, distance2axis: npt.NDArray, distance2plane: npt.NDArray, tangent: npt.NDArray, cellVolume: npt.NDArray, totalVolume: float) -> npt.NDArray:
        heat = np.zeros_like(distance2axis)
        return heat



class CircularActuatorRegion(AbstractActuatorRegion):
    def __init__(self, 
                 centerPoint: npt.NDArray, 
                 thrustVector: npt.NDArray, 
                 innerDiameter: float, 
                 outerDiameter: float, 
                 regionDepth: float, 
                 thrust: float = 0, 
                 heat: float = 0):
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

        indices =  np.logical_and(
                np.logical_and(
                    distance2axis >= self._innerDiameter/2,
                    distance2axis <= self._outerDiameter/2,
                    ),
                np.logical_and(
                    distance2plane >= 0,
                    distance2plane <= self._regionDepth
                    )
                )

        flags[indices] = 1

        return flags





class UniformActuatorRegion(CircularActuatorRegion):
    def computeCellForceVector(self, distance2axis, distance2plane, tangent, cellVolume, totalVolume):
        force = np.zeros((len(distance2axis), 3))
        thrustBC = self._boundaryConditions[self._thrustBcName]

        force[:, :] = np.outer(
                thrustBC * cellVolume / totalVolume, 
                self.thrustVector
                )

        return force

    def computeCellHeatVector(self, distance2axis, distance2plane, tangent, cellVolume, totalVolume):
        heat = np.zeros_like(distance2axis)
        heatBC = self._boundaryConditions[self._heatBcName]

        heat = heatBC * cellVolume / totalVolume

        return heat





class BSplineActuatorRegion(CircularActuatorRegion):
    def __init__(self, 
                 centerPoint: npt.NDArray, 
                 thrustVector: npt.NDArray, 
                 innerDiameter: float, 
                 outerDiameter: float, 
                 regionDepth: float, 
                 thrustDistribution: Tuple[npt.NDArray, npt.NDArray, float] = (np.array([0., 0., 1., 1.]), np.array([0., 0., 0., 0.]), 1),
                 tangentDistribution: Tuple[npt.NDArray, npt.NDArray, float] = (np.array([0., 0., 1., 1.]), np.array([0., 0., 0., 0.]), 1),
                 radialDistribution: Tuple[npt.NDArray, npt.NDArray, float] = (np.array([0., 0., 1., 1.]), np.array([0., 0., 0., 0.]), 1),
                 heatDistribution: Tuple[npt.NDArray, npt.NDArray, float] = (np.array([0., 0., 1., 1.]), np.array([0., 0., 0., 0.]), 1),
                 thrust: float = 0, 
                 heat: float = 0
                 ):
        super().__init__(centerPoint, thrustVector, innerDiameter, outerDiameter, regionDepth, thrust, heat)

        self._thrustSpline = interpolate.BSpline(*thrustDistribution) 
        self._tangentSpline = interpolate.BSpline(*tangentDistribution) 
        self._radialSpline = interpolate.BSpline(*radialDistribution) 
        self._heatSpline = interpolate.BSpline(*heatDistribution) 

    def _computeNormalizedRadius(self, distance2axis):
        normalizedRadius = (distance2axis - self._innerDiameter/2) / (self._outerDiameter/2 - self._innerDiameter)
        return normalizedRadius


    def computeCellForceVector(self, distance2axis, distance2plane, tangent, cellVolume, totalVolume):
        force = np.zeros((len(distance2axis), 3))
        thrustBC = self._boundaryConditions[self._thrustBcName]

        normalizedRadius = self._computeNormalizedRadius(distance2axis) 
        thrustFactor = thrustBC * cellVolume / totalVolume

        # add thrust
        force[:, :] = np.outer(
                self._thrustSpline(normalizedRadius) * thrustFactor, 
                self.thrustVector
                )

        # add swirl (tangential contribution)
        force += np.multiply(
                tangent.T, 
                np.array([self._tangentSpline(normalizedRadius) * thrustFactor])
                ).T

        # add radial contribution
        radius = np.cross(tangent, self.thrustVector)
        force += np.multiply(
                radius.T, 
                np.array([self._radialSpline(normalizedRadius) * thrustFactor])
                ).T

        return force

    def computeCellHeatVector(self, distance2axis, distance2plane, tangent, cellVolume, totalVolume):
        heat = np.zeros_like(distance2axis)
        heatBC = self._boundaryConditions[self._heatBcName]

        normalizedRadius = self._computeNormalizedRadius(distance2axis) 

        heat = self._heatSpline(normalizedRadius) * heatBC * cellVolume / totalVolume

        return heat

