from abc import ABC, abstractmethod
from typing import Tuple

import numpy as np
import numpy.typing as npt


class AbstractActuatorRegion(ABC):
    def __init__(self, centerPoint: npt.NDArray, thrustVector: npt.NDArray):
        if centerPoint.shape != (3,):
            raise ValueError('"centerPoint" must have shape "(3,)"') 
        if thrustVector.shape != (3,):
            raise ValueError('"thrustVector" must have shape "(3,)"') 
        if np.linalg.norm(thrustVector) == 0:
            raise ValueError('"trustVector" can not have a length of "0"')

        self.centerPoint = centerPoint
        self.thrustVector = thrustVector / np.linalg.norm(thrustVector)
        self.iRegion = -1
        self.nLocalCells = -1
        self.familyName = None
        self._boundaryConditions = dict()

    def setBCValue(self, bcName, bcValue):
        if bcName not in self._boundaryConditions.keys():
            raise ValueError(f'"{bcName}" is not a valid boundary Condition. Valid Conditions are: {list(self._boundaryConditions.keys())}')
        
        self._boundaryConditions[bcName] = bcValue

    @abstractmethod
    def tagActiveCells(self, distance2axis: npt.NDArray, distance2plane: npt.NDArray) -> npt.NDArray:
        flags = np.zeros_like(distance2axis)
        return flags

    @abstractmethod
    def computeCellForceVector(self, distance2axis: npt.NDArray, distance2plane: npt.NDArray, cellVolume: npt.NDArray, totalVolume: float) -> npt.NDArray:
        force = np.zeros((3, len(distance2axis)))
        return force

    @abstractmethod
    def computeCellHeatVector(self, distance2axis: npt.NDArray, distance2plane: npt.NDArray, cellVolume: npt.NDArray, totalVolume: float) -> npt.NDArray:
        heat = np.zeros_like(distance2axis)
        return heat





class UniformActuatorRegion(AbstractActuatorRegion):
    def __init__(self, centerPoint: npt.NDArray, thrustVector: npt.NDArray, outerDiameter: float, regionDepth: float, thrust: float = 0., heat: float = 0.):
        super().__init__(centerPoint, thrustVector)

        if regionDepth <= 0:
            raise ValueError('"depth" must be greater than 0.')
        if outerDiameter <= 0:
            raise ValueError('"outerDiameter" must be greater than 0.')
        if thrust < 0:
            raise ValueError('"thrust" must not be smaller than 0.')
        if heat < 0:
            raise ValueError('"heat" must not be smaller than 0.')

        self._outerDiameter = outerDiameter
        self._regionDepth = regionDepth

        self._heatBcName = 'Heat'
        self._thrustBcName = 'Thrust'

        self._boundaryConditions = {
                self._thrustBcName: thrust,
                self._heatBcName: heat,
                }

    def tagActiveCells(self, distance2axis, distance2plane):
        flags = np.zeros_like(distance2axis)

        indices =  np.logical_and(
                distance2axis <= self._outerDiameter/2,
                np.logical_and(
                    distance2plane >= 0,
                    distance2plane <= self._regionDepth
                    )
                )

        flags[indices] = 1

        return flags

    def computeCellForceVector(self, distance2axis: npt.NDArray, distance2plane: npt.NDArray, cellVolume: npt.NDArray, totalVolume: float):
        force = np.zeros((len(distance2axis), 3))
        thrustBC = self._boundaryConditions[self._thrustBcName]

        force[:, :] = np.outer(
                thrustBC * cellVolume / totalVolume, 
                self.thrustVector
                )

        return force

    def computeCellHeatVector(self, distance2axis: npt.NDArray, distance2plane: npt.NDArray, cellVolume: npt.NDArray, totalVolume: float):
        heat = np.zeros_like(distance2axis)
        heatBC = self._boundaryConditions[self._heatBcName]

        heat = heatBC * cellVolume / totalVolume

        return heat

