from abc import ABC, abstractmethod
from typing import Tuple

import numpy as np
import numpy.typing as npt


class AbstractActuatorRegion(ABC):
    def __init__(self, centerPoint: npt.NDArray, thrustVector: npt.NDArray):
        # todo add error checking

        self.centerPoint = centerPoint
        self.thrustVector = thrustVector
        self.iRegion = -1
        self.nLocalCells = -1

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
    def __init__(self, centerPoint: npt.NDArray, thrustVector: npt.NDArray, outerDiameter: float, regionDepth: float, thrust: float, heat: float):
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
        self._thrust = thrust
        self._heat = heat


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

        force[:, :] = np.outer(
                self._thrust * cellVolume / totalVolume, 
                self.thrustVector
                )

        return force

    def computeCellHeatVector(self, distance2axis: npt.NDArray, distance2plane: npt.NDArray, cellVolume: npt.NDArray, totalVolume: float):
        heat = np.zeros_like(distance2axis)

        heat = self._heat * cellVolume / totalVolume

        return heat

