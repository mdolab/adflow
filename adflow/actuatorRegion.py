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

