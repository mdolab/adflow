module costFunctions
  !
  !      This module contains the paramter values for the solver and
  !      adjoint cost functions
  !
  use constants, only : intType, realType, maxCGNSNameLen
  implicit none
  !
  !      Cost functions.

  integer(kind=intType), parameter :: nCostFunction = 47
  integer(kind=intType), parameter :: &
       costFuncLift       = 1,&
       costFuncDrag       = 2,&
       costFuncLiftCoef   = 3,&
       costFuncDragCoef   = 4,&
       costFuncForceX     = 5,&
       costFuncForceY     = 6,&
       costFuncForceZ     = 7,&
       costFuncForceXCoef = 8,&
       costFuncForceYCoef = 9,&
       costFuncForceZCoef = 10,&
       costFuncMomX       = 11,&
       costFuncMomY       = 12,&
       costFuncMomZ       = 13,&
       costFuncMomXCoef   = 14,&
       costFuncMomYCoef   = 15,&
       costFuncMomZCoef   = 16,&
       costFuncCm0        = 17,&
       costFuncCmzAlpha   = 18,&
       costFuncCmzAlphaDot= 19,&
       costFuncCmzq       = 20,&
       costFuncCmzqDot    = 21,&
       costFuncCl0        = 22,&
       costFuncClAlpha    = 23,&
       costFuncClAlphaDot = 24,&
       costFuncClq        = 25,&
       costFuncClqDot     = 26,&
       costFuncCd0        = 27,&
       costFuncCdAlpha    = 28,&
       costFuncCdAlphadot = 29,&
       costFuncCdq        = 30,&
       costFuncCdqDot     = 31,&
       costFuncCfy0       = 32,&
       costFuncCfyAlpha   = 33,&
       costFuncCfyAlphadot= 34,&
       costFuncCfyq       = 35,&
       costFuncCfyqDot    = 36,&
       costFuncBendingCoef= 37,&
       costFuncSepSensor  = 38,&
       costFuncSepSensorAvgX = 39, &
       costFuncSepSensorAvgY = 40, &
       costFuncSepSensorAvgZ = 41, &
       costFuncCavitation = 42, &
       costFuncMdot = 43, &
       costFuncMavgPtot = 44, &
       costFuncMavgTtot = 45, &
       costFuncMavgPs = 46, & 
       costFuncMavgMN = 47

  integer(kind=intType), parameter :: &
       iFp =  1, &
       iFv =  4, &
       iMp =  7, &
       iMv = 10, &
       iSepSensor  = 13, &
       iSepAvg     = 14, &
       iCavitation = 17, &
       iyPlus    = 18, &
       iMassFlow = 19, &
       iMassPTot = 20, &
       iMassTtot = 21, &
       iMassPs   = 22, &
       iFlowFp   = 23, & 
       iFlowMp   = 26, & 
       iFlowFm   = 29, & 
       iFlowMm   = 32, & 
       iMassMN   = 35 


  integer(kind=intType), parameter :: nLocalValues=36

  real(kind=realType), dimension(nCostFunction) ::  funcValues
#ifndef USE_TAPENADE
  real(kind=realType), dimension(nCostFunction) ::  funcValuesd
#endif

  real(kind=realtype) :: sepSensorOffset, sepSensorSharpness
  character(len=maxCGNSNameLen), dimension(:), allocatable :: maskFams
end module costFunctions
