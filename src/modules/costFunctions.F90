!
!      File:          costFunctions.F90                               
!      Author:        Gaetan Kenway                                   
!      Starting date: 11-07-2010                                      
!      Last modified: 11-07-2010                                      
!
module costFunctions
  !
  !      This module contains the paramter values for the solver and    
  !      adjoint cost functions                                         
  !
  use constants, only : intType, realType, maxCGNSNameLen
  implicit none
  !
  !      Cost functions.                                                

  integer(kind=intType), parameter :: nCostFunction = 42
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
       costFuncCavitation = 42

  integer(kind=intType), parameter :: &
       icFp1 = 1, &
       icFp2 = 3, &
       icFv1 = 4, &
       icFv2 = 6, &
       icMp1 = 7, &
       icMp2 = 9, &
       icMv1 =10, &
       icMv2 =12, &
       iSepSensor = 13, &
       iSepAvg1 = 14, &
       iSepAvg2 = 16, &
       iFp1 =17, &
       iFp2 =19, &
       iFv1 =20, &
       iFv2 =22, &
       iMp1 =23, &
       iMp2 =25, &
       iMv1 =26, &
       iMv2 =28, &
       iCavitation=29, &
       iyPlus = 30
  integer(kind=intType), parameter :: nLocalValues=30
  real(kind=realType), dimension(nCostFunction) ::  funcValues
#ifndef USE_TAPENADE 
  real(kind=realType), dimension(nCostFunction) ::  funcValuesd
#endif
  
  real(kind=realtype) :: sepSensorOffset, sepSensorSharpness
  character(len=maxCGNSNameLen), dimension(:), allocatable :: maskFams
end module costFunctions


