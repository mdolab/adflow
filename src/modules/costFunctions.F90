!
!     ******************************************************************
!     *                                                                *
!     * File:          costFunctions.F90                               *
!     * Author:        Gaetan Kenway                                   *
!     * Starting date: 11-07-2010                                      *
!     * Last modified: 11-07-2010                                      *
!     *                                                                *
!     ******************************************************************
!
module costFunctions
  !
  !     ******************************************************************
  !     *                                                                *
  !     * This module contains the paramter values for the solver and    *
  !     * adjoint cost functions                                         *
  !     *                                                                *
  !     ******************************************************************
  !
  use constants
  implicit none
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Cost functions.                                                *
  !     *                                                                *
  !     ******************************************************************

  integer(kind=intType), parameter :: nCostFunction = 37
  integer(kind=intType), parameter :: costFuncLift       = 1,&
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
       costFuncBendingCoef= 37

  real(kind=realType), allocatable, dimension(:)   :: functionValue
  real(kind=realType), dimension(6, nCostFunction) :: costFuncMat
#ifndef USE_TAPENADE
  REAL(kind=realtype), DIMENSION(6, ncostfunction) :: costfuncmatd
  real(kind=realType), dimension(:,:,:), allocatable :: dCostFuncmatdExtra
  real(kind=realType), dimension(6) :: FMExtra
#endif
end module costFunctions


