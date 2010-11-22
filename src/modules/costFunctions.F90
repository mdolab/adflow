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

  integer(kind=intType), parameter :: nCostFunction = 31_intType

  integer(kind=intType), parameter :: costFuncLift       = 1_intType,&
       costFuncDrag       = 2_intType,&
       costFuncLiftCoef   = 3_intType,&
       costFuncDragCoef   = 4_intType,&
       costFuncForceX     = 5_intType,&
       costFuncForceY     = 6_intType,&
       costFuncForceZ     = 7_intType,&
       costFuncForceXCoef = 8_intType,&
       costFuncForceYCoef = 9_intType,&
       costFuncForceZCoef = 10_intType,&
       costFuncMomX       = 11_intType,&
       costFuncMomY       = 12_intType,&
       costFuncMomZ       = 13_intType,&
       costFuncMomXCoef   = 14_intType,&
       costFuncMomYCoef   = 15_intType,&
       costFuncMomZCoef   = 16_intType,&

       costFuncCm0        = 17_intType,&
       costFuncCmzAlpha   = 18_intType,&
       costFuncCmzAlphaDot= 19_intType,&
       costFuncCmzq       = 20_intType,&
       costFuncCmzqDot    = 21_intType,&

       costFuncCl0        = 22_intType,&
       costFuncClAlpha    = 23_intType,&
       costFuncClAlphaDot = 24_intType,&
       costFuncClq        = 25_intType,&
       costFuncClqDot     = 26_intType,&

       costFuncCd0        = 27_intType,&
       costFuncCdAlpha    = 28_intType,& 
       costFuncCdAlphadot = 29_intType,&
       costFuncCdq        = 30_intType,&
       costFuncCdqDot     = 31_intType

  real(kind=realType), allocatable, dimension(:)   :: functionValue
end module costFunctions


