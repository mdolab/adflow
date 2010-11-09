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

  integer(kind=intType), parameter :: nCostFunction = 24_intType

  integer(kind=intType), parameter :: costFuncLift = 1_intType,&
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
       costFuncCmzAlpha   = 17_intType,&
       costFuncCm0        = 18_intType,&
       costFuncClAlpha    = 19_intType,&
       costFuncCl0        = 20_intType,&
       costFuncCdAlpha    = 21_intType,&
       costFuncCd0        = 22_intType,&
       costFuncCmzAlphaDot= 23_intType,&
       costFuncCmzq       = 24_intType

  real(kind=realType), allocatable, dimension(:)   :: functionValue
end module costFunctions
