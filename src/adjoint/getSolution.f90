!
!     ******************************************************************
!     *                                                                *
!     * File:          getSolution.f90                                *
!     * Authors:       C.A(Sandy) Mader                                *
!     * Starting date: 23-07-2008                                      *
!     * Last modified: 23-07-2008                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine getSolution(sps)

  use costFunctions
  use inputTSStabDeriv !TSStability
  use communication
  implicit none

  integer(kind=intType) :: sps,ierr
  !
  !     Local variables.
  !
  real(kind=realType)   :: alpha, beta

  real(kind=realType) :: cl0,cd0,cmz0
  real(kind=realType) :: dcldalpha,dcddalpha,dcmzdalpha
  real(kind=realType) :: dcldalphaDot,dcddalphaDot,dcmzdalphaDot
  real(kind=realType) :: dcldq,dcddq,dcmzdq
  real(kind=realType) :: dcldqdot,dcddqdot,dcmzdqdot

  real(kind=realType),dimension(:),allocatable :: localVal,globalVal
  ! Function values

  if (.not. allocated(functionValue)) then
     allocate(functionValue(nCostFunction))
  end if

  functionValue(:) = 0.0
  call computeAeroCoef(sps)
  
  if(TSStability)then
     call computeTSDerivatives(cl0,cd0,cmz0,dcldalpha,dcddalpha,&
          dcmzdalpha,dcldalphadot,dcddalphadot,dcmzdalphadot,dcldq,&
          dcddq,dcmzdq,dcldqdot,dcddqdot,dcmzdqdot)

     functionValue( costFuncCl0  )       = cl0
     functionValue( costFuncCd0 )        = cd0
     functionValue( costFuncCm0 )        = cmz0

     functionValue( costFuncClAlpha)     = dcldalpha
     functionValue( costFuncCdAlpha)     = dcddalpha
     functionValue( costFuncCmzAlpha)    = dcmzdalpha

     functionValue( costFuncClAlphaDot)     = dcldalphadot
     functionValue( costFuncCdAlphaDot)     = dcddalphadot
     functionValue( costFuncCmzAlphaDot)    = dcmzdalphadot
     
     functionValue( costFuncClq)         = dcldq
     functionValue( costFuncCdq)         = dcddq
     functionValue( costFuncCmzq)        = dcmzdq

     functionValue( costFuncClqDot)         = dcldqDot
     functionValue( costFuncCdqDot)         = dcddqDot
     functionValue( costFuncCmzqDot)        = dcmzdqDot

  end if

  ! Now we will mpi_allReduce them
  allocate(localVal(nCostFunction),globalVal(nCostFunction))
  ! Copy Everything into the list

  localVal(:) = functionValue(:)
  call mpi_allreduce(localVal, globalVal, nCostFunction, sumb_real, &
       mpi_sum, SUmb_comm_world, ierr)
 
  functionValue(:) = globalVal(:)

  deallocate(localVal,globalVal)
  
end subroutine getSolution
