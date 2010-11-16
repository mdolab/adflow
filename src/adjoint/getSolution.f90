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
  real(kind=realType) :: cl0,cd0,cmz0,dcldalpha,dcddalpha,dcmzdalpha
  real(kind=realType) :: dcmzdalphadot,dcmzdq
  real(kind=realType),dimension(:),allocatable :: localVal,globalVal
  ! Function values

  if (.not. allocated(functionValue)) then
     allocate(functionValue(nCostFunction))
  end if

  functionValue(:) = 0.0
  call computeAeroCoef(sps)
  
  if(TSStability)then

     call computeTSDerivatives(cl0,cd0,cmz0,dcldalpha,dcddalpha,&
          dcmzdalpha,dcmzdalphadot,dcmzdq)
     functionValue(costFuncCmzAlpha)     = dcmzdalpha
     functionValue( costFuncCm0)         = cmz0
     functionValue( costFuncClAlpha)     = dcldalpha
     functionValue( costFuncCl0  )       = cl0
     functionValue( costFuncCdAlpha )    = dcmzdalpha
     functionValue( costFuncCd0 )        = cd0
     functionValue( costFuncCmzAlphaDot) = dcmzdalphadot
     functionValue( costFuncCmzq)         = dcmzdq
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
