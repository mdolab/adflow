!
!     ******************************************************************
!     *                                                                *
!     * File:          verifyBendingDerivatives.f90                    *
!     * Authors:       C.A(Sandy) Mader                                *
!     * Starting date: 16-07-2011                                      *
!     * Last modified: 16-07-2011                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine verifyBendingDerivatives()

  use costFunctions
  use inputPhysics
  use inputTimeSpectral
  use communication
  implicit none

  integer(kind=intType) :: sps,ierr,i
  !
  !     Local variables.
  !
  real(kind=realType), dimension(nCostFunction)::globalCFVals,globalCFValsRef,GlobalCFValsb
  real(kind=realType)::bendingMoment,bendingMomentRef,bendingMomentb
  real(kind=realType)::deltax = 1e-6
  real(kind=realType),dimension(3)::pointrefRef
  ! Function values

  do sps = 1,nTimeIntervalsSpectral
     call computeAeroCoef(globalCFVals,sps)
     
     call computeRootBendingMoment(globalCFVals,bendingMoment)
     globalCFValsRef = GlobalCFVals
     bendingMomentRef = bendingMoment
     if(myid==0)then
        print *,'Bending Coefficient',bendingMoment,sps
     end if
     pointrefb =0.0
     bendingmomentb = 1.0
     call COMPUTEROOTBENDINGMOMENT_B(globalCFVals, globalCFValsb, bendingmoment, &
          &  bendingmomentb)
     if(myid==0)then
        print *,'AD derivatives',globalCFValsb,pointrefb
     endif
     
     do i =1,nCostfunction
        globalCFVals = globalCFValsRef
        globalCFVals(i) = globalCFValsRef(i)+deltax
        call computerootBendingMoment(globalCFVals,bendingMoment)
        if(myid==0)then
           print *,'derivative',i,(bendingMoment-bendingMomentRef)/deltax
        endif
     end do
     pointrefref = pointref
     do i =1,3
        pointref = pointrefref
        pointref(i) = pointrefRef(i)+deltax
        call computerootBendingMoment(globalCFVals,bendingMoment)
        if(myid==0)then
           print *,'Pointref derivative',i,(bendingMoment-bendingMomentRef)/deltax
        endif
     end do
  end do

end subroutine verifyBendingDerivatives
