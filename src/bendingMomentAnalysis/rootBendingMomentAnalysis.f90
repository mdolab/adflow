!
!****************************************************
!
!   Filename: rootBendingMomentAnalysis.f90
!   Author: C.A.(Sandy) Mader
!   Date Started: July 10, 2011
!   Date Modified: July 10,2011
!
!****************************************************

subroutine computeRootBendingMoment(cf, cm, liftIndex, bendingMoment)

  !*******************************************************
  !                                                      *
  ! Compute a normalized bending moment coefficient from *
  ! the force and moment coefficient. At the moment this *
  ! Routine only works for a half body. Additional logic *
  ! would be needed for a full body.                     *
  !                                                      *
  !*******************************************************

  use inputPhysics  
  use costFunctions
  use communication
  use inputIteration
  implicit none

  !input/output variables
  real(kind=realType), intent(in), dimension(3) :: cf, cm
  integer(kind=intType), intent(in) :: liftIndex
  real(kind=realType), intent(out) :: bendingMoment

  !Subroutine Variables
  real(kind=realType):: elasticMomentx, elasticMomenty, elasticMomentz
  bendingMoment = zero
  if (liftIndex == 2) then
     !z out wing sum momentx,momentz
     elasticMomentx = cm(1) + cf(2)*(pointRefEC(3)-pointRef(3))/lengthref-cf(3)*(pointRefEC(2)-pointRef(2))/lengthref
     elasticMomentz = cm(3) - cf(2)*(pointRefEC(1)-pointref(1))/lengthref+cf(1)*(pointRefEC(2)-pointRef(2))/lengthref
     bendingMoment = sqrt(elasticMomentx**2+elasticMomentz**2)
  elseif (liftIndex == 3) then
     !y out wing sum momentx,momenty
     elasticMomentx = cm(1) + cf(3)*(pointrefEC(2)-pointRef(2))/lengthref+cf(3)*(pointrefEC(3)-pointref(3))/lengthref
     elasticMomenty = cm(2) + cf(3)*(pointRefEC(1)-pointRef(1))/lengthref+cf(1)*(pointrefEC(3)-pointRef(3))/lengthref
     bendingMoment = sqrt(elasticMomentx**2+elasticMomenty**2)
  end if

end subroutine computeRootBendingMoment

