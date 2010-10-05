

subroutine getSurfaceNormalsCouplingAdj(pts,normAdj,rightHanded)

  
  use constants
  implicit none

  ! Subroutine Arguments
  
  real(kind=realType), intent(in) :: pts(3,3,3)
  logical            , intent(in) :: rightHanded
  real(kind=realType), intent(out):: normAdj(3,2,2)

  ! Local Variables

  integer(kind=intType) :: i,j
  real(kind=realType) :: v1(3),v2(3),fact2

  if (rightHanded) then
     fact2 = half
  else
     fact2 = -half
  end if

  do j=1,2
     do i=1,2

        v1(:) = pts(:,i+1,j+1)-pts(:,i,j)
        v2(:) = pts(:,i,j+1)-pts(:,i+1,j)

        ! The face normal, which is the cross product of the two
        ! diagonal vectors times fact; remember that fact2 is
        ! either -0.5 or 0.5.

        normAdj(1,i,j) = fact2*(v1(2)*v2(3) - v1(3)*v2(2))
        normAdj(2,i,j) = fact2*(v1(3)*v2(1) - v1(1)*v2(3))
        normAdj(3,i,j) = fact2*(v1(1)*v2(2) - v1(2)*v2(1))
     end do
  end do


end subroutine getSurfaceNormalsCouplingAdj
