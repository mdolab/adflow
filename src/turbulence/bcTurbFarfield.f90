subroutine bcTurbFarfield(nn)
  !
  !       bcTurbFarfield applies the implicit treatment of the           
  !       farfield boundary condition to subface nn. As the farfield     
  !       boundary condition is independent of the turbulence model,     
  !       this routine is valid for all models. It is assumed that the   
  !       pointers in blockPointers are already set to the correct       
  !       block on the correct grid level.                               
  !
  use constants
  use blockPointers
  use flowVarRefState
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: nn
  !
  !      Local variables.
  !
  integer(kind=intType) :: i, j, l

  real(kind=realType) :: nnx, nny, nnz, dot

  ! Loop over the faces of the subfaces and set the values of
  ! bmt and bvt for an implicit treatment.

  do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
     do i=BCData(nn)%icBeg, BCData(nn)%icEnd

        ! Determine the dot product between the outward pointing
        ! normal and the free stream velocity direction and add the
        ! possible grid velocity.

        dot = BCData(nn)%norm(i,j,1)*wInf(ivx) + &
             BCData(nn)%norm(i,j,2)*wInf(ivy) + &
             BCData(nn)%norm(i,j,3)*wInf(ivz) - BCData(nn)%rface(i,j)

        ! Determine whether we are dealing with an inflow or
        ! outflow boundary here.

        if(dot > zero) then

           ! Outflow. Simply extrapolation or zero Neumann BC
           ! of the turbulent variables.

           do l=nt1,nt2
              select case (BCFaceID(nn))
              case (iMin)
                 bmti1(i,j,l,l) = -one
              case (iMax)
                 bmti2(i,j,l,l) = -one
              case (jMin)
                 bmtj1(i,j,l,l) = -one
              case (jMax)
                 bmtj2(i,j,l,l) = -one
              case (kMin)
                 bmtk1(i,j,l,l) = -one
              case (kMax)
                 bmtk2(i,j,l,l) = -one
              end select
           end do

        else

           ! Inflow. Turbulent variables are prescribed.

           do l=nt1,nt2
              select case(BCFaceID(nn))
              case (iMin)
                 bvti1(i,j,l) = wInf(l)
              case (iMax)
                 bvti2(i,j,l) = wInf(l)
              case (jMin)
                 bvtj1(i,j,l) = wInf(l)
              case (jMax)
                 bvtj2(i,j,l) = wInf(l)
              case (kMin)
                 bvtk1(i,j,l) = wInf(l)
              case (kMax)
                 bvtk2(i,j,l) = wInf(l)
              end select
           enddo
        endif
     enddo
  enddo
end subroutine bcTurbFarfield
