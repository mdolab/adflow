!        generated by tapenade     (inria, tropics team)
!  tapenade 3.10 (r5363) -  9 sep 2014 09:53
!
!  differentiation of bcturbfarfield in forward (tangent) mode (with options i4 dr8 r8):
!   variations   of useful results: *bvtj1 *bvtj2 *bvtk1 *bvtk2
!                *bvti1 *bvti2
!   with respect to varying inputs: winf *bvtj1 *bvtj2 *bvtk1 *bvtk2
!                *bvti1 *bvti2
!   plus diff mem management of: bvtj1:in bvtj2:in bvtk1:in bvtk2:in
!                bvti1:in bvti2:in bcdata:in
!
!      ******************************************************************
!      *                                                                *
!      * file:          bcturbfarfield.f90                              *
!      * author:        georgi kalitzin, edwin van der weide            *
!      * starting date: 06-15-2003                                      *
!      * last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine bcturbfarfield_d(nn)
!
!      ******************************************************************
!      *                                                                *
!      * bcturbfarfield applies the implicit treatment of the           *
!      * farfield boundary condition to subface nn. as the farfield     *
!      * boundary condition is independent of the turbulence model,     *
!      * this routine is valid for all models. it is assumed that the   *
!      * pointers in blockpointers are already set to the correct       *
!      * block on the correct grid level.                               *
!      *                                                                *
!      ******************************************************************
!
  use constants
  use blockpointers
  use flowvarrefstate
  implicit none
!
!      subroutine arguments.
!
  integer(kind=inttype), intent(in) :: nn
!
!      local variables.
!
  integer(kind=inttype) :: i, j, l
  real(kind=realtype) :: nnx, nny, nnz, dot
!
!      ******************************************************************
!      *                                                                *
!      * begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
! loop over the faces of the subfaces and set the values of
! bmt and bvt for an implicit treatment.
  do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
    do i=bcdata(nn)%icbeg,bcdata(nn)%icend
! determine the dot product between the outward pointing
! normal and the free stream velocity direction and add the
! possible grid velocity.
      dot = bcdata(nn)%norm(i, j, 1)*winf(ivx) + bcdata(nn)%norm(i, j, 2&
&       )*winf(ivy) + bcdata(nn)%norm(i, j, 3)*winf(ivz) - bcdata(nn)%&
&       rface(i, j)
! determine whether we are dealing with an inflow or
! outflow boundary here.
      if (dot .gt. zero) then
! outflow. simply extrapolation or zero neumann bc
! of the turbulent variables.
        do l=nt1,nt2
          select case  (bcfaceid(nn)) 
          case (imin) 
            bmti1(i, j, l, l) = -one
          case (imax) 
            bmti2(i, j, l, l) = -one
          case (jmin) 
            bmtj1(i, j, l, l) = -one
          case (jmax) 
            bmtj2(i, j, l, l) = -one
          case (kmin) 
            bmtk1(i, j, l, l) = -one
          case (kmax) 
            bmtk2(i, j, l, l) = -one
          end select
        end do
      else
! inflow. turbulent variables are prescribed.
        do l=nt1,nt2
          select case  (bcfaceid(nn)) 
          case (imin) 
            bvti1d(i, j, l) = winfd(l)
            bvti1(i, j, l) = winf(l)
          case (imax) 
            bvti2d(i, j, l) = winfd(l)
            bvti2(i, j, l) = winf(l)
          case (jmin) 
            bvtj1d(i, j, l) = winfd(l)
            bvtj1(i, j, l) = winf(l)
          case (jmax) 
            bvtj2d(i, j, l) = winfd(l)
            bvtj2(i, j, l) = winf(l)
          case (kmin) 
            bvtk1d(i, j, l) = winfd(l)
            bvtk1(i, j, l) = winf(l)
          case (kmax) 
            bvtk2d(i, j, l) = winfd(l)
            bvtk2(i, j, l) = winf(l)
          end select
        end do
      end if
    end do
  end do
end subroutine bcturbfarfield_d
