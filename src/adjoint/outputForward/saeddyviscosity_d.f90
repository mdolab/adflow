!        generated by tapenade     (inria, tropics team)
!  tapenade 3.10 (r5363) -  9 sep 2014 09:53
!
!  differentiation of saeddyviscosity in forward (tangent) mode (with options i4 dr8 r8):
!   variations   of useful results: *rev
!   with respect to varying inputs: *w *rlv
!   plus diff mem management of: rev:in w:in rlv:in
!      ==================================================================
!      ==================================================================
subroutine saeddyviscosity_d()
!
!      ******************************************************************
!      *                                                                *
!      * saeddyviscosity computes the eddy-viscosity according to the   *
!      * spalart-allmaras model for the block given in blockpointers.   *
!      * this routine for both the original version as well as the      *
!      * modified version according to edwards.                         *
!      *                                                                *
!      ******************************************************************
!
  use constants
  use blockpointers
  use constants
  use paramturb
  implicit none
!
!      local variables.
!
  integer(kind=inttype) :: i, j, k, ii
  real(kind=realtype) :: chi, chi3, fv1, rnusa, cv13
  real(kind=realtype) :: chid, chi3d, fv1d, rnusad
!
!      ******************************************************************
!      *                                                                *
!      * begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
! store the cv1^3; cv1 is a constant of the spalart-allmaras model.
  cv13 = rsacv1**3
  revd = 0.0_8
! loop over the cells of this block and compute the eddy viscosity.
! do not include halo's.
  do k=1,ke
    do j=1,je
      do i=1,ie
        rnusad = wd(i, j, k, itu1)*w(i, j, k, irho) + w(i, j, k, itu1)*&
&         wd(i, j, k, irho)
        rnusa = w(i, j, k, itu1)*w(i, j, k, irho)
        chid = (rnusad*rlv(i, j, k)-rnusa*rlvd(i, j, k))/rlv(i, j, k)**2
        chi = rnusa/rlv(i, j, k)
        chi3d = 3*chi**2*chid
        chi3 = chi**3
        fv1d = (chi3d*(chi3+cv13)-chi3*chi3d)/(chi3+cv13)**2
        fv1 = chi3/(chi3+cv13)
        revd(i, j, k) = fv1d*rnusa + fv1*rnusad
        rev(i, j, k) = fv1*rnusa
      end do
    end do
  end do
end subroutine saeddyviscosity_d
