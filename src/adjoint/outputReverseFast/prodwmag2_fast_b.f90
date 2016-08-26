!        generated by tapenade     (inria, tropics team)
!  tapenade 3.10 (r5363) -  9 sep 2014 09:53
!
!  differentiation of prodwmag2 in reverse (adjoint) mode (with options i4 dr8 r8 noisize):
!   gradient     of useful results: *w *scratch
!   with respect to varying inputs: *w *scratch
!   rw status of diff variables: *w:incr *scratch:in-out
!   plus diff mem management of: w:in scratch:in
!
!       file:          prodwmag2.f90                                   
!       author:        georgi kalitzin, edwin van der weide            
!       starting date: 06-23-2003                                      
!       last modified: 06-12-2005                                      
!
subroutine prodwmag2_fast_b()
!
!       prodwmag2 computes the term:                                   
!          2*oij*oij  with oij=0.5*(duidxj - dujdxi).                  
!       this is equal to the magnitude squared of the vorticity.       
!       it is assumed that the pointer vort, stored in turbmod, is     
!       already set to the correct entry.                              
!
  use constants
  use blockpointers
  use flowvarrefstate
  use section
  use turbmod
  implicit none
!
!      local variables.
!
  integer :: i, j, k, ii
  real(kind=realtype) :: uuy, uuz, vvx, vvz, wwx, wwy
  real(kind=realtype) :: uuyd, uuzd, vvxd, vvzd, wwxd, wwyd
  real(kind=realtype) :: fact, vortx, vorty, vortz
  real(kind=realtype) :: vortxd, vortyd, vortzd
  real(kind=realtype) :: omegax, omegay, omegaz
  intrinsic mod
!
!       begin execution                                                
!
! determine the non-dimensional wheel speed of this block.
  omegax = timeref*sections(sectionid)%rotrate(1)
  omegay = timeref*sections(sectionid)%rotrate(2)
  omegaz = timeref*sections(sectionid)%rotrate(3)
  do ii=0,nx*ny*nz-1
    i = mod(ii, nx) + 2
    j = mod(ii/nx, ny) + 2
    k = ii/(nx*ny) + 2
! compute the necessary derivatives of u in the cell center.
! use is made of the fact that the surrounding normals sum up
! to zero, such that the cell i,j,k does not give a
! contribution. the gradient is scaled by a factor 2*vol.
    uuy = w(i+1, j, k, ivx)*si(i, j, k, 2) - w(i-1, j, k, ivx)*si(i-1, j&
&     , k, 2) + w(i, j+1, k, ivx)*sj(i, j, k, 2) - w(i, j-1, k, ivx)*sj(&
&     i, j-1, k, 2) + w(i, j, k+1, ivx)*sk(i, j, k, 2) - w(i, j, k-1, &
&     ivx)*sk(i, j, k-1, 2)
    uuz = w(i+1, j, k, ivx)*si(i, j, k, 3) - w(i-1, j, k, ivx)*si(i-1, j&
&     , k, 3) + w(i, j+1, k, ivx)*sj(i, j, k, 3) - w(i, j-1, k, ivx)*sj(&
&     i, j-1, k, 3) + w(i, j, k+1, ivx)*sk(i, j, k, 3) - w(i, j, k-1, &
&     ivx)*sk(i, j, k-1, 3)
! idem for the gradient of v.
    vvx = w(i+1, j, k, ivy)*si(i, j, k, 1) - w(i-1, j, k, ivy)*si(i-1, j&
&     , k, 1) + w(i, j+1, k, ivy)*sj(i, j, k, 1) - w(i, j-1, k, ivy)*sj(&
&     i, j-1, k, 1) + w(i, j, k+1, ivy)*sk(i, j, k, 1) - w(i, j, k-1, &
&     ivy)*sk(i, j, k-1, 1)
    vvz = w(i+1, j, k, ivy)*si(i, j, k, 3) - w(i-1, j, k, ivy)*si(i-1, j&
&     , k, 3) + w(i, j+1, k, ivy)*sj(i, j, k, 3) - w(i, j-1, k, ivy)*sj(&
&     i, j-1, k, 3) + w(i, j, k+1, ivy)*sk(i, j, k, 3) - w(i, j, k-1, &
&     ivy)*sk(i, j, k-1, 3)
! and for the gradient of w.
    wwx = w(i+1, j, k, ivz)*si(i, j, k, 1) - w(i-1, j, k, ivz)*si(i-1, j&
&     , k, 1) + w(i, j+1, k, ivz)*sj(i, j, k, 1) - w(i, j-1, k, ivz)*sj(&
&     i, j-1, k, 1) + w(i, j, k+1, ivz)*sk(i, j, k, 1) - w(i, j, k-1, &
&     ivz)*sk(i, j, k-1, 1)
    wwy = w(i+1, j, k, ivz)*si(i, j, k, 2) - w(i-1, j, k, ivz)*si(i-1, j&
&     , k, 2) + w(i, j+1, k, ivz)*sj(i, j, k, 2) - w(i, j-1, k, ivz)*sj(&
&     i, j-1, k, 2) + w(i, j, k+1, ivz)*sk(i, j, k, 2) - w(i, j, k-1, &
&     ivz)*sk(i, j, k-1, 2)
! compute the three components of the vorticity vector.
! substract the part coming from the rotating frame.
    fact = half/vol(i, j, k)
    vortx = fact*(wwy-vvz) - two*omegax
    vorty = fact*(uuz-wwx) - two*omegay
    vortz = fact*(vvx-uuy) - two*omegaz
! compute the magnitude squared of the vorticity.
    vortxd = 2*vortx*scratchd(i, j, k, ivort)
    vortyd = 2*vorty*scratchd(i, j, k, ivort)
    vortzd = 2*vortz*scratchd(i, j, k, ivort)
    scratchd(i, j, k, ivort) = 0.0_8
    vvxd = fact*vortzd
    uuyd = -(fact*vortzd)
    uuzd = fact*vortyd
    wwxd = -(fact*vortyd)
    wwyd = fact*vortxd
    vvzd = -(fact*vortxd)
    wd(i+1, j, k, ivz) = wd(i+1, j, k, ivz) + si(i, j, k, 2)*wwyd
    wd(i-1, j, k, ivz) = wd(i-1, j, k, ivz) - si(i-1, j, k, 2)*wwyd
    wd(i, j+1, k, ivz) = wd(i, j+1, k, ivz) + sj(i, j, k, 2)*wwyd
    wd(i, j, k+1, ivz) = wd(i, j, k+1, ivz) + sk(i, j, k, 2)*wwyd
    wd(i, j-1, k, ivz) = wd(i, j-1, k, ivz) - sj(i, j-1, k, 2)*wwyd
    wd(i, j, k-1, ivz) = wd(i, j, k-1, ivz) - sk(i, j, k-1, 2)*wwyd
    wd(i+1, j, k, ivz) = wd(i+1, j, k, ivz) + si(i, j, k, 1)*wwxd
    wd(i-1, j, k, ivz) = wd(i-1, j, k, ivz) - si(i-1, j, k, 1)*wwxd
    wd(i, j+1, k, ivz) = wd(i, j+1, k, ivz) + sj(i, j, k, 1)*wwxd
    wd(i, j, k+1, ivz) = wd(i, j, k+1, ivz) + sk(i, j, k, 1)*wwxd
    wd(i, j-1, k, ivz) = wd(i, j-1, k, ivz) - sj(i, j-1, k, 1)*wwxd
    wd(i, j, k-1, ivz) = wd(i, j, k-1, ivz) - sk(i, j, k-1, 1)*wwxd
    wd(i+1, j, k, ivy) = wd(i+1, j, k, ivy) + si(i, j, k, 3)*vvzd
    wd(i-1, j, k, ivy) = wd(i-1, j, k, ivy) - si(i-1, j, k, 3)*vvzd
    wd(i, j+1, k, ivy) = wd(i, j+1, k, ivy) + sj(i, j, k, 3)*vvzd
    wd(i, j, k+1, ivy) = wd(i, j, k+1, ivy) + sk(i, j, k, 3)*vvzd
    wd(i, j-1, k, ivy) = wd(i, j-1, k, ivy) - sj(i, j-1, k, 3)*vvzd
    wd(i, j, k-1, ivy) = wd(i, j, k-1, ivy) - sk(i, j, k-1, 3)*vvzd
    wd(i+1, j, k, ivy) = wd(i+1, j, k, ivy) + si(i, j, k, 1)*vvxd
    wd(i-1, j, k, ivy) = wd(i-1, j, k, ivy) - si(i-1, j, k, 1)*vvxd
    wd(i, j+1, k, ivy) = wd(i, j+1, k, ivy) + sj(i, j, k, 1)*vvxd
    wd(i, j, k+1, ivy) = wd(i, j, k+1, ivy) + sk(i, j, k, 1)*vvxd
    wd(i, j-1, k, ivy) = wd(i, j-1, k, ivy) - sj(i, j-1, k, 1)*vvxd
    wd(i, j, k-1, ivy) = wd(i, j, k-1, ivy) - sk(i, j, k-1, 1)*vvxd
    wd(i+1, j, k, ivx) = wd(i+1, j, k, ivx) + si(i, j, k, 3)*uuzd
    wd(i-1, j, k, ivx) = wd(i-1, j, k, ivx) - si(i-1, j, k, 3)*uuzd
    wd(i, j+1, k, ivx) = wd(i, j+1, k, ivx) + sj(i, j, k, 3)*uuzd
    wd(i, j, k+1, ivx) = wd(i, j, k+1, ivx) + sk(i, j, k, 3)*uuzd
    wd(i, j-1, k, ivx) = wd(i, j-1, k, ivx) - sj(i, j-1, k, 3)*uuzd
    wd(i, j, k-1, ivx) = wd(i, j, k-1, ivx) - sk(i, j, k-1, 3)*uuzd
    wd(i+1, j, k, ivx) = wd(i+1, j, k, ivx) + si(i, j, k, 2)*uuyd
    wd(i-1, j, k, ivx) = wd(i-1, j, k, ivx) - si(i-1, j, k, 2)*uuyd
    wd(i, j+1, k, ivx) = wd(i, j+1, k, ivx) + sj(i, j, k, 2)*uuyd
    wd(i, j, k+1, ivx) = wd(i, j, k+1, ivx) + sk(i, j, k, 2)*uuyd
    wd(i, j-1, k, ivx) = wd(i, j-1, k, ivx) - sj(i, j-1, k, 2)*uuyd
    wd(i, j, k-1, ivx) = wd(i, j, k-1, ivx) - sk(i, j, k-1, 2)*uuyd
  end do
end subroutine prodwmag2_fast_b
