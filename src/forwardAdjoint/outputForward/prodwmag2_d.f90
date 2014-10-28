   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.10 (r5363) -  9 Sep 2014 09:53
   !
   !  Differentiation of prodwmag2 in forward (tangent) mode (with options i4 dr8 r8):
   !   variations   of useful results: *dw
   !   with respect to varying inputs: *w *vol *si *sj *sk
   !   Plus diff mem management of: dw:in w:in vol:in si:in sj:in
   !                sk:in
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          prodWmag2.f90                                   *
   !      * Author:        Georgi Kalitzin, Edwin van der Weide            *
   !      * Starting date: 06-23-2003                                      *
   !      * Last modified: 06-12-2005                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   SUBROUTINE PRODWMAG2_D()
   !
   !      ******************************************************************
   !      *                                                                *
   !      * prodWmag2 computes the term:                                   *
   !      *    2*oij*oij  with oij=0.5*(duidxj - dujdxi).                  *
   !      * This is equal to the magnitude squared of the vorticity.       *
   !      * It is assumed that the pointer vort, stored in turbMod, is     *
   !      * already set to the correct entry.                              *
   !      *                                                                *
   !      ******************************************************************
   !
   USE BLOCKPOINTERS_D
   USE FLOWVARREFSTATE
   USE SECTION
   USE TURBMOD
   IMPLICIT NONE
   !
   !      Local variables.
   !
   INTEGER :: i, j, k
   REAL(kind=realtype) :: uy, uz, vx, vz, wx, wy
   REAL(kind=realtype) :: uyd, uzd, vxd, vzd, wxd, wyd
   REAL(kind=realtype) :: fact, vortx, vorty, vortz
   REAL(kind=realtype) :: factd, vortxd, vortyd, vortzd
   REAL(kind=realtype) :: omegax, omegay, omegaz
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Determine the non-dimensional wheel speed of this block.
   omegax = timeref*sections(sectionid)%rotrate(1)
   omegay = timeref*sections(sectionid)%rotrate(2)
   omegaz = timeref*sections(sectionid)%rotrate(3)
   dwd = 0.0_8
   ! Loop over the cell centers of the given block. It may be more
   ! efficient to loop over the faces and to scatter the gradient,
   ! but in that case the gradients for u, v and w must be stored.
   ! In the current approach no extra memory is needed.
   DO k=2,kl
   DO j=2,jl
   DO i=2,il
   ! Compute the necessary derivatives of u in the cell center.
   ! Use is made of the fact that the surrounding normals sum up
   ! to zero, such that the cell i,j,k does not give a
   ! contribution. The gradient is scaled by a factor 2*vol.
   uyd = wd(i+1, j, k, ivx)*si(i, j, k, 2) + w(i+1, j, k, ivx)*sid(&
   &         i, j, k, 2) - wd(i-1, j, k, ivx)*si(i-1, j, k, 2) - w(i-1, j, &
   &         k, ivx)*sid(i-1, j, k, 2) + wd(i, j+1, k, ivx)*sj(i, j, k, 2) &
   &         + w(i, j+1, k, ivx)*sjd(i, j, k, 2) - wd(i, j-1, k, ivx)*sj(i&
   &         , j-1, k, 2) - w(i, j-1, k, ivx)*sjd(i, j-1, k, 2) + wd(i, j, &
   &         k+1, ivx)*sk(i, j, k, 2) + w(i, j, k+1, ivx)*skd(i, j, k, 2) -&
   &         wd(i, j, k-1, ivx)*sk(i, j, k-1, 2) - w(i, j, k-1, ivx)*skd(i&
   &         , j, k-1, 2)
   uy = w(i+1, j, k, ivx)*si(i, j, k, 2) - w(i-1, j, k, ivx)*si(i-1&
   &         , j, k, 2) + w(i, j+1, k, ivx)*sj(i, j, k, 2) - w(i, j-1, k, &
   &         ivx)*sj(i, j-1, k, 2) + w(i, j, k+1, ivx)*sk(i, j, k, 2) - w(i&
   &         , j, k-1, ivx)*sk(i, j, k-1, 2)
   uzd = wd(i+1, j, k, ivx)*si(i, j, k, 3) + w(i+1, j, k, ivx)*sid(&
   &         i, j, k, 3) - wd(i-1, j, k, ivx)*si(i-1, j, k, 3) - w(i-1, j, &
   &         k, ivx)*sid(i-1, j, k, 3) + wd(i, j+1, k, ivx)*sj(i, j, k, 3) &
   &         + w(i, j+1, k, ivx)*sjd(i, j, k, 3) - wd(i, j-1, k, ivx)*sj(i&
   &         , j-1, k, 3) - w(i, j-1, k, ivx)*sjd(i, j-1, k, 3) + wd(i, j, &
   &         k+1, ivx)*sk(i, j, k, 3) + w(i, j, k+1, ivx)*skd(i, j, k, 3) -&
   &         wd(i, j, k-1, ivx)*sk(i, j, k-1, 3) - w(i, j, k-1, ivx)*skd(i&
   &         , j, k-1, 3)
   uz = w(i+1, j, k, ivx)*si(i, j, k, 3) - w(i-1, j, k, ivx)*si(i-1&
   &         , j, k, 3) + w(i, j+1, k, ivx)*sj(i, j, k, 3) - w(i, j-1, k, &
   &         ivx)*sj(i, j-1, k, 3) + w(i, j, k+1, ivx)*sk(i, j, k, 3) - w(i&
   &         , j, k-1, ivx)*sk(i, j, k-1, 3)
   ! Idem for the gradient of v.
   vxd = wd(i+1, j, k, ivy)*si(i, j, k, 1) + w(i+1, j, k, ivy)*sid(&
   &         i, j, k, 1) - wd(i-1, j, k, ivy)*si(i-1, j, k, 1) - w(i-1, j, &
   &         k, ivy)*sid(i-1, j, k, 1) + wd(i, j+1, k, ivy)*sj(i, j, k, 1) &
   &         + w(i, j+1, k, ivy)*sjd(i, j, k, 1) - wd(i, j-1, k, ivy)*sj(i&
   &         , j-1, k, 1) - w(i, j-1, k, ivy)*sjd(i, j-1, k, 1) + wd(i, j, &
   &         k+1, ivy)*sk(i, j, k, 1) + w(i, j, k+1, ivy)*skd(i, j, k, 1) -&
   &         wd(i, j, k-1, ivy)*sk(i, j, k-1, 1) - w(i, j, k-1, ivy)*skd(i&
   &         , j, k-1, 1)
   vx = w(i+1, j, k, ivy)*si(i, j, k, 1) - w(i-1, j, k, ivy)*si(i-1&
   &         , j, k, 1) + w(i, j+1, k, ivy)*sj(i, j, k, 1) - w(i, j-1, k, &
   &         ivy)*sj(i, j-1, k, 1) + w(i, j, k+1, ivy)*sk(i, j, k, 1) - w(i&
   &         , j, k-1, ivy)*sk(i, j, k-1, 1)
   vzd = wd(i+1, j, k, ivy)*si(i, j, k, 3) + w(i+1, j, k, ivy)*sid(&
   &         i, j, k, 3) - wd(i-1, j, k, ivy)*si(i-1, j, k, 3) - w(i-1, j, &
   &         k, ivy)*sid(i-1, j, k, 3) + wd(i, j+1, k, ivy)*sj(i, j, k, 3) &
   &         + w(i, j+1, k, ivy)*sjd(i, j, k, 3) - wd(i, j-1, k, ivy)*sj(i&
   &         , j-1, k, 3) - w(i, j-1, k, ivy)*sjd(i, j-1, k, 3) + wd(i, j, &
   &         k+1, ivy)*sk(i, j, k, 3) + w(i, j, k+1, ivy)*skd(i, j, k, 3) -&
   &         wd(i, j, k-1, ivy)*sk(i, j, k-1, 3) - w(i, j, k-1, ivy)*skd(i&
   &         , j, k-1, 3)
   vz = w(i+1, j, k, ivy)*si(i, j, k, 3) - w(i-1, j, k, ivy)*si(i-1&
   &         , j, k, 3) + w(i, j+1, k, ivy)*sj(i, j, k, 3) - w(i, j-1, k, &
   &         ivy)*sj(i, j-1, k, 3) + w(i, j, k+1, ivy)*sk(i, j, k, 3) - w(i&
   &         , j, k-1, ivy)*sk(i, j, k-1, 3)
   ! And for the gradient of w.
   wxd = wd(i+1, j, k, ivz)*si(i, j, k, 1) + w(i+1, j, k, ivz)*sid(&
   &         i, j, k, 1) - wd(i-1, j, k, ivz)*si(i-1, j, k, 1) - w(i-1, j, &
   &         k, ivz)*sid(i-1, j, k, 1) + wd(i, j+1, k, ivz)*sj(i, j, k, 1) &
   &         + w(i, j+1, k, ivz)*sjd(i, j, k, 1) - wd(i, j-1, k, ivz)*sj(i&
   &         , j-1, k, 1) - w(i, j-1, k, ivz)*sjd(i, j-1, k, 1) + wd(i, j, &
   &         k+1, ivz)*sk(i, j, k, 1) + w(i, j, k+1, ivz)*skd(i, j, k, 1) -&
   &         wd(i, j, k-1, ivz)*sk(i, j, k-1, 1) - w(i, j, k-1, ivz)*skd(i&
   &         , j, k-1, 1)
   wx = w(i+1, j, k, ivz)*si(i, j, k, 1) - w(i-1, j, k, ivz)*si(i-1&
   &         , j, k, 1) + w(i, j+1, k, ivz)*sj(i, j, k, 1) - w(i, j-1, k, &
   &         ivz)*sj(i, j-1, k, 1) + w(i, j, k+1, ivz)*sk(i, j, k, 1) - w(i&
   &         , j, k-1, ivz)*sk(i, j, k-1, 1)
   wyd = wd(i+1, j, k, ivz)*si(i, j, k, 2) + w(i+1, j, k, ivz)*sid(&
   &         i, j, k, 2) - wd(i-1, j, k, ivz)*si(i-1, j, k, 2) - w(i-1, j, &
   &         k, ivz)*sid(i-1, j, k, 2) + wd(i, j+1, k, ivz)*sj(i, j, k, 2) &
   &         + w(i, j+1, k, ivz)*sjd(i, j, k, 2) - wd(i, j-1, k, ivz)*sj(i&
   &         , j-1, k, 2) - w(i, j-1, k, ivz)*sjd(i, j-1, k, 2) + wd(i, j, &
   &         k+1, ivz)*sk(i, j, k, 2) + w(i, j, k+1, ivz)*skd(i, j, k, 2) -&
   &         wd(i, j, k-1, ivz)*sk(i, j, k-1, 2) - w(i, j, k-1, ivz)*skd(i&
   &         , j, k-1, 2)
   wy = w(i+1, j, k, ivz)*si(i, j, k, 2) - w(i-1, j, k, ivz)*si(i-1&
   &         , j, k, 2) + w(i, j+1, k, ivz)*sj(i, j, k, 2) - w(i, j-1, k, &
   &         ivz)*sj(i, j-1, k, 2) + w(i, j, k+1, ivz)*sk(i, j, k, 2) - w(i&
   &         , j, k-1, ivz)*sk(i, j, k-1, 2)
   ! Compute the three components of the vorticity vector.
   ! Substract the part coming from the rotating frame.
   factd = -(half*vold(i, j, k)/vol(i, j, k)**2)
   fact = half/vol(i, j, k)
   vortxd = factd*(wy-vz) + fact*(wyd-vzd)
   vortx = fact*(wy-vz) - two*omegax
   vortyd = factd*(uz-wx) + fact*(uzd-wxd)
   vorty = fact*(uz-wx) - two*omegay
   vortzd = factd*(vx-uy) + fact*(vxd-uyd)
   vortz = fact*(vx-uy) - two*omegaz
   ! Compute the magnitude squared of the vorticity.
   dwd(i, j, k, ivort) = 2*vortx*vortxd + 2*vorty*vortyd + 2*vortz*&
   &         vortzd
   dw(i, j, k, ivort) = vortx**2 + vorty**2 + vortz**2
   END DO
   END DO
   END DO
   END SUBROUTINE PRODWMAG2_D
