   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.10 (r5363) -  9 Sep 2014 09:53
   !
   !  Differentiation of resetxxssrhodd2wallbwd in reverse (adjoint) mode (with options i4 dr8 r8 noISIZE):
   !   gradient     of useful results: *w *x *si *sj *sk rho1 rho2
   !                xx ss
   !   with respect to varying inputs: *w *x *si *sj *sk rho1 rho2
   !                xx ss
   !   Plus diff mem management of: w:in x:in si:in sj:in sk:in
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          resetxxssrhodd2WallBwd.f90                      *
   !      * Author:        Peter Zhoujie Lyu                               *
   !      * Starting date: 10-28-2014                                      *
   !      * Last modified: 10-28-2014                                      *
   !      *                                                                *
   !      ******************************************************************
   SUBROUTINE RESETXXSSRHODD2WALLBWD_B(nn, xx, xxd, ss, ssd, rho1, rho1d, &
   & rho2, rho2d, dd2wall)
   USE BCTYPES
   USE BLOCKPOINTERS
   USE FLOWVARREFSTATE
   USE INPUTPHYSICS
   IMPLICIT NONE
   !
   !      Subroutine arguments.
   !
   INTEGER(kind=inttype), INTENT(IN) :: nn
   REAL(kind=realtype), DIMENSION(imaxdim, jmaxdim) :: rho2, rho1
   REAL(kind=realtype), DIMENSION(imaxdim, jmaxdim) :: rho2d, rho1d
   REAL(kind=realtype), DIMENSION(imaxdim - 2, jmaxdim - 2) :: dd2wall
   REAL(kind=realtype), DIMENSION(imaxdim, jmaxdim, 3) :: ss
   REAL(kind=realtype), DIMENSION(imaxdim, jmaxdim, 3) :: ssd
   REAL(kind=realtype), DIMENSION(imaxdim + 1, jmaxdim + 1, 3) :: xx
   REAL(kind=realtype), DIMENSION(imaxdim+1, jmaxdim+1, 3) :: xxd
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Determine the face id on which the subface is located and set
   ! the pointers accordinly.
   SELECT CASE  (bcfaceid(nn)) 
   CASE (imin) 
   xxd(1:je+1, 1:ke+1, :) = xxd(1:je+1, 1:ke+1, :) + xd(1, 0:je, 0:ke, &
   &     :)
   xd(1, 0:je, 0:ke, :) = 0.0_8
   ssd(1:je, 1:ke, :) = ssd(1:je, 1:ke, :) + sid(1, 1:je, 1:ke, :)
   sid(1, 1:je, 1:ke, :) = 0.0_8
   rho1d(1:je, 1:ke) = rho1d(1:je, 1:ke) + wd(1, 1:je, 1:ke, irho)
   wd(1, 1:je, 1:ke, irho) = 0.0_8
   rho2d(1:je, 1:ke) = rho2d(1:je, 1:ke) + wd(2, 1:je, 1:ke, irho)
   wd(2, 1:je, 1:ke, irho) = 0.0_8
   CASE (imax) 
   xxd(1:je+1, 1:ke+1, :) = xxd(1:je+1, 1:ke+1, :) + xd(il, 0:je, 0:ke&
   &     , :)
   xd(il, 0:je, 0:ke, :) = 0.0_8
   ssd(1:je, 1:ke, :) = ssd(1:je, 1:ke, :) + sid(il, 1:je, 1:ke, :)
   sid(il, 1:je, 1:ke, :) = 0.0_8
   rho1d(1:je, 1:ke) = rho1d(1:je, 1:ke) + wd(ie, 1:je, 1:ke, irho)
   wd(ie, 1:je, 1:ke, irho) = 0.0_8
   rho2d(1:je, 1:ke) = rho2d(1:je, 1:ke) + wd(il, 1:je, 1:ke, irho)
   wd(il, 1:je, 1:ke, irho) = 0.0_8
   CASE (jmin) 
   xxd(1:ie+1, 1:ke+1, :) = xxd(1:ie+1, 1:ke+1, :) + xd(0:ie, 1, 0:ke, &
   &     :)
   xd(0:ie, 1, 0:ke, :) = 0.0_8
   ssd(1:ie, 1:ke, :) = ssd(1:ie, 1:ke, :) + sjd(1:ie, 1, 1:ke, :)
   sjd(1:ie, 1, 1:ke, :) = 0.0_8
   rho1d(1:ie, 1:ke) = rho1d(1:ie, 1:ke) + wd(1:ie, 1, 1:ke, irho)
   wd(1:ie, 1, 1:ke, irho) = 0.0_8
   rho2d(1:ie, 1:ke) = rho2d(1:ie, 1:ke) + wd(1:ie, 2, 1:ke, irho)
   wd(1:ie, 2, 1:ke, irho) = 0.0_8
   CASE (jmax) 
   xxd(1:ie+1, 1:ke+1, :) = xxd(1:ie+1, 1:ke+1, :) + xd(0:ie, jl, 0:ke&
   &     , :)
   xd(0:ie, jl, 0:ke, :) = 0.0_8
   ssd(1:ie, 1:ke, :) = ssd(1:ie, 1:ke, :) + sjd(1:ie, jl, 1:ke, :)
   sjd(1:ie, jl, 1:ke, :) = 0.0_8
   rho1d(1:ie, 1:ke) = rho1d(1:ie, 1:ke) + wd(1:ie, je, 1:ke, irho)
   wd(1:ie, je, 1:ke, irho) = 0.0_8
   rho2d(1:ie, 1:ke) = rho2d(1:ie, 1:ke) + wd(1:ie, jl, 1:ke, irho)
   wd(1:ie, jl, 1:ke, irho) = 0.0_8
   CASE (kmin) 
   xxd(1:ie+1, 1:je+1, :) = xxd(1:ie+1, 1:je+1, :) + xd(0:ie, 0:je, 1, &
   &     :)
   xd(0:ie, 0:je, 1, :) = 0.0_8
   ssd(1:ie, 1:je, :) = ssd(1:ie, 1:je, :) + skd(1:ie, 1:je, 1, :)
   skd(1:ie, 1:je, 1, :) = 0.0_8
   rho1d(1:ie, 1:je) = rho1d(1:ie, 1:je) + wd(1:ie, 1:je, 1, irho)
   wd(1:ie, 1:je, 1, irho) = 0.0_8
   rho2d(1:ie, 1:je) = rho2d(1:ie, 1:je) + wd(1:ie, 1:je, 2, irho)
   wd(1:ie, 1:je, 2, irho) = 0.0_8
   CASE (kmax) 
   xxd(1:ie+1, 1:je+1, :) = xxd(1:ie+1, 1:je+1, :) + xd(0:ie, 0:je, kl&
   &     , :)
   xd(0:ie, 0:je, kl, :) = 0.0_8
   ssd(1:ie, 1:je, :) = ssd(1:ie, 1:je, :) + skd(1:ie, 1:je, kl, :)
   skd(1:ie, 1:je, kl, :) = 0.0_8
   rho1d(1:ie, 1:je) = rho1d(1:ie, 1:je) + wd(1:ie, 1:je, ke, irho)
   wd(1:ie, 1:je, ke, irho) = 0.0_8
   rho2d(1:ie, 1:je) = rho2d(1:ie, 1:je) + wd(1:ie, 1:je, kl, irho)
   wd(1:ie, 1:je, kl, irho) = 0.0_8
   END SELECT
   END SUBROUTINE RESETXXSSRHODD2WALLBWD_B
