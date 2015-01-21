   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.10 (r5363) -  9 Sep 2014 09:53
   !
   !  Differentiation of resetssbwd in reverse (adjoint) mode (with options i4 dr8 r8 noISIZE):
   !   gradient     of useful results: *si *sj *sk
   !   with respect to varying inputs: *si *sj *sk ssi ssj ssk
   !   Plus diff mem management of: si:in sj:in sk:in
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          resetSSBwd.f90                                  *
   !      * Author:        Peter Zhoujie Lyu                               *
   !      * Starting date: 10-21-2014                                      *
   !      * Last modified: 10-21-2014                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   SUBROUTINE RESETSSBWD_B(nn, ssi, ssid, ssj, ssjd, ssk, sskd, ss)
   USE BCTYPES
   USE BLOCKPOINTERS
   USE FLOWVARREFSTATE
   IMPLICIT NONE
   !
   !      Subroutine arguments.
   !
   INTEGER(kind=inttype), INTENT(IN) :: nn
   REAL(kind=realtype), DIMENSION(imaxdim, jmaxdim, 3) :: ssi, ssj, ssk
   REAL(kind=realtype), DIMENSION(imaxdim, jmaxdim, 3) :: ssid, ssjd, &
   & sskd
   REAL(kind=realtype), DIMENSION(imaxdim, jmaxdim, 3) :: ss
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
   sskd = 0.0_8
   sskd(1:je, 0:ke, :) = sskd(1:je, 0:ke, :) + skd(2, 1:je, 0:ke, :)
   skd(2, 1:je, 0:ke, :) = 0.0_8
   ssjd = 0.0_8
   ssjd(0:je, 1:ke, :) = ssjd(0:je, 1:ke, :) + sjd(2, 0:je, 1:ke, :)
   sjd(2, 0:je, 1:ke, :) = 0.0_8
   ssid = 0.0_8
   ssid(1:je, 1:ke, :) = ssid(1:je, 1:ke, :) + sid(1, 1:je, 1:ke, :)
   sid(1, 1:je, 1:ke, :) = 0.0_8
   CASE (imax) 
   sskd = 0.0_8
   sskd(1:je, 0:ke, :) = sskd(1:je, 0:ke, :) + skd(il, 1:je, 0:ke, :)
   skd(il, 1:je, 0:ke, :) = 0.0_8
   ssjd = 0.0_8
   ssjd(0:je, 1:ke, :) = ssjd(0:je, 1:ke, :) + sjd(il, 0:je, 1:ke, :)
   sjd(il, 0:je, 1:ke, :) = 0.0_8
   ssid = 0.0_8
   ssid(1:je, 1:ke, :) = ssid(1:je, 1:ke, :) + sid(il, 1:je, 1:ke, :)
   sid(il, 1:je, 1:ke, :) = 0.0_8
   CASE (jmin) 
   sskd = 0.0_8
   sskd(1:ie, 0:ke, :) = sskd(1:ie, 0:ke, :) + skd(1:ie, 2, 0:ke, :)
   skd(1:ie, 2, 0:ke, :) = 0.0_8
   ssjd = 0.0_8
   ssjd(1:ie, 1:ke, :) = ssjd(1:ie, 1:ke, :) + sid(1:ie, 2, 1:ke, :)
   sid(1:ie, 2, 1:ke, :) = 0.0_8
   ssid = 0.0_8
   ssid(0:ie, 1:ke, :) = ssid(0:ie, 1:ke, :) + sjd(0:ie, 1, 1:ke, :)
   sjd(0:ie, 1, 1:ke, :) = 0.0_8
   CASE (jmax) 
   sskd = 0.0_8
   sskd(1:ie, 0:ke, :) = sskd(1:ie, 0:ke, :) + skd(1:ie, jl, 0:ke, :)
   skd(1:ie, jl, 0:ke, :) = 0.0_8
   ssjd = 0.0_8
   ssjd(1:ie, 1:ke, :) = ssjd(1:ie, 1:ke, :) + sid(1:ie, jl, 1:ke, :)
   sid(1:ie, jl, 1:ke, :) = 0.0_8
   ssid = 0.0_8
   ssid(0:ie, 1:ke, :) = ssid(0:ie, 1:ke, :) + sjd(0:ie, jl, 1:ke, :)
   sjd(0:ie, jl, 1:ke, :) = 0.0_8
   CASE (kmin) 
   sskd = 0.0_8
   sskd(1:ie, 1:je, :) = sskd(1:ie, 1:je, :) + sjd(1:ie, 1:je, 2, :)
   sjd(1:ie, 1:je, 2, :) = 0.0_8
   ssjd = 0.0_8
   ssjd(1:ie, 0:je, :) = ssjd(1:ie, 0:je, :) + sid(1:ie, 0:je, 2, :)
   sid(1:ie, 0:je, 2, :) = 0.0_8
   ssid = 0.0_8
   ssid(0:ie, 1:je, :) = ssid(0:ie, 1:je, :) + skd(0:ie, 1:je, 1, :)
   skd(0:ie, 1:je, 1, :) = 0.0_8
   CASE (kmax) 
   sskd = 0.0_8
   sskd(1:ie, 1:je, :) = sskd(1:ie, 1:je, :) + sjd(1:ie, 1:je, kl, :)
   sjd(1:ie, 1:je, kl, :) = 0.0_8
   ssjd = 0.0_8
   ssjd(1:ie, 0:je, :) = ssjd(1:ie, 0:je, :) + sid(1:ie, 0:je, kl, :)
   sid(1:ie, 0:je, kl, :) = 0.0_8
   ssid = 0.0_8
   ssid(0:ie, 1:je, :) = ssid(0:ie, 1:je, :) + skd(0:ie, 1:je, kl, :)
   skd(0:ie, 1:je, kl, :) = 0.0_8
   CASE DEFAULT
   ssid = 0.0_8
   ssjd = 0.0_8
   sskd = 0.0_8
   END SELECT
   END SUBROUTINE RESETSSBWD_B
