!
!      ******************************************************************
!      *                                                                *
!      * File:          vfScale.f90                                     *
!      * Author:        Georgi Kalitzin                                 *
!      * Starting date: 04-19-2004                                      *
!      * Last modified: 04-11-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine vfScale
!
!      ******************************************************************
!      *                                                                *
!      * time and length scale definition for v2f turbulence model. The *
!      * upper bound can be switched on by setting rvfB to .true.       *
!      * The strain squared is defined as: strain2 = 2 sij sij          *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use inputPhysics
       use paramTurb
       use turbMod
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k

       real(kind=realType) :: sqrt3
       real(kind=realType) :: tkea, tepa, tv2a, supi, rn2
       real(kind=realType) :: rsct, rscl2, rnu, rstrain
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Some constants in the model.

       sqrt3 = sqrt(three)

       ! Set the pointer for dvt in dw, such that the code is more
       ! readable. Also set the pointers for the production term,
       ! vorticity, strain and the time and lenght scale of v2f.

       dvt     => dw(1:,1:,1:,idvt:)
       prod    => dw(1:,1:,1:,iprod)
       sct     => dw(1:,1:,1:,isct)
       scl2    => dw(1:,1:,1:,iscl2)
       vort    => prod
       strain2 => prod
!
!      ******************************************************************
!      *                                                                *
!      * Production term.                                               *
!      *                                                                *
!      ******************************************************************
!
       select case (turbProd)
         case (strain)
           call prodSmag2

         case (vorticity)
           call prodWmag2

         case (katoLaunder)
           call prodKatoLaunder

       end select
!
!      ******************************************************************
!      *                                                                *
!      * Compute the length and time scale for all internal cells.      *
!      *                                                                *
!      ******************************************************************
!
       if( rvfB ) then

         do k=2,kl
           do j=2,jl
             do i=2,il

               ! Compute the time and length scale with upper bound

               rstrain     = sqrt(strain2(i,j,k))
               rnu         = rlv(i,j,k)/w(i,j,k,irho)
               tkea        = abs(w(i,j,k,itu1))
               tepa        = abs(w(i,j,k,itu2))
               tv2a        = abs(w(i,j,k,itu3))
               supi        = tepa*tkea/max(sqrt3*tv2a*rvfCmu*rstrain,eps)
               rn2         = rvfCn**2*(rnu*tepa)**1.5_realType

               rsct        = max(tkea,six*sqrt(rnu*tepa))
               sct(i,j,k)  = min(rsct,0.6_realType*supi)
               rscl2       = tkea*min(tkea**2,supi**2)
               scl2(i,j,k) = rvfCl**2*max(rscl2,rn2)

             enddo
           enddo
         enddo

       else

         do k=2,kl
           do j=2,jl
             do i=2,il

               ! Compute the time and length scale without upper bound
 
               rnu         = rlv(i,j,k)/w(i,j,k,irho)
               tkea        = abs(w(i,j,k,itu1))
               tepa        = abs(w(i,j,k,itu2))
               rn2         = rvfCn**2*(rnu*tepa)**1.5_realType

               rsct        = max(tkea,six*sqrt(rnu*tepa))
               sct(i,j,k)  = rsct
               rscl2       = tkea**3
               scl2(i,j,k) = rvfCl**2*max(rscl2,rn2)
             enddo
           enddo
         enddo
       endif

       end subroutine vfScale
