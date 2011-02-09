!
!      ******************************************************************
!      *                                                                *
!      * File:          writeLoglaw.f90                                 *
!      * Author:        Georgi Kalitzin                                 *
!      * Starting date: 08-24-2004                                      *
!      * Last modified: 04-12-2004                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine writeLoglaw
!
!      ******************************************************************
!      *                                                                *
!      * writeLoglaw writes a profile of the velocity and turbulence    *
!      * variables in plus units for a flat plate with i in x, j in y   *
!      * and k in z direction at a hardcoded location in x direction    *
!      * and j=2, k=2. It is assumed that the pointers in               *
!      * blockPointers already point to the correct block.              *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use constants
       use inputPhysics
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: ii, i, j, k, nn

       real(kind=realType) :: utauF, rdF, ruF
       real(kind=realType) :: rtu1F, rtu2F, rtu3F, rtu4F, rrevF
       real(kind=realType) :: rarglog, rfunlog, rarglin, rfunlin
       real(kind=realType) :: t1, t2

       ! Use only for flat plate

       return

       k=2
       i=0
       j=2
       nn=1
       do ii=2,il
         if( half*(x(ii  ,j,k,1)+x(ii-1,j,k,1)) > 0.9_realType .and. &
             half*(x(ii-1,j,k,1)+x(ii-2,j,k,1)) < 0.9_realType ) i=ii
       enddo

!      I=137
!      i=98
!      i=110
!      i=118
!      if(x(i,2,2,1).Ge.1.) return

       if(i .ne. 0 .and. jl.gt.2) then
          write(8,3001)
 3001     format('#variables="y+","u+","ylog+","ulog+","ylin+","ulin+",&
                                   &"rev+","tu1+","tu2+","tu3+","tu4+"')
          if(wallFunctions) then
            utauF = viscSubface(nn)%utau(i,k)
          else
            utauF = sqrt( rlv(i,2,k)/w(i,2,k,irho)*abs(w(i,2,k,ivx)) &
                   /       d2Wall(i,2,k) )
          endif
          write(*,*) '# x(i,2,2,1),utauF=',x(i,2,2,1),utauF
          do j=1,jl
             rdF    = w(i,j,k,irho)*d2Wall(i,j,k)*utauF/rlv(i,j,k)
             if(j.eq.1) rdF = - w(i,2,k,irho)*d2Wall(i,2,k)*utauF/rlv(i,2,k)
             ruF    = w(i,j,k,ivx)/utauF
             rrevF  = rev(i,j,k)/rlv(i,j,k)
             select case (turbModel)

               case (spalartAllmaras)
                 rtu1F  = w(i,j,k,itu1)/rlv(i,j,k)*w(i,j,k,irho)
                 rtu2F  = zero
                 rtu3F  = zero
                 rtu4F  = zero

               case (komegaWilcox, komegaModified)
                 rtu1F  = w(i,j,k,itu1)/utauF**2
                 rtu2F  = w(i,j,k,itu2)*rlv(i,j,k)/w(i,j,k,irho)/utauF**2
                 rtu3F  = zero
                 rtu4F  = zero

               case (menterSST)
                 rtu1F  = w(i,j,k,itu1)/utauF**2
                 rtu2F  = w(i,j,k,itu2)*rlv(i,j,k)/w(i,j,k,irho)/utauF**2
                 rtu3F  = zero
                 rtu4F  = zero

               case (ktau)
                 rtu1F  = w(i,j,k,itu1)/utauF**2
                 rtu2F  = w(i,j,k,itu2)/rlv(i,j,k)*w(i,j,k,irho)*utauF**2
                 rtu3F  = zero
                 rtu4F  = zero

               case (v2f)
                 rtu1F  = w(i,j,k,itu1)/utauF**2
                 rtu2F  = w(i,j,k,itu2)*rlv(i,j,k)/w(i,j,k,irho)/utauF**4
                 rtu3F  = w(i,j,k,itu3)/utauF**2
                 rtu4F  = w(i,j,k,itu4)*rlv(i,j,k)/w(i,j,k,irho)/utauF**2
 
               case default
                 call terminate("writeLoglaw", &
                                "Turbulence model not implemented yet")

             end select

             t1 = exp((10.0_realType - 5.5_realType)*0.418_realType)
             t2 = exp((35.0_realType - 5.5_realType)*0.418_realType)

             rarglog = max(t1,min(t2,rdF))
             rfunlog = one/0.418_realType*log(rarglog) + 5.5_realType
             rarglin = min(rdF,15.0_realType)
             rfunlin = rarglin
             write(8,'(14e14.6)') rdF,ruF,rarglog,rfunlog,rarglin,rfunlin,&
                                  rrevF,rtu1F,rtu2F,rtu3F,rtu4F
          enddo
          close(8)
       endif

       end subroutine writeLoglaw
