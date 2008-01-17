!
!      ******************************************************************
!      *                                                                *
!      * File:          spectralInterpolCoef.f90                        *
!      * Author:        Edwin van der Weide, Arathi K. Gopinath.        *
!      * Starting date: 08-15-2004                                      *
!      * Last modified: 04-29-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine spectralInterpolCoef(nsps, t, alpScal, alpMat)
!
!      ******************************************************************
!      *                                                                *
!      * spectralInterpolCoef determines the scalar and matrix          *
!      * spectral interpolation coefficients for the given number of    *
!      * spectral solutions for the given t, where t is the ratio of    *
!      * the time and the periodic interval time. Note that the index   *
!      * of the spectral solutions of both alpScal and alpMat start     *
!      * at 0. In this way these coefficients are easier to determine.  *
!      *                                                                *
!      ******************************************************************
!
       use constants
       use inputTimeSpectral
       use section
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nsps
       real(kind=realType),   intent(in) :: t

       real(kind=realType), dimension(0:nsps-1), intent(out) :: alpScal
       real(kind=realType), dimension(nSections,0:nsps-1,3,3), &
                                                  intent(out) :: alpMat
!
!      Local variables.
!
       integer(kind=intType) :: jj, nn, j, p, r, nhalfM1, m, mhalfM1

       real(kind=realType) :: nspsInv, mInv, tm, alp

       real(kind=realType), dimension(3,3) :: rp, tmp
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
!      ******************************************************************
!      *                                                                *
!      * Scalar coefficients.                                           *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of spectral solutions to compute the
       ! coefficients. Note that the loop starts at 0.

       if (mod(nsps,2).eq.0) then
         nhalfM1 = nsps/2 - 1
       else
         nhalfM1 = (nsps-1)/2
       endif

       nspsInv = one/real(nsps,realType)

       do j=0,(nsps-1)
         if (mod(nsps,2).eq.0) then
           alpScal(j) = one + cos(j*pi)*cos(nsps*pi*t)
         else 
           alpScal(j) = one + cos(j*pi*(nsps+1)/nsps)*cos((nsps+1)*pi*t)
         endif

         do r=1,nhalfM1
           alpScal(j) = alpScal(j)                                  &
                      + two*cos(r*j*two*pi*nspsInv)*cos(r*two*pi*t) &
                      + two*sin(r*j*two*pi*nspsInv)*sin(r*two*pi*t)
         enddo

         alpScal(j) = alpScal(j)*nspsInv

       enddo
!
!      ******************************************************************
!      *                                                                *
!      * Matrix coefficients. These are (can be) different for every    *
!      * section and they must therefore be determined for every        *
!      * section.                                                       *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of sections in the grid.

       sectionLoop: do nn=1,nSections

         ! Compute the numbers for the entire wheel for this section.
         ! Note that also t must be adapted, because t is a ratio between
         ! the actual time and the periodic time.

         m        = nsps*sections(nn)%nSlices
         if (mod(m,2).eq.0) then
           mhalfM1 = m/2 - 1
         else 
           mhalfM1 = (m-1)/2 
         endif
         mInv    = one/real(m,realType)
         tm       = t/real(sections(nn)%nSlices,realType)

         ! Loop over the number of spectral solutions.

         spectralLoop: do jj=0,(nsps-1)

           ! Initialize the matrix coefficients to zero and the matrix
           ! rp to the identity matrix. Rp is the rotation matrix of this
           ! section to the power p, which starts at 0, i.e. rp = i.

           alpMat(nn,jj,1,1) = zero
           alpMat(nn,jj,1,2) = zero
           alpMat(nn,jj,1,3) = zero

           alpMat(nn,jj,2,1) = zero
           alpMat(nn,jj,2,2) = zero
           alpMat(nn,jj,2,3) = zero

           alpMat(nn,jj,3,1) = zero
           alpMat(nn,jj,3,2) = zero
           alpMat(nn,jj,3,3) = zero

           rp(1,1) = one
           rp(1,2) = zero
           rp(1,3) = zero

           rp(2,1) = zero
           rp(2,2) = one
           rp(2,3) = zero

           rp(3,1) = zero
           rp(3,2) = zero
           rp(3,3) = one

           ! Loop over the number of slices of this section. Note that
           ! this loop starts at zero, which simplifies the formulas.

           slicesLoop: do p=0,(sections(nn)%nSlices-1)

             ! Determine the index j, the index of alp in the entire
             ! wheel.

             j = jj + p*nsps

             ! Compute the scalar coefficient alp of the index j in
             ! the entire wheel.

             if (mod(m,2).eq.0) then
               alp = one + cos(j*pi)*cos(m*pi*tm)
             else
               alp = one + cos(j*pi*(m+1)/m)*cos((m+1)*pi*tm)
             endif
             do r=1,mhalfM1
               alp = alp + two*cos(r*j*two*pi*mInv)*cos(r*two*pi*tm) &
                   +       two*sin(r*j*two*pi*mInv)*sin(r*two*pi*tm)
             enddo

             alp = alp*mInv

             ! Update the matrix coefficient.

             do r=1,3
               do j=1,3
                 alpMat(nn,jj,r,j) = alpMat(nn,jj,r,j) + alp*rp(r,j)
               enddo
             enddo

             ! Multiply rp by the rotation matrix to obtain the correct
             ! matrix for the next slice. Use tmp as temporary storage.

             do r=1,3
               do j=1,3
                 tmp(r,j) = rp(r,1)*rotMatrixSpectral(nn,1,j) &
                          + rp(r,2)*rotMatrixSpectral(nn,2,j) &
                          + rp(r,3)*rotMatrixSpectral(nn,3,j)
               enddo
             enddo

             rp = tmp

           enddo slicesLoop
         enddo spectralLoop
       enddo sectionLoop

       end subroutine spectralInterpolCoef
