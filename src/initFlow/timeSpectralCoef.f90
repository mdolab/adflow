!
!      ******************************************************************
!      *                                                                *
!      * File:          timeSpectralCoef.f90                            *
!      * Author:        Edwin van der Weide, Arathi K. Gopinath         *
!      * Starting date: 07-29-2004                                      *
!      * Last modified: 06-28-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine timeSpectralCoef(coefSpectral, matrixCoefSpectral, &
                                   diagMatCoefSpectral)
!
!      ******************************************************************
!      *                                                                *
!      * timeSpectralCoef computes the time integration coefficients    *
!      * for the time spectral method. As it is possible that sections  *
!      * have different periodic times these coefficients are           *
!      * determined for all the sections. For vector quantities, such   *
!      * as momentum, these coefficients can also be different due to   *
!      * rotation and the fact that only a part of the wheel is         *
!      * simulated.                                                     *
!      *                                                                *
!      ******************************************************************
!
       use flowVarRefState
       use inputTimeSpectral
       use section
       implicit none
!
!      Subroutine arguments.
!
       real(kind=realType),                               &
           dimension(nSections,nTimeIntervalsSpectral-1), &
                                     intent(out) :: coefSpectral
       real(kind=realType),                                    &
            dimension(nSections,nTimeIntervalsSpectral-1,3,3), &
                                     intent(out) :: matrixCoefSpectral
       real(kind=realType), dimension(nSections,3,3), &
                                   intent(out) :: diagMatCoefSpectral
!
!      Local variables.
!
       integer(kind=intType) :: pp, nn, mm, ii, i, j, ntot
       real(kind=realType)   :: coef, dAngle, angle, fact, slicesFact

       real(kind=realType), dimension(3,3) :: rotMat, tmp
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of sections.

       sectionLoop: do mm=1,nSections

         ! Initialize dAngle (smallest angle in the cotangent function)
         ! and coef, which is the multiplication factor in front of the
         ! cotangent/cosecant function. Coef is a combination of the 1/2
         ! and the 2*pi/timePeriod

         dAngle = pi/real(nTimeIntervalsSpectral,realType)
         coef   = pi*timeRef/sections(mm)%timePeriod

         ! Computation of the scalar coefficients.

         scalarLoop: do nn=1,(nTimeIntervalsSpectral-1)

           angle = nn*dAngle

           ! The coefficient for an odd and even number of time
           ! instances are different; the former is 1/sin, the
           ! latter cos/sin or 1/tan.

           coefSpectral(mm,nn) = coef/sin(angle)

           if (mod(nTimeIntervalsSpectral,2_intType) == 0) & 
             coefSpectral(mm,nn) = coefSpectral(mm,nn)*cos(angle)

           ! Negate coef for the next spectral coefficient.

           coef = -coef

         enddo scalarLoop

         ! Initialize dAngle to the smallest angle in the cotangent
         ! or cosecant function. Now this angle is for the entire wheel, 
         ! i.e. the number of slices must be taken into account.

         ntot   = nTimeIntervalsSpectral*sections(mm)%nSlices
         dAngle = pi/real(ntot,realType)

         ! Initialize the rotation matrix to the unity matrix.

         rotMat(1,1) = one;  rotMat(1,2) = zero; rotMat(1,3) = zero
         rotMat(2,1) = zero; rotMat(2,2) = one;  rotMat(2,3) = zero
         rotMat(3,1) = zero; rotMat(3,2) = zero; rotMat(3,3) = one

         ! Loop over the number of spectral coefficient to initialize the
         ! matrix coefficients; this is basically pp == 0 in the loop
         ! over the number slices. Use is made of the fact that the
         ! rotation matrix is the identity for pp == 0.
         ! coef changes sign at every time instance
 
         slicesFact = one/real(sections(mm)%nSlices,realType)
         fact = one

         do nn=1,(nTimeIntervalsSpectral-1)

          ! Determine the scalar coefficient. This value depends now
          ! whether the total number of time instances in the wheel is
          ! odd or even.

           angle = nn*dAngle
           coef  = one/sin(angle)

           if (mod(ntot,2_intType) == 0) &
             coef = coef*cos(angle)

           coef = coef*fact*slicesFact

           ! The first part of matrixCoefSpectral is a diagonal matrix,
           ! because this indicates the contribution of the current
           ! slice to the time derivative.

           matrixCoefSpectral(mm,nn,1,1) = coef
           matrixCoefSpectral(mm,nn,1,2) = zero
           matrixCoefSpectral(mm,nn,1,3) = zero

           matrixCoefSpectral(mm,nn,2,1) = zero
           matrixCoefSpectral(mm,nn,2,2) = coef
           matrixCoefSpectral(mm,nn,2,3) = zero

           matrixCoefSpectral(mm,nn,3,1) = zero
           matrixCoefSpectral(mm,nn,3,2) = zero
           matrixCoefSpectral(mm,nn,3,3) = coef

           fact = -fact

         enddo

         ! Initialize diagMatCoefSpectral to zero, because the
         ! starting index in the loop over the number of slices -1 is
         ! 1, i.e. the slice where the actual computation takes places
         ! does not contribute to diagMatCoefSpectral.

         do j=1,3
           do i=1,3
             diagMatCoefSpectral(mm,i,j) = zero
           enddo
         enddo

         ! Loop over the additional slices which complete an entire
         ! revolution. To be able to compute the coefficients a bit
         ! easier the loop runs from 1 to nSlices-1 and not from
         ! 2 to nSlices.

         slicesLoop: do pp=1,(sections(mm)%nSlices-1)

           ! Compute the rotation matrix for this slice. This is the
           ! old one multiplied by the transformation matrix going from
           ! one slices to the next. Use tmp as temporary storage.

           do j=1,3
             do i=1,3
               tmp(i,j) = rotMatrixSpectral(mm,i,1)*rotMat(1,j) &
                        + rotMatrixSpectral(mm,i,2)*rotMat(2,j) &
                        + rotMatrixSpectral(mm,i,3)*rotMat(3,j)
             enddo
           enddo

           rotMat = tmp

           slicesFact = one/real(sections(mm)%nSlices,realType)

           ! Loop over the number of spectral coefficients and update
           ! matrixCoefSpectral. The multiplication with (-1)**nn 
           ! takes place here too.

           ! Multiply also by the term (-1)**(pN+1) 

           fact = one
           if (mod(pp*nTimeIntervalsSpectral,2_intType) /= 0) &
             fact = -one
           slicesFact = fact*slicesFact

           fact = one
           ii = pp*nTimeIntervalsSpectral
           do nn=1,(nTimeIntervalsSpectral-1)

             ! Compute the coefficient multiplying the rotation matrix.
             ! Again make a distinction between an odd and an even
             ! number of time instances for the entire wheel.

             angle = (nn+ii)*dAngle
             coef  = one/sin(angle)

             if (mod(ntot,2_intType) == 0) &
               coef = coef*cos(angle)

             coef = coef*fact*slicesFact

             ! Update matrixCoefSpectral.

             do j=1,3
               do i=1,3
                 matrixCoefSpectral(mm,nn,i,j) = &
                   matrixCoefSpectral(mm,nn,i,j) + coef*rotMat(i,j)
               enddo
             enddo

             fact = -fact

           enddo

           ! Update diagMatCoefSpectral. Also here the distinction
           ! between odd and even number of time instances.
           
           angle = ii*dAngle
           coef  = one/sin(angle)

           if (mod(ntot,2_intType) == 0) &
             coef = coef*cos(angle)

           coef = coef*slicesFact

           do j=1,3
             do i=1,3
               diagMatCoefSpectral(mm,i,j) = &
                 diagMatCoefSpectral(mm,i,j) - coef*rotMat(i,j)
             enddo
           enddo

         enddo slicesLoop

         ! The matrix coefficients must be multiplied by the leading
         ! coefficient, which depends on the actual periodic time.

         coef = pi*timeRef/sections(mm)%timePeriod     
 
         do j=1,3
           do i=1,3
             diagMatCoefSpectral(mm,i,j) = &
                         coef*diagMatCoefSpectral(mm,i,j) 
           enddo
         enddo

         do nn=1,(nTimeIntervalsSpectral-1)

           do j=1,3
             do i=1,3
               matrixCoefSpectral(mm,nn,i,j) = &
                         coef*matrixCoefSpectral(mm,nn,i,j) 
             enddo
           enddo

         end do

       enddo sectionLoop

       end subroutine timeSpectralCoef
