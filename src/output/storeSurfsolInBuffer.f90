!
!      ******************************************************************
!      *                                                                *
!      * File:          storeSurfsolInBuffer.f90                        *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 05-19-2003                                      *
!      * Last modified: 07-14-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine storeSurfsolInBuffer(sps, buffer, nn, blockID,   &
                                       faceID, cellRange, solName, &
                                       viscousSubface)
!
!      ******************************************************************
!      *                                                                *
!      * storeSurfsolInBuffer stores the variable indicated by          *
!      * solName of the given block ID in the buffer. As the solution   *
!      * must be stored in the center of the boundary face the average  *
!      * value of the first internal cell and its corresponding halo is *
!      * computed. The counter nn is updated in this routine. However   *
!      * it is not initialized, because multiple contributions may be   *
!      * stored in buffer.                                              *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use BCTypes
       use cgnsNames
       use constants
       use flowVarRefState
       use inputPhysics
       use inputIO
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in)    :: sps, blockID, faceID
       integer(kind=intType), intent(inout) :: nn
       integer(kind=intType), dimension(3,2), intent(in) :: cellRange
       real(kind=realType), dimension(*), intent(out) :: buffer
       character(len=*), intent(in) :: solName
       logical, intent(in) :: viscousSubface
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k
       integer(kind=intType) :: ii, jj, mm, iiMax, jjMax, offVis

       integer(kind=intType), dimension(2,2) :: rangeFace
       integer(kind=intType), dimension(3,2) :: rangeCell

       integer(kind=intType), dimension(:,:), pointer :: viscPointer
       integer(kind=intType), dimension(:,:), pointer :: iblank2

       real(kind=realType) :: fact, gm1, ptotInf, ptot, psurf, rsurf
       real(kind=realType) :: usurf, vsurf, wsurf, m2surf, musurf
       real(kind=realType) :: fx, fy, fz, fn, a2Tot, a2, qw
       real(kind=realType) :: tauxx, tauyy, tauzz
       real(kind=realType) :: tauxy, tauxz, tauyz

       real(kind=realType), dimension(3) :: norm

       real(kind=realType), dimension(:,:,:), pointer :: ww1, ww2
       real(kind=realType), dimension(:,:),   pointer :: pp1, pp2
       real(kind=realType), dimension(:,:),   pointer :: gamma1, gamma2
       real(kind=realType), dimension(:,:),   pointer :: rlv1, rlv2
       real(kind=realType), dimension(:,:),   pointer :: dd2Wall
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set the pointers to this block.

       call setPointers(blockID, 1_intType, sps)

       ! Set the offset for the viscous data, such that the range is
       ! limited to the actual physical face. Viscous data, like skin
       ! friction, need gradIent information, which is not available
       ! in the halo's.

       offVis = 0
       if( storeRindLayer ) offVis = 1

       ! CellRange contains the range of the current block in the
       ! original cgns block. Substract the offset and store the local
       ! range in rangeCell.

       rangeCell(1,1) = cellRange(1,1) - iBegor + 1
       rangeCell(1,2) = cellRange(1,2) - iBegor + 1

       rangeCell(2,1) = cellRange(2,1) - jBegor + 1
       rangeCell(2,2) = cellRange(2,2) - jBegor + 1

       rangeCell(3,1) = cellRange(3,1) - kBegor + 1
       rangeCell(3,2) = cellRange(3,2) - kBegor + 1
!
!      ******************************************************************
!      *                                                                *
!      *          Viscous variables for a non-viscous wall.             *
!      *          Simply set the variables to zero and return.          *
!      *                                                                *
!      ******************************************************************
!
       if(.not. viscousSubface) then

         select case (solName)

           case (cgnsSkinFmag, cgnsStanton, cgnsYplus, &
                 cgnsSkinFx, cgnsSkinFy, cgnsSkinFz)

             ! Update the counter and set this entry of buffer to 0.

             do k=rangeCell(3,1),rangeCell(3,2)
               do j=rangeCell(2,1),rangeCell(2,2)
                 do i=rangeCell(1,1),rangeCell(1,2)
                   nn = nn + 1
                   buffer(nn) = zero
                 enddo
               enddo
             enddo

             ! Work has been done for this variable. So return.

             return

         end select
       endif
!
!      ******************************************************************
!      *                                                                *
!      * Determine the face on which the subface is located and set     *
!      * a couple of variables accordingly. In this way a generic       *
!      * treatment is possible and there is no need to repeat the code  *
!      * for each of the six block faces.                               *
!      * Note that for dd2Wall a slightly different notation must be    *
!      * used. Reason is that d2Wall starts at index 2, rather than 0.  *
!      *                                                                *
!      ******************************************************************
!
       select case (faceID)

         case (iMin)
           rangeFace(1,1:2) = rangeCell(2,1:2)
           rangeFace(2,1:2) = rangeCell(3,1:2)
           iiMax = jl; jjMax = kl

           ww1    => w(1,1:,1:,:);   ww2    => w(2,1:,1:,:)
           pp1    => p(1,1:,1:);     pp2    => p(2,1:,1:)
           gamma1 => gamma(1,1:,1:); gamma2 => gamma(2,1:,1:)

           iblank2 => iblank(2,1:,1:)
           viscPointer => viscIminPointer

           if( viscous ) then
             rlv1 => rlv(1,1:,1:); rlv2 => rlv(2,1:,1:)
           endif

           if(equations == RANSEquations) dd2Wall => d2Wall(2,:,:)

         !===============================================================

         case (iMax)
           rangeFace(1,1:2) = rangeCell(2,1:2)
           rangeFace(2,1:2) = rangeCell(3,1:2)
           iiMax = jl; jjMax = kl

           ww1    => w(ie,1:,1:,:);   ww2    => w(il,1:,1:,:)
           pp1    => p(ie,1:,1:);     pp2    => p(il,1:,1:)
           gamma1 => gamma(ie,1:,1:); gamma2 => gamma(il,1:,1:)

           iblank2 => iblank(il,1:,1:)
           viscPointer => viscImaxPointer

           if( viscous ) then
             rlv1 => rlv(ie,1:,1:); rlv2 => rlv(il,1:,1:)
           endif

           if(equations == RANSEquations) dd2Wall => d2Wall(il,:,:)

         !===============================================================

         case (jMin)
           rangeFace(1,1:2) = rangeCell(1,1:2)
           rangeFace(2,1:2) = rangeCell(3,1:2)
           iiMax = il; jjMax = kl

           ww1    => w(1:,1,1:,:);   ww2    => w(1:,2,1:,:)
           pp1    => p(1:,1,1:);     pp2    => p(1:,2,1:)
           gamma1 => gamma(1:,1,1:); gamma2 => gamma(1:,2,1:)

           iblank2 => iblank(1:,2,1:)
           viscPointer => viscJminPointer

           if( viscous ) then
             rlv1 => rlv(1:,1,1:); rlv2 => rlv(1:,2,1:)
           endif

           if(equations == RANSEquations) dd2Wall => d2Wall(:,2,:)

         !===============================================================

         case (jMax)
           rangeFace(1,1:2) = rangeCell(1,1:2)
           rangeFace(2,1:2) = rangeCell(3,1:2)
           iiMax = il; jjMax = kl

           ww1    => w(1:,je,1:,:);   ww2    => w(1:,jl,1:,:)
           pp1    => p(1:,je,1:);     pp2    => p(1:,jl,1:)
           gamma1 => gamma(1:,je,1:); gamma2 => gamma(1:,jl,1:)

           iblank2 => iblank(1:,jl,1:)
           viscPointer => viscJmaxPointer

           if( viscous ) then
             rlv1 => rlv(1:,je,1:); rlv2 => rlv(1:,jl,1:)
           endif

           if(equations == RANSEquations) dd2Wall => d2Wall(:,jl,:)

         !===============================================================

         case (kMin)
           rangeFace(1,1:2) = rangeCell(1,1:2)
           rangeFace(2,1:2) = rangeCell(2,1:2)
           iiMax = il; jjMax = jl

           ww1    => w(1:,1:,1,:);   ww2    => w(1:,1:,2,:)
           pp1    => p(1:,1:,1);     pp2    => p(1:,1:,2)
           gamma1 => gamma(1:,1:,1); gamma2 => gamma(1:,1:,2)

           iblank2 => iblank(1:,1:,2)
           viscPointer => viscKminPointer

           if( viscous ) then
             rlv1 => rlv(1:,1:,1); rlv2 => rlv(1:,1:,2)
           endif

           if(equations == RANSEquations) dd2Wall => d2Wall(:,:,2)

         !===============================================================

         case (kMax)
           rangeFace(1,1:2) = rangeCell(1,1:2)
           rangeFace(2,1:2) = rangeCell(2,1:2)
           iiMax = il; jjMax = jl

           ww1    => w(1:,1:,ke,:);   ww2    => w(1:,1:,kl,:)
           pp1    => p(1:,1:,ke);     pp2    => p(1:,1:,kl)
           gamma1 => gamma(1:,1:,ke); gamma2 => gamma(1:,1:,kl)

           iblank2 => iblank(1:,1:,kl)
           viscPointer => viscKmaxPointer

           if( viscous ) then
             rlv1 => rlv(1:,1:,ke); rlv2 => rlv(1:,1:,kl)
           endif

           if(equations == RANSEquations) dd2Wall => d2Wall(:,:,kl)

       end select
!
!      ******************************************************************
!      *                                                                *
!      * The actual part for storing the data. Determine the variable   *
!      * to be written and loop over the boundary faces of the subface. *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the variable to be written.

       varName: select case (solName)

         case (cgnsDensity)

           do j=rangeFace(2,1), rangeFace(2,2)
             do i=rangeFace(1,1), rangeFace(1,2)
               nn = nn + 1
               buffer(nn) = half*(ww1(i,j,irho) + ww2(i,j,irho))
             enddo
           enddo

         !===============================================================

         case (cgnsPressure)

           do j=rangeFace(2,1), rangeFace(2,2)
             do i=rangeFace(1,1), rangeFace(1,2)
               nn = nn + 1
               buffer(nn) = half*(pp1(i,j) + pp2(i,j))
             enddo
           enddo

         !===============================================================

         case (cgnsTemp)

           do j=rangeFace(2,1), rangeFace(2,2)
             do i=rangeFace(1,1), rangeFace(1,2)
               nn = nn + 1
               buffer(nn) = (pp1(i,j) + pp2(i,j)) &
                          / (RGas*(ww1(i,j,irho) + ww2(i,j,irho)))
             enddo
           enddo

         !===============================================================

         case (cgnsVelx)

           do j=rangeFace(2,1), rangeFace(2,2)
             do i=rangeFace(1,1), rangeFace(1,2)
               nn = nn + 1
               buffer(nn) = half*(ww1(i,j,ivx) + ww2(i,j,ivx))
             enddo
           enddo

         !===============================================================

         case (cgnsVely)

           do j=rangeFace(2,1), rangeFace(2,2)
             do i=rangeFace(1,1), rangeFace(1,2)
               nn = nn + 1
               buffer(nn) = half*(ww1(i,j,ivy) + ww2(i,j,ivy))
             enddo
           enddo

         !===============================================================

         case (cgnsVelz)

           do j=rangeFace(2,1), rangeFace(2,2)
             do i=rangeFace(1,1), rangeFace(1,2)
               nn = nn + 1
               buffer(nn) = half*(ww1(i,j,ivz) + ww2(i,j,ivz))
             enddo
           enddo

         !================================================================

         case (cgnsCp)

           ! Factor multiplying p-pInf

           fact = two/(gammaInf*pInf*MachCoef*MachCoef)

           do j=rangeFace(2,1), rangeFace(2,2)
             do i=rangeFace(1,1), rangeFace(1,2)
               nn = nn + 1
               buffer(nn) = fact*(half*(pp1(i,j) + pp2(i,j)) - pInf)
             enddo
           enddo

         !===============================================================

         case (cgnsPtotloss)

           ! First compute the total pressure of the free stream.

           call computePtot(rhoInf, uInf, zero, zero, &
                             pInf, ptotInf, 1_intType)
           ptotInf = one/ptotInf

           ! Loop over the faces and compute the total pressure loss.

           do j=rangeFace(2,1), rangeFace(2,2)
             do i=rangeFace(1,1), rangeFace(1,2)

               psurf = half*(pp1(i,j) + pp2(i,j))
               rsurf = half*(ww1(i,j,irho) + ww2(i,j,irho))
               usurf = half*(ww1(i,j,ivx)  + ww2(i,j,ivx))
               vsurf = half*(ww1(i,j,ivy)  + ww2(i,j,ivy))
               wsurf = half*(ww1(i,j,ivz)  + ww2(i,j,ivz))

               call computePtot(rsurf, usurf, vsurf, wsurf, &
                                psurf, ptot, 1_intType)

               nn = nn + 1
               buffer(nn) = one - ptot*ptotInf
             enddo
           enddo

         !===============================================================

         case (cgnsMach)

           do j=rangeFace(2,1), rangeFace(2,2)
             do i=rangeFace(1,1), rangeFace(1,2)

               psurf  = half*(pp1(i,j) + pp2(i,j))
               rsurf  = half*(ww1(i,j,irho) + ww2(i,j,irho))
               usurf  = half*(ww1(i,j,ivx)  + ww2(i,j,ivx))
               vsurf  = half*(ww1(i,j,ivy)  + ww2(i,j,ivy))
               wsurf  = half*(ww1(i,j,ivz)  + ww2(i,j,ivz))
               m2surf = rsurf*(usurf**2 + vsurf**2 + wsurf**2) &
                      / (half*(gamma1(i,j) + gamma2(i,j))*psurf)

               nn = nn + 1
               buffer(nn) = sqrt(m2surf)
             enddo
           enddo

!        ================================================================

         case (cgnsSkinFmag, cgnsYplus, &
               cgnsSkinFx, cgnsSkinFy, cgnsSkinFz)

           ! To avoid a lot of code duplication these 5 variables are
           ! treated together.

           ! Multiplication factor to obtain the skin friction from
           ! the wall shear stress.

           fact = two/(gammaInf*pInf*MachCoef*MachCoef)

           ! Loop over the given range of faces. As the viscous data is
           ! only present in the owned faces, the values of the halo's
           ! are set equal to the nearest physical face. Therefore the
           ! working indices are ii and jj.

           do j=rangeFace(2,1), rangeFace(2,2)
             if(j == rangeFace(2,1)) then
               jj = j + offVis
             else if(j == rangeFace(2,2)) then
               jj = j - offVis
             else
               jj = j
             endif

             do i=rangeFace(1,1), rangeFace(1,2)
               if(i == rangeFace(1,1)) then
                 ii = i + offVis
               else if(i == rangeFace(1,2)) then
                 ii = i - offVis
               else
                 ii = i
               endif

               ! Determine the viscous subface on which this
               ! face is located.

               mm = viscPointer(ii,jj)

               ! Store the 6 components of the viscous stress tensor
               ! a bit easier.

               tauxx = viscSubface(mm)%tau(ii,jj,1)
               tauyy = viscSubface(mm)%tau(ii,jj,2)
               tauzz = viscSubface(mm)%tau(ii,jj,3)
               tauxy = viscSubface(mm)%tau(ii,jj,4)
               tauxz = viscSubface(mm)%tau(ii,jj,5)
               tauyz = viscSubface(mm)%tau(ii,jj,6)

               ! Compute the "unit" force on this face. The unit normal
               ! is outward pointing per definition. A minus sign is
               ! present, because of the definition of the viscous
               ! stress tensor. Note that in the normal the indices i
               ! and j could be used. However this is not done.

               norm(1) = BCData(mm)%norm(ii,jj,1)
               norm(2) = BCData(mm)%norm(ii,jj,2)
               norm(3) = BCData(mm)%norm(ii,jj,3)

               fx = -(tauxx*norm(1) + tauxy*norm(2) + tauxz*norm(3))
               fy = -(tauxy*norm(1) + tauyy*norm(2) + tauyz*norm(3))
               fz = -(tauxz*norm(1) + tauyz*norm(2) + tauzz*norm(3))
 
               fn = fx*norm(1) + fy*norm(2) + fz*norm(3)

               fx = fx - fn*norm(1)
               fy = fy - fn*norm(2)
               fz = fz - fn*norm(3)

               ! Determine the variable to be stored and compute it.
               ! Note that an offset of -1 must be used in dd2Wall,
               ! because the original array, d2Wall, starts at 2.
               ! First update the counter nn.

               nn = nn + 1

               select case (solName)
                 case (cgnsSkinFmag)
                   buffer(nn) = fact*sqrt(fx*fx + fy*fy + fz*fz)

                 case (cgnsSkinFx)
                   buffer(nn) = fact*fx

                 case (cgnsSkinFy)
                   buffer(nn) = fact*fy

                 case (cgnsSkinFz)
                   buffer(nn) = fact*fz

                 case (cgnsYplus)
                   rsurf      = half*(ww1(ii,jj,irho) + ww2(ii,jj,irho))
                   musurf     = half*(rlv1(ii,jj)     + rlv2(ii,jj))
                   buffer(nn) = sqrt(rsurf*sqrt(fx*fx + fy*fy + fz*fz)) &
                              * dd2Wall(ii-1,jj-1)/musurf
               end select

             enddo
           enddo

!        ================================================================

         case (cgnsStanton)

           ! Some constants needed to compute the stanton number.

           gm1   = gammaInf - one
           a2Tot = gammaInf*pInf*(one + half*gm1*MachCoef*MachCoef) &
                  / rhoInf
           fact   = MachCoef*sqrt(gammaInf*pInf*rhoInf)/gm1

           ! Loop over the given range of faces. As the viscous data is
           ! only present in the owned faces, the values of the halo's
           ! are set equal to the nearest physical face. Therefore the
           ! working indices are ii and jj.

           do j=rangeFace(2,1), rangeFace(2,2)
             if(j == rangeFace(2,1)) then
               jj = j + offVis
             else if(j == rangeFace(2,2)) then
               jj = j - offVis
             else
               jj = j
             endif

             do i=rangeFace(1,1), rangeFace(1,2)
               if(i == rangeFace(1,1)) then
                 ii = i + offVis
               else if(i == rangeFace(1,2)) then
                 ii = i - offVis
               else
                 ii = i
               endif

               ! Determine the viscous subface on which this
               ! face is located.

               mm = viscPointer(ii,jj)

               ! Compute the heat flux. Multipy with the sign of the
               ! normal to obtain the correct value.

               qw = viscSubface(mm)%q(ii,jj,1)*BCData(mm)%norm(ii,jj,1) &
                  + viscSubface(mm)%q(ii,jj,2)*BCData(mm)%norm(ii,jj,2) &
                  + viscSubface(mm)%q(ii,jj,3)*BCData(mm)%norm(ii,jj,3)

               ! Compute the speed of sound squared at the wall and
               ! the stanton number, which is stored in buffer.

               a2 = half*(gamma1(ii,jj)   + gamma2(ii,jj)) &
                  *      (pp1(ii,jj)      + pp2(ii,jj))    &
                  /      (ww1(ii,jj,irho) + ww2(ii,jj,irho))

               nn = nn + 1
               buffer(nn) = qw/(fact*(a2Tot-a2))

             enddo
           enddo

!        ================================================================

         case (cgnsBlank)

           ! Loop over the given range of faces. Since iblanks are set
           ! to 2 for boundary conditions and >= 10 for the boundary,
           ! take the minimum of the value and 1, so that cells with
           ! valid data always have an iblank of 1.

           do j=rangeFace(2,1), rangeFace(2,2)
             do i=rangeFace(1,1), rangeFace(1,2)
               nn = nn + 1
               buffer(nn) = real(min(iblank2(i,j), 1_intType), realType)
             enddo
           enddo

       end select varName

       end subroutine storeSurfsolInBuffer
