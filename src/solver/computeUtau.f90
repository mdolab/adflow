!
!      ******************************************************************
!      *                                                                *
!      * File:          computeUtau.f90                                 *
!      * Author:        Georgi Kalitzin, Edwin van der Weide            *
!      * Starting date: 03-03-2004                                      *
!      * Last modified: 03-25-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine computeUtau
!
!      ******************************************************************
!      *                                                                *
!      * computeUtau computes the skin friction velocity for the        *
!      * viscous subfaces. This data is only needed if wall functions   *
!      * are used.                                                      *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use inputPhysics
       use inputTimeSpectral
       use iteration
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: sps, nn, mm, i, j

       real(kind=realType) :: re, vx, vy, vz, veln, veltmag

       real(kind=realType), dimension(:,:,:), pointer :: ww, norm, uSlip
       real(kind=realType), dimension(:,:),   pointer :: dd2Wall, rrlv
       real(kind=realType), dimension(:,:),   pointer :: utau
!
!      Function definition.
!
       real(kind=realType) :: curveUpRe
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if no wall functions must be used.

       if(.not. wallFunctions) return

       ! Loop over the number of spectral solutions and local blocks.

       spectralLoop: do sps=1,nTimeIntervalsSpectral
         domains: do nn=1,nDom

           ! Set the pointers to this block.

           call setPointers(nn, groundLevel, sps)

           ! Loop over the viscous subfaces of this block.

           viscSubfaces: do mm=1,nViscBocos

             ! Set a bunch of pointers depending on the face id to make
             ! a generic treatment possible.

             select case (BCFaceID(mm))

               case (iMin)
                 ww => w(2,1:,1:,:);
                 dd2Wall => d2Wall(2,:,:); rrlv => rlv(2,1:,1:)

               !=========================================================

               case (iMax)
                 ww => w(il,1:,1:,:)
                 dd2Wall => d2Wall(il,:,:); rrlv => rlv(il,1:,1:)

               !=========================================================

               case (jMin)
                 ww => w(1:,2,1:,:)
                 dd2Wall => d2Wall(:,2,:); rrlv => rlv(1:,2,1:)

               !=========================================================

               case (jMax)
                 ww => w(1:,jl,1:,:)
                 dd2Wall => d2Wall(:,jl,:); rrlv => rlv(1:,jl,1:)

               !=========================================================

               case (kMin)
                 ww => w(1:,1:,2,:)
                 dd2Wall => d2Wall(:,:,2); rrlv => rlv(1:,1:,2)

               !=========================================================

               case (kMax)
                 ww => w(1:,1:,kl,:)
                 dd2Wall => d2Wall(:,:,kl); rrlv => rlv(1:,1:,kl)

             end select

             ! Set the pointers for the unit outward normals, uSlip
             ! and utau to make the code more readable.

             norm  => BCData(mm)%norm
             uSlip => BCData(mm)%uSlip
             utau  => viscSubface(mm)%utau

             ! Loop over the quadrilateral faces of the subface. Note
             ! that the nodal range of BCData must be used and not the
             ! cell range, because the latter may include the halo's in i
             ! and j-direction. The offset +1 is there, because inBeg and
             ! jnBeg refer to nodal ranges and not to cell ranges.
             ! Note that an offset of -1 must be used in dd2Wall, because
             ! the original array d2Wall starts at 2.

             do j=(BCData(mm)%jnBeg+1),BCData(mm)%jnEnd
               do i=(BCData(mm)%inBeg+1),BCData(mm)%inEnd

                 ! Compute the velocity difference between the internal
                 ! cell and the wall.

                 vx = ww(i,j,ivx) - uSlip(i,j,1)
                 vy = ww(i,j,ivy) - uSlip(i,j,2)
                 vz = ww(i,j,ivz) - uSlip(i,j,3)

                 ! Compute the normal velocity of the internal cell.

                 veln  = vx*norm(i,j,1) + vy*norm(i,j,2) + vz*norm(i,j,3)

                 ! Compute the magnitude of the tangential velocity.

                 veltmag = max(eps,sqrt(vx*vx + vy*vy + vz*vz - veln*veln))

                 ! Compute the Reynolds number. Note that an offset of -1
                 ! must be used in dd2Wall, because the original array
                 ! d2Wall starts at 2.
                 ! Afterwards compute utau.

                 re = ww(i,j,irho)*veltmag*dd2Wall(i-1,j-1)/rrlv(i,j)
                 utau(i,j) = veltmag/max(curveUpRe(re),eps)

               enddo
             enddo

           enddo viscSubfaces
         enddo domains
       enddo spectralLoop

       end subroutine computeUtau
