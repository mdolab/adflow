!
!      ******************************************************************
!      *                                                                *
!      * File:          normalVelocities.f90                            *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-23-2004                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine normalVelocitiesAllLevels(sps)
!
!      ******************************************************************
!      *                                                                *
!      * normalVelocitiesAllLevels computes the normal grid             *
!      * velocities of some boundary faces of the moving blocks for     *
!      * spectral mode sps. All grid levels from ground level to the    *
!      * coarsest level are considered.                                 *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use iteration
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: sps
!
!      Local variables.
!
       integer(kind=intType) :: nLevels, level, nn, mm
       integer(kind=intType) :: i, j

       real(kind=realType) :: weight, mult

       real(kind=realType), dimension(:,:),   pointer :: sFace
       real(kind=realType), dimension(:,:,:), pointer :: ss
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of grid levels, starting at groundLevel,
       ! the currently finest mesh.

       nLevels = ubound(flowDoms,2)
       levelLoop: do level=groundLevel,nLevels

         ! Loop over the number of local blocks.

         domains: do nn=1,nDom

           ! Set the pointers for this block.

           call setPointers(nn, level, sps)

           ! Check for a moving block. As it is possible that in a
           ! multidisicplinary environment additional grid velocities
           ! are set, the test should be done on addGridVelocities
           ! and not on blockIsMoving.

           testMoving: if( addGridVelocities ) then
!
!            ************************************************************
!            *                                                          *
!            * Determine the normal grid velocities of the boundaries.  *
!            * As these values are based on the unit normal. A division *
!            * by the length of the normal is needed.                   *
!            * Furthermore the boundary unit normals are per definition *
!            * outward pointing, while on the iMin, jMin and kMin       *
!            * boundaries the face normals are inward pointing. This    *
!            * is taken into account by the factor mult.                *
!            *                                                          *
!            ************************************************************
!
             ! Loop over the boundary subfaces.

             bocoLoop: do mm=1,nBocos

               ! Check whether rFace is allocated.

               testAssoc: if( associated(BCData(mm)%rFace) ) then

                 ! Determine the block face on which the subface is
                 ! located and set some variables accordingly.

                 select case (BCFaceID(mm))

                   case (iMin)
                     mult = -one
                     ss => si(1,:,:,:);  sFace => sFaceI(1,:,:)
                   case (iMax)
                     mult = one
                     ss => si(il,:,:,:); sFace => sFaceI(il,:,:)
                   case (jMin)
                     mult = -one
                     ss => sj(:,1,:,:);  sFace => sFaceJ(:,1,:)
                   case (jMax)
                     mult = one
                     ss => sj(:,jl,:,:); sFace => sFaceJ(:,jl,:)
                   case (kMin)
                     mult = -one
                     ss => sk(:,:,1,:);  sFace => sFaceK(:,:,1)
                   case (kMax)
                     mult = one
                     ss => sk(:,:,kl,:); sFace => sFaceK(:,:,kl)

                 end select

                 ! Loop over the faces of the subface.

                 do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                   do i=BCData(mm)%icBeg, BCData(mm)%icEnd

                     ! Compute the inverse of the length of the normal
                     ! vector and possibly correct for inward pointing.

                     weight = sqrt(ss(i,j,1)**2 + ss(i,j,2)**2 &
                            +      ss(i,j,3)**2)
                     if(weight > zero) weight = mult/weight

                     ! Compute the normal velocity based on the outward
                     ! pointing unit normal.

                     BCData(mm)%rFace(i,j) = weight*sFace(i,j)

                   enddo
                 enddo
 
               endif testAssoc
             enddo bocoLoop

           else testMoving

             ! Block is not moving. Loop over the boundary faces and set
             ! the normal grid velocity to zero if allocated.

             do mm=1,nBocos
               if( associated(BCData(mm)%rFace) ) &
                 BCData(mm)%rFace = zero
             enddo

           endif testMoving
         enddo domains

       enddo levelLoop

       end subroutine normalVelocitiesAllLevels
