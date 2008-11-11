!
!      ******************************************************************
!      *                                                                *
!      * File:          normalVelocitiesforcesAdj.f90                   *
!      * Author:        Edwin van der Weide, C.A.(Sandy) Mader          *
!      * Starting date: 02-23-2004                                      *
!      * Last modified: 10-25-2008                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine normalVelocitiesAllLevelsforcesAdj(sps,mm,sFaceIAdj,&
            iiBeg,iiEnd,jjBeg,jjEnd,i2Beg,i2End,j2Beg,j2End,&
            sFaceJAdj,sFaceKAdj,siAdj, sjAdj, skAdj,rFaceAdj)
!(sps)
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
       integer(kind=intType), intent(in) :: sps,mm

       real(kind=realType), dimension(0:2,iiBeg:iiEnd,jjBeg:jjEnd,3) :: siAdj 
       ! notice the range of y dim is set 1:2 which corresponds to 1/jl
       real(kind=realType), dimension(iiBeg:iiEnd,0:2,jjBeg:jjEnd,3) :: sjAdj
       ! notice the range of z dim is set 1:2 which corresponds to 1/kl
       real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd,0:2,3) :: skAdj

       real(kind=realType), dimension(1:2,iiBeg:iiEnd,jjBeg:jjEnd) :: sFaceiAdj
       real(kind=realType), dimension(iiBeg:iiEnd,1:2,jjBeg:jjEnd) :: sFacejAdj
       real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd,1:2) :: sFacekAdj
       real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd) :: rFaceAdj

       integer(kind=intType), intent(in) :: iiBeg,iiEnd,jjBeg,jjEnd
       integer(kind=intType), intent(in) :: i2Beg,i2End,j2Beg,j2End
!
!      Local variables.
!
       integer(kind=intType) :: nLevels, level, nn!, mm
       integer(kind=intType) :: i, j

       real(kind=realType) :: weight, mult

       real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd,3) ::ssAdj
       real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd) :: sfaceAdj
       !real(kind=realType), dimension(:,:),   pointer :: sFace
       !real(kind=realType), dimension(:,:,:), pointer :: ss
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
!moved outside this function
!!$       ! Loop over the number of grid levels, starting at groundLevel,
!!$       ! the currently finest mesh.
!!$
!!$       nLevels = ubound(flowDoms,2)
!!$       levelLoop: do level=groundLevel,nLevels
!!$
!!$         ! Loop over the number of local blocks.
!!$
!!$         domains: do nn=1,nDom
!!$
!!$           ! Set the pointers for this block.
!!$
!!$           call setPointers(nn, level, sps)

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

!also moved outside....
!!$             ! Loop over the boundary subfaces.
!!$
!!$             bocoLoop: do mm=1,nBocos
!!$
!!$               ! Check whether rFace is allocated.
!!$
!!$               testAssoc: if( associated(BCData(mm)%rFace) ) then

                 ! Determine the block face on which the subface is
                 ! located and set some variables accordingly.

                 select case (BCFaceID(mm))

                   case (iMin)
                     mult = -one
                     ssAdj = siAdj(1,:,:,:);  sFaceAdj = sFaceIAdj(1,:,:)
                   case (iMax)
                     mult = one
                     ssAdj = siAdj(2,:,:,:)! which was si(il,:,:,:)
                     sFaceAdj = sFaceIAdj(2,:,:)! which was sFaceI(il,:,:)
                   case (jMin)
                     mult = -one
                     ssAdj = sjAdj(:,1,:,:);  sFaceAdj = sFaceJAdj(:,1,:)
                   case (jMax)
                     mult = one
                     ssAdj = sjAdj(:,2,:,:)! which was sj(:,jl,:,:)
                     sFaceAdj = sFaceJAdj(:,2,:)! which was sFaceJ(:,jl,:)
                   case (kMin)
                     mult = -one
                     ssAdj = skAdj(:,:,1,:);  sFaceAdj = sFaceKAdj(:,:,1)
                   case (kMax)
                     mult = one
                     ssAdj = skAdj(:,:,2,:)! which was sk(:,:,kl,:)
                     sFaceAdj = sFaceKAdj(:,:,2)! which was sFacek(:,:,kl)

                 end select

                 ! Loop over the faces of the subface.
                 do j=jjBeg,jjEnd
                    do i=iiBeg,iiEnd
                 !do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                 !  do i=BCData(mm)%icBeg, BCData(mm)%icEnd

                     ! Compute the inverse of the length of the normal
                     ! vector and possibly correct for inward pointing.

                     weight = sqrt(ssAdj(i,j,1)**2 + ssAdj(i,j,2)**2 &
                            +      ssAdj(i,j,3)**2)
                     if(weight > zero) weight = mult/weight

                     ! Compute the normal velocity based on the outward
                     ! pointing unit normal.

                     !BCData(mm)%rFace(i,j) = weight*sFace(i,j)
                     rFaceAdj(i,j) = weight*sFaceAdj(i,j)

                   enddo
                 enddo
 
!               endif testAssoc
!             enddo bocoLoop

           else testMoving

             ! Block is not moving. Loop over the boundary faces and set
             ! the normal grid velocity to zero if allocated.

             !do mm=1,nBocos
               !if( associated(BCData(mm)%rFace) ) &
               !  BCData(mm)%rFace = zero
             rFaceAdj = zero
             !enddo

           endif testMoving
!         enddo domains

!       enddo levelLoop

         end subroutine normalVelocitiesAllLevelsforcesAdj
