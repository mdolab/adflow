!
!      ******************************************************************
!      *                                                                *
!      * File:          hqpSetInducedVelocities.f90                     *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 11-02-2004                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine hqpSetInducedVelocities
!
!      ******************************************************************
!      *                                                                *
!      * hqpSetInducedVelocities sets or adds the additional grid       *
!      * velocities to model the wake of the helicopter blade.          *
!      * This is done for the current ground level of the computation.  *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use inputTimeSpectral
       use iteration
       use monitor
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: i, j, k, sps, nn, mm, ll
       integer(kind=intType) :: level, nLevels

       real(kind=realType) :: tCurrent

       real(kind=realType), dimension(3) :: vf

       real(kind=realType), dimension(:,:,:,:), allocatable :: coor, vel
       real(kind=realType), dimension(:,:),     allocatable :: coorFace
       real(kind=realType), dimension(:,:),     allocatable :: velFace

       real(kind=realType), dimension(:,:,:), pointer :: uSlip, xFace
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Store the current time in tCurrent and determine the total
       ! number of grid levels stored.

       tCurrent = timeUnsteady + timeUnsteadyRestart
       nLevels  = ubound(flowDoms,2)

       ! Loop over the number of spectral solutions and local blocks.

       spectralLoop: do sps=1,nTimeIntervalsSpectral
         domains: do nn=1,nDom

           ! Set addGridVelocities to .True. For this block and all
           ! its coarser blocks. Note that this must be done for the
           ! flowDoms array, because addGridVelocities is not a pointer.

           do level=groundLevel,nLevels
             flowDoms(nn,level,sps)%addGridVelocities = .true.
           enddo

           ! Have the pointers point to this block.

           call setPointers(nn,groundLevel,sps)

           ! Allocate the memory for coor and vel.

           allocate(coor(3,0:ie,0:je,0:ke), vel(3,0:ie,0:je,0:ke), &
                    stat=ierr)
           if(ierr /= 0)                               &
             call terminate("hqpSetInducedVelocities", &
                            "Memory allocation failure for coor and vel")

           ! Copy the coordinates into coor.

           do k=0,ke
             do j=0,je
               do i=0,ie
                 coor(1,i,j,k) = x(i,j,k,1)
                 coor(2,i,j,k) = x(i,j,k,2)
                 coor(3,i,j,k) = x(i,j,k,3)
               enddo
             enddo
           enddo

           ! Determine the induced velocities.

           ll = (ke+1)*(je+1)*(ie+1)
           call hqpPert(coor, vel, ll, tCurrent)

           ! Check if the block is flagged as not moving. In that
           ! case some additional things need to be done.

           if(.not. blockIsMoving) then

             ! Check if the memory for the mesh velocities is allocated.
             ! If not, do so, including the coarser grids. As these
             ! variables are pointers, associated should be used.

             if(.not. associated(s) ) &
               call hqpAllocMeshVelocities(nLevels)

             ! Initialize the mesh velocities to zero.

             s      = zero
             sFaceI = zero
             sFaceJ = zero
             sFaceK = zero

           endif

           ! Loop over the cell centers and update the grid velocities
           ! of the cell centers, including the 1st level halo cells.

           do ll=1,3
             do k=1,ke
               do j=1,je
                 do i=1,ie
                   s(i,j,k,ll) = s(i,j,k,ll) + eighth &
                               * (vel(ll,i-1,j-1,k-1) + vel(ll,i,j-1,k-1) &
                               +  vel(ll,i-1,j,  k-1) + vel(ll,i,j,  k-1) &
                               +  vel(ll,i-1,j-1,k  ) + vel(ll,i,j-1,k  ) &
                               +  vel(ll,i-1,j,  k  ) + vel(ll,i,j,  k  ))
                 enddo
               enddo
             enddo
           enddo

           ! Normal velocities in i-direction.

           do k=1,ke
             do j=1,je
               do i=0,ie

                 ! Compute four times the velocity in the face center.

                 vf(1) = vel(1,i,j-1,k-1) + vel(1,i,j,k-1) &
                       + vel(1,i,j-1,k  ) + vel(1,i,j,k  )
                 vf(2) = vel(2,i,j-1,k-1) + vel(2,i,j,k-1) &
                       + vel(2,i,j-1,k  ) + vel(2,i,j,k  )
                 vf(3) = vel(3,i,j-1,k-1) + vel(3,i,j,k-1) &
                       + vel(3,i,j-1,k  ) + vel(3,i,j,k  )

                 ! Update the normal velocity.

                 sFaceI(i,j,k) = sFaceI(i,j,k)             &
                               + fourth*(vf(1)*si(i,j,k,1) &
                               +         vf(2)*si(i,j,k,2) &
                               +         vf(3)*si(i,j,k,3))
               enddo
             enddo
           enddo

           ! Normal velocities in j-direction.

           do k=1,ke
             do j=0,je
               do i=1,ie

                 ! Compute four times the velocity in the face center.

                 vf(1) = vel(1,i-1,j,k-1) + vel(1,i,j,k-1) &
                       + vel(1,i-1,j,k  ) + vel(1,i,j,k  )
                 vf(2) = vel(2,i-1,j,k-1) + vel(2,i,j,k-1) &
                       + vel(2,i-1,j,k  ) + vel(2,i,j,k  )
                 vf(3) = vel(3,i-1,j,k-1) + vel(3,i,j,k-1) &
                       + vel(3,i-1,j,k  ) + vel(3,i,j,k  )

                 ! Update the normal velocity.

                 sFaceJ(i,j,k) = sFaceJ(i,j,k)             &
                               + fourth*(vf(1)*sj(i,j,k,1) &
                               +         vf(2)*sj(i,j,k,2) &
                               +         vf(3)*sj(i,j,k,3))
               enddo
             enddo
           enddo

           ! Normal velocities in k-direction.

           do k=0,ke
             do j=1,je
               do i=1,ie

                 ! Compute four times the velocity in the face center.

                 vf(1) = vel(1,i-1,j-1,k) + vel(1,i,j-1,k) &
                       + vel(1,i-1,j,  k) + vel(1,i,j,  k)
                 vf(2) = vel(2,i-1,j-1,k) + vel(2,i,j-1,k) &
                       + vel(2,i-1,j,  k) + vel(2,i,j,  k)
                 vf(3) = vel(3,i-1,j-1,k) + vel(3,i,j-1,k) &
                       + vel(3,i-1,j,  k) + vel(3,i,j,  k)

                 ! Update the normal velocity.

                 sFaceK(i,j,k) = sFaceK(i,j,k)             &
                               + fourth*(vf(1)*sk(i,j,k,1) &
                               +         vf(2)*sk(i,j,k,2) &
                               +         vf(3)*sk(i,j,k,3))
               enddo
             enddo
           enddo

           ! Release the memory of the coordinates and induced
           ! velocities.

           deallocate(coor, vel, stat=ierr)
           if(ierr /= 0)                               &
             call terminate("hqpSetInducedVelocities", &
                            "Deallocation failure for coor and vel")

           ! Loop over the viscous subfaces of this block.

           bocoLoop: do mm=1,nViscBocos

             ! Determine the grid face on which the subface is located
             ! and set the pointer for the coordinates accordingly.

             select case (BCFaceID(mm))
               case (iMin)
                 xFace => x(1, :,:,:)
               case (iMax)
                 xFace => x(il,:,:,:)
               case (jMin)
                 xFace => x(:,1, :,:)
               case (jMax)
                 xFace => x(:,jl,:,:)
               case (kMin)
                 xFace => x(:,:,1, :)
               case (kMax)
                 xFace => x(:,:,kl,:)
             end select

             ! Allocate the memory for the coordinates and induced
             ! velocities of the face center.

             ll = (BCData(mm)%jcEnd - BCData(mm)%jcBeg + 1) &
                * (BCData(mm)%icEnd - BCData(mm)%icBeg + 1)
             allocate(coorFace(3,ll), velFace(3,ll), stat=ierr)
             if(ierr /= 0)                                  &
               call terminate("hqpSetInducedVelocities", &
                              "Memory allocation failure for coorFace &
                              &and velFace")

             ! Loop over the face centers and determine their
             ! coordinates. Note that due to the usage of the pointer
             ! xFace an offset of +1 must be used in the coordinate
             ! arrays, because x starts at 0 for i, j and k.

             ll = 0
             do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
               do i=BCData(mm)%icBeg, BCData(mm)%icEnd

                 ll = ll + 1
                 coorFace(1,ll) = fourth*(xFace(i+1,j+1,1) &
                                +         xFace(i+1,j,  1) &
                                +         xFace(i,  j+1,1) &
                                +         xFace(i,  j,  1))
                 coorFace(2,ll) = fourth*(xFace(i+1,j+1,2) &
                                +         xFace(i+1,j,  2) &
                                +         xFace(i,  j+1,2) &
                                +         xFace(i,  j,  2))
                 coorFace(3,ll) = fourth*(xFace(i+1,j+1,3) &
                                +         xFace(i+1,j,  3) &
                                +         xFace(i,  j+1,3) &
                                +         xFace(i,  j,  3))
               enddo
             enddo

             ! Determine the induced velocities.

             call hqpPert(coorFace, velFace, ll, tCurrent)

             ! Loop again over the face centers and add the
             ! induced velocities to the slip velocities.

             uSlip => BCData(mm)%uSlip
             ll = 0

             do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
               do i=BCData(mm)%icBeg, BCData(mm)%icEnd

                 ll = ll + 1
                 uSlip(i,j,1) = uSlip(i,j,1) + velFace(1,ll)
                 uSlip(i,j,2) = uSlip(i,j,2) + velFace(2,ll)
                 uSlip(i,j,3) = uSlip(i,j,3) + velFace(3,ll)

               enddo
             enddo

             ! Deallocate the memory for coorFace and velFace.

             deallocate(coorFace, velFace, stat=ierr)
             if(ierr /= 0)                               &
               call terminate("hqpSetInducedVelocities", &
                              "Deallocation failure for coorFace &
                              &and velFace")

           enddo bocoLoop

         enddo domains
       enddo spectralLoop

       end subroutine hqpSetInducedVelocities

!      ==================================================================

       subroutine hqpAllocMeshVelocities(nLevels)
!
!      ******************************************************************
!      *                                                                *
!      * HqpAllocMeshVelocities allocates the mesh velocities for       *
!      * the currently active block as well as its coarser blocks.      *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use iteration
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nLevels
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: sps, nn, level
       integer(kind=intType) :: iie, jje, kke
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Store the current block id and spectral solution a bit easier.

       nn  = nbkLocal
       sps = spectralSol

       ! Loop over the number of multigrid levels.

       levelLoop: do level=groundLevel,nLevels

         ! Check if the mesh velocities are not allocated. It is
         ! possible that the coarser levels are allocated due to
         ! grid sequencing.

         if(.not. associated(flowDoms(nn,level,sps)%s) ) then

           ! Easier storage of the array boundaries.

           iie = flowDoms(nn,level,sps)%ie
           jje = flowDoms(nn,level,sps)%je
           kke = flowDoms(nn,level,sps)%ke

           ! Allocate the memory.

           allocate(flowDoms(nn,level,sps)%s(iie,jje,kke,3),      &
                    flowDoms(nn,level,sps)%sFaceI(0:iie,jje,kke), &
                    flowDoms(nn,level,sps)%sFaceJ(iie,0:jje,kke), &
                    flowDoms(nn,level,sps)%sFaceK(iie,jje,0:kke), &
                    stat=ierr)
           if(ierr /= 0)                              &
             call terminate("hqpAllocMeshVelocities", &
                            "Memory allocation failure for s, &
                            &sFaceI, sFaceJ and sFaceK.")

         endif

       enddo levelLoop

       ! Reset the grid velocity pointers of blockPointers, such that
       ! they correspond to the just allocated memory on groundLevel.

       s      => flowDoms(nn,groundLevel,sps)%s
       sFaceI => flowDoms(nn,groundLevel,sps)%sFaceI
       sFaceJ => flowDoms(nn,groundLevel,sps)%sFaceJ
       sFaceK => flowDoms(nn,groundLevel,sps)%sFaceK

       end subroutine hqpAllocMeshVelocities
