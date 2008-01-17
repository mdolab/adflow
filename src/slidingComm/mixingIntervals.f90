!
!      ******************************************************************
!      *                                                                *
!      * File:          mixingIntervals.f90                             *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-31-2005                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine mixingIntervals(slideID, level, color, nSlices)
!
!      ******************************************************************
!      *                                                                *
!      * MixingIntervals determines the point distribution used as      *
!      * interpolation intervals on the given grid level for the given  *
!      * sliding interface.                                             *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use interfaceGroups
       use localSubfacesMod
       use mixingData
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: slideID
       integer(kind=intType), intent(in) :: level
       integer(kind=intType), intent(in) :: color
       integer(kind=intType), intent(in) :: nSlices
!
!      Local variables.
!
       integer :: comm, nProc, myID, ierr, nl

       integer, dimension(myInterfaces(color)%nProcSlide) :: recvcounts
       integer, dimension(myInterfaces(color)%nProcSlide) :: displs

       integer(kind=intType) :: nn, mm, kk, jj, ii, i, j
       integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
       integer(kind=intType) :: nSingularQuad, nnzq, nLocal, nAlloc
       integer(kind=intType) :: nAngles, nCuts, nGlobal

       integer(kind=intType), dimension(4) :: nk, edgeCut
       integer(kind=intType), dimension(5) :: quadrants

       integer(kind=intType), dimension(:), allocatable :: intBuf1
       integer(kind=intType), dimension(:), allocatable :: intBuf2

       real(kind=realType) :: rMin, rMax, axMin, axMax, ttMin
       real(kind=realType) :: rMinLocal, ttMinLocal, r, tt
       real(kind=realType) :: tMinLocal, t, psi, perAngle
       real(kind=realType) :: thetapMin, thetapMax
       real(kind=realType) :: thetanMin, thetanMax
       real(kind=realType) :: ww1, ww2, dr, theta1, theta2

       real(kind=realType), dimension(4) :: theta, rad, radCut

       real(kind=realType), dimension(:),     allocatable :: angles
       real(kind=realType), dimension(:,:),   allocatable :: realBuf1
       real(kind=realType), dimension(:,:),   allocatable :: realBuf2
       real(kind=realType), dimension(:,:,:), allocatable :: xFace

       logical :: addMinHalo, thetaPiCrossed

       type(mixingIntervalType), dimension(:), allocatable :: edges
       type(mixingIntervalType), dimension(:), pointer     :: intervals
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Store the communicator, the number of processors and my
       ! processor ID for this processor group a bit easier.

       comm  = myInterfaces(color)%commSlide
       nProc = myInterfaces(color)%nProcSlide
       myID  = myInterfaces(color)%myIDSlide

       ! Determine the periodic angle for the interface and initialize
       ! the minimum and maximum positive and negative angles as well
       ! as thetaPiCrossed. The latter logical indicates whether or not
       ! the sliding interface crosses the line theta == pi.

       perAngle = two*pi/nSlices

       thetapMin =  pi
       thetapMax = -pi
       thetanMin =  pi
       thetanMax = -pi

       thetaPiCrossed = .false.

       ! Determine the type of interface we are having here. Either
       ! a radial or an axial interface. The decision is based on
       ! geometrical reasons.

       rMin  = myInterfaces(color)%rMin
       rMax  = myInterfaces(color)%rMax
       axMin = myInterfaces(color)%axMin
       axMax = myInterfaces(color)%axMax

       if(abs(rMax-rMin) >= abs(axMax-axMin)) then
         radialInterface = .true.
       else
         radialInterface = .false.
       endif
!
!      ******************************************************************
!      *                                                                *
!      * Determine the local candidate of the starting point of the     *
!      * interpolation interval. If a polar singular quad is present    *
!      * this is automatically the starting point. Otherwise a point on *
!      * the minimum boundary will be chosen. Note that the algorithm   *
!      * below is valid for both radial and axial interface.            *
!      *                                                                *
!      ******************************************************************
!
       ! Some initializations.

       nSingularQuad = 0
       rMinLocal     = large

       ! Loop over the local subfaces of the sliding interface and
       ! determine a possible local starting position for the
       ! interpolation interval.

       domainLoop1: do nn=1,nDom

         ! Set the pointers for this block. Only the 1st spectral
         ! level needs to be considered.

         call setPointers(nn, level, 1_intType)

         ! Loop over the boundary conditions.

         bocoLoop1: do mm=1,nBocos

           ! Check if the subface belongs to the sliding mesh part
           ! currently investigated.

           testInterface1: if(BCType(mm)   == slidingInterface .and. &
                              groupNum(mm) == slideID) then

             ! Store the nodal range of the subface a bit easier.

             jBeg = BCData(mm)%jnBeg
             jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg
             iEnd = BCData(mm)%inEnd

             ! Allocate the memory for the nodal cylindrical coordinates
             ! of the subface and determine them.

             allocate(xFace(iBeg:iEnd,jBeg:jEnd,3), stat=ierr)
             if(ierr /= 0) &
               call terminate("mixingIntervals", &
                              "Memory allocation failure for xFace.")

             call cylCoorNodesOnBlockFace(xFace, iBeg, iEnd,  &
                                          jBeg,  jEnd, color, &
                                          BCFaceID(mm))

             ! If this is a radial interface, store the radii in the
             ! first position such that the same piece of code can be
             ! used for both radial and axial interfaces.

             if( radialInterface ) then
               do j=jBeg,jEnd
                 do i=iBeg,iEnd
                   xFace(i,j,1) = xFace(i,j,2)
                 enddo
               enddo
             endif

             ! Loop over the nodes and update the minimum and maximum
             ! positive and negative angles.

             do j=jBeg,jEnd
               do i=iBeg,iEnd

                 if(xFace(i,j,3) >= zero) then
                   thetapMax = max(thetapMax, xFace(i,j,3))
                   thetapMin = min(thetapMin, xFace(i,j,3))
                 else
                   thetanMax = max(thetanMax, xFace(i,j,3))
                   thetanMin = min(thetanMin, xFace(i,j,3))
                 endif

               enddo
             enddo

             ! Check if a polar singular quad is present or if the lines
             ! theta == zero or theta == pi are crossed. In the latter
             ! case the minimum and maximum positive and negative angles
             ! need to be adapted.

             do j=(jBeg+1),jEnd
               do i=(iBeg+1),iEnd

                 ! Store the 4 polar angles of the quad a bit easier.

                 theta(1) = xFace(i-1,j-1,3)
                 theta(2) = xFace(i,  j-1,3)
                 theta(3) = xFace(i,  j  ,3)
                 theta(4) = xFace(i-1,j,  3)

                 ! Determine the number of nodes in the 4 quadrants of
                 ! the polar plane.

                 nk(1) = 0; nk(2) = 0; nk(3) = 0; nk(4) = 0

                 do jj=1,4
                   if(theta(jj) > zero) then
                     if(theta(jj) > half*pi) then
                       nk(2) = nk(2) + 1
                     else
                       nk(1) = nk(1) + 1
                     endif
                   else
                     if(theta(jj) > -half*pi) then
                       nk(4) = nk(4) + 1
                     else
                       nk(3) = nk(3) + 1
                     endif
                   endif
                 enddo

                 ! Determine the number of nonzero quadrants.

                 nnzq = 0
                 if(nk(1) > 0) nnzq = nnzq + 1
                 if(nk(2) > 0) nnzq = nnzq + 1
                 if(nk(3) > 0) nnzq = nnzq + 1
                 if(nk(4) > 0) nnzq = nnzq + 1

                 ! Check if this is a polar singular quad.
                 ! If so update nSingularQuad

                 select case (nnzq)
                   case (4_intType, 3_intType)

                     ! Rotation center is inside the quad, i.e. it
                     ! is polar singular.

                     nSingularQuad = nSingularQuad + 1

                   !=====================================================

                   case (2_intType)

                     ! Quad has nodes in two quadrants. Investigate the
                     ! situation a bit further.

                     if(nk(2) > 0 .and. nk(3) > 0) then

                       ! The quad crosses the line theta == pi.
                       ! Change thetapMax and thetanMin and set
                       ! thetaPiCrossed to .true.

                       thetapMax =  pi
                       thetanMin = -pi

                       thetaPiCrossed = .true.

                     else if(nk(1) > 0 .and. nk(4) > 0) then

                       ! The quad crosses the line theta == 0.
                       ! Change thetapMin and thetanMax.

                       thetapMin = zero
                       thetapMax = zero

                     else if((nk(1) > 0 .and. nk(3) > 0) .or. &
                             (nk(2) > 0 .and. nk(4) > 0)) then

                       ! Quad has nodes in opposite quadrants, which
                       ! means that this is a polar singular quad.

                       nSingularQuad = nSingularQuad + 1

                     endif

                 end select

               enddo
             enddo

             ! If no singular quads are present, check the 4 edges
             ! of the subface for the minimum coordinate.

             testNoSingularQuad: if(nSingularQuad == 0) then

               ! iMin boundary. To make sure that the point is unique
               ! use the edge center rather than the nodes. Nodes can
               ! be present on multiple processors. In the case of equal
               ! radii, the one is taken whose angle is closest to zero.
               ! This is just a way to make sure that a unique starting
               ! point is obtained in parallel mode when blocks could
               ! have been split.

               do j=(jBeg+1),jEnd
                 r  = half*(xFace(iBeg,j-1,1) + xFace(iBeg,j,1))
                 t  = half*(xFace(iBeg,j-1,3) + xFace(iBeg,j,3))
                 tt = abs(t)

                 if(r < rMinLocal) then
                   rMinLocal = r; tMinLocal = t; ttMinLocal = tt
                 else if(r == rMinLocal .and. tt < ttMinLocal) then
                   rMinLocal = r; tMinLocal = t; ttMinLocal = tt
                 endif
               enddo

               ! iMax boundary.

               do j=(jBeg+1),jEnd
                 r  = half*(xFace(iEnd,j-1,1) + xFace(iEnd,j,1))
                 t  = half*(xFace(iEnd,j-1,3) + xFace(iEnd,j,3))
                 tt = abs(t)

                 if(r < rMinLocal) then
                   rMinLocal = r; tMinLocal = t; ttMinLocal = tt
                 else if(r == rMinLocal .and. tt < ttMinLocal) then
                   rMinLocal = r; tMinLocal = t; ttMinLocal = tt
                 endif
               enddo

               ! jMin boundary.

               do i=(iBeg+1),iEnd
                 r  = half*(xFace(i-1,jBeg,1) + xFace(i,jBeg,1))
                 t  = half*(xFace(i-1,jBeg,3) + xFace(i,jBeg,3))
                 tt = abs(t)

                 if(r < rMinLocal) then
                   rMinLocal = r; tMinLocal = t; ttMinLocal = tt
                 else if(r == rMinLocal .and. tt < ttMinLocal) then
                   rMinLocal = r; tMinLocal = t; ttMinLocal = tt
                 endif
               enddo

               ! jMax boundary.

               do i=(iBeg+1),iEnd
                 r  = half*(xFace(i-1,jEnd,1) + xFace(i,jEnd,1))
                 t  = half*(xFace(i-1,jEnd,3) + xFace(i,jEnd,3))
                 tt = abs(t)

                 if(r < rMinLocal) then
                   rMinLocal = r; tMinLocal = t; ttMinLocal = tt
                 else if(r == rMinLocal .and. tt < ttMinLocal) then
                   rMinLocal = r; tMinLocal = t; ttMinLocal = tt
                 endif
               enddo

             endif testNoSingularQuad

             ! Release the memory of the nodal coordinates again.

             deallocate(xFace, stat=ierr)
             if(ierr /= 0)                       &
               call terminate("mixingIntervals", &
                              "Deallocation failure for xFace.")

           endif testInterface1
         enddo bocoLoop1
       enddo domainLoop1
!
!      ******************************************************************
!      *                                                                *
!      * Determine the interpolation angle. If a polar singular quad is *
!      * present it means that the entire wheel is present and thus any *
!      * angle can be taken. Here zero is chosen. If no polar singular  *
!      * quad is present the angle corresponding to the minimum radius  *
!      * it taken.                                                      *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the global number of polar singularities. If more
       ! than one is found print an error message and exit.

       nn = nSingularQuad
       call mpi_allreduce(nn, nSingularQuad, 1, sumb_integer, &
                          mpi_sum, comm, ierr)

       if(nSingularQuad > 1) then
         if(myID == 0)                       &
           call terminate("mixingIntervals", &
                          "More than one polar singular quad found &
                          &in sliding mesh interface.")
         call mpi_barrier(comm, ierr)
       endif

       ! Set addMinHalo to .false. to indicate that no halo needs to
       ! be added at the min boundary.

       addMinHalo = .false.

       ! Test for the presence of a polar singular quad.

       testPolarQuadPresent: if(nSingularQuad == 1) then

         ! Set the angle psi, the interpolation direction, to zero.

         psi = zero

       else testPolarQuadPresent

         ! No polar singular quad is present. Some more work needs to
         ! be done to find the starting processor and the interpolation
         ! angle psi.

         ! Determine the global minimum of the radius.

         call mpi_allreduce(rMinLocal, rMin, 1, sumb_real, &
                            mpi_min, comm, ierr)

         ! It is possible that the minimum occurs on multiple processors.
         ! Therefore take the point with the polar angle closest to zero.
         ! To avoid problems set the local minimum angle to a large
         ! value such that it does not interfer.

         if(rminLocal > rMin) ttMinLocal = large

         call mpi_allreduce(ttMinLocal, ttMin, 1, sumb_real, &
                            mpi_min, comm, ierr)

         ! If I contain the minimum of the angles (and thus the minimum
         ! of the radii) I'm the processor where the interval starts.
         ! Set the interpolation angle, for the moment stored in tt;
         ! on other processors set it to a value which is always larger.

         if(ttMin == ttMinLocal) then
           tt = tMinLocal
         else
           tt = two*pi
         endif

         ! Perform a global reduce such that the interpolation direction
         ! psi is known on all processors.

         call mpi_allreduce(tt, psi, 1, sumb_real, mpi_min, comm, ierr)

         ! If the minimum radius does not correspond to the rotation
         ! point (or if this is an axial interface) a halo must be added
         ! to the lower boundary later on.

         if((rMin > eps) .or. (.not. radialInterface)) &
            addMinHalo = .true.

       endif testPolarQuadPresent
!
!      ******************************************************************
!      *                                                                *
!      * Determine the local edges, i.e. determine for every            *
!      * quadrilateral of the sliding mesh subfaces whether or not it   *
!      * is cut but the line of the interpolation angle. Rotational     *
!      * periodicity is taken into account, if appropriate.             *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the global minima and maxima of the polar angles.
       ! Use theta1 as a temporary buffer.

       theta1 = thetapMin
       call mpi_allreduce(theta1, thetapMin, 1, sumb_real, &
                          mpi_min, comm, ierr)
       theta1 = thetanMin
       call mpi_allreduce(theta1, thetanMin, 1, sumb_real, &
                          mpi_min, comm, ierr)

       theta1 = thetapMax
       call mpi_allreduce(theta1, thetapMax, 1, sumb_real, &
                          mpi_max, comm, ierr)
       theta1 = thetanMax
       call mpi_allreduce(theta1, thetanMax, 1, sumb_real, &
                          mpi_max, comm, ierr)

       ! Check if the line theta == pi is crossed for this subface.

       nn = 0
       if( thetaPiCrossed ) nn = 1
       call mpi_allreduce(nn, mm, 1, sumb_integer, mpi_max, comm, ierr)

       thetaPiCrossed = .false.
       if(mm == 1) thetaPiCrossed = .true.

       ! Determine the number of angles to check for. If this a
       ! rotational periodic interface not only psi, but some of its
       ! periodic neighbors must be taken into account.

       if(nSlices == 1) then

         ! No rotational periodicity. Only one angle needs to be checked.

         nAngles = 1

       else

         ! Rotational periodicity is used. Determine the number of
         ! sections to cover the entire range of angles present.

         ! Store the angle range in theta1 to theta2.

         if( thetaPiCrossed ) then

           ! The line theta == pi is crossed. Make sure that the angle
           ! range is continuous by adding 2*pi to thetanMax and
           ! possibly to psi.

           theta1 = thetapMin
           theta2 = thetanMax + two*pi
           if(psi <= thetanMax) psi = psi + two*pi

         else

           ! Regular continuous range. The min and max function can be
           ! used because of the initializations of thetapMin,
           ! thetapMax, thetanMin and thetanMax.

           theta1 = min(thetapMin, thetanMin)
           theta2 = max(thetapMax, thetanMax)

         endif

         ! Determine the number of rotations, starting from psi, in
         ! positive and negative direction to reach theta1 and theta2.
         ! Determine the number of angles needed accordingly.

         nn = (theta1 - psi)/perAngle
         mm = (theta2 - psi)/perAngle

         nAngles = mm - nn + 1
         nAngles = min(nAngles,nSlices)

         ! Adapt psi such that it corresponds to the smallest periodic
         ! equivalent of the original value that still intersects with
         ! the sliding interface.

         psi = psi + nn*perAngle

       endif

       ! Allocate the memory for angles.

       allocate(angles(nAngles), stat=ierr)
       if(ierr /= 0) &
         call terminate("mixingIntervals", &
                        "Memory allocation failure for angles.")

       ! Set the angles to be checked.

       angles(1) = psi

       do nn=2,nAngles
         angles(nn) = psi + (nn-1)*perAngle
       enddo

       ! Make sure that all angles are in the range -pi .. pi.

       do nn=1,nAngles
         if(angles(nn) > pi) then
           angles(nn) = angles(nn) - two*pi
         else if(angles(nn) <= -pi) then
           angles(nn) = angles(nn) + two*pi
         endif
       enddo

       ! Determine the quadrants of the angles.

       do nn=1,nAngles
         if(angles(nn) > zero) then
           if(angles(nn) > half*pi) then
             quadrants(nn) = 2
           else
             quadrants(nn) = 1
           endif
         else
           if(angles(nn) > -half*pi) then
             quadrants(nn) = 4
           else
             quadrants(nn) = 3
           endif
         endif
       enddo

       ! Do an initial allocation for the local intervals.

       nLocal = 0
       nAlloc = 1000

       allocate(intervals(nAlloc), stat=ierr)
       if(ierr /= 0)                       &
         call terminate("mixingIntervals", &
                        "Memory allocation failure for intervals.")

       ! Loop over the local subfaces of the sliding interface and store
       ! the edges obtained by quadrilateral faces cut by the
       ! interpolation direction.

       domainLoop2: do nn=1,nDom

         ! Set the pointers for this block. Only the 1st spectral
         ! level needs to be considered.

         call setPointers(nn, level, 1_intType)

         ! Loop over the boundary conditions.

         bocoLoop2: do mm=1,nBocos

           ! Check if the subface belongs to the sliding mesh part
           ! currently investigated.

           testInterface2: if(BCType(mm)   == slidingInterface .and. &
                              groupNum(mm) == slideID) then

             ! Store the nodal range of the subface a bit easier.

             jBeg = BCData(mm)%jnBeg
             jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg
             iEnd = BCData(mm)%inEnd

             ! Allocate the memory for the nodal cylindrical coordinates
             ! of the subface and determine them.

             allocate(xFace(iBeg:iEnd,jBeg:jEnd,3), stat=ierr)
             if(ierr /= 0) &
               call terminate("mixingIntervals", &
                              "Memory allocation failure for xFace.")

             call cylCoorNodesOnBlockFace(xFace, iBeg, iEnd,  &
                                          jBeg,  jEnd, color, &
                                          BCFaceID(mm))

             ! If this is a radial interface, store the radii in the
             ! first position such that the same piece of code can be
             ! used for both radial and axial interfaces.

             if( radialInterface ) then
               do j=jBeg,jEnd
                 do i=iBeg,iEnd
                   xFace(i,j,1) = xFace(i,j,2)
                 enddo
               enddo
             endif

             ! Loop over the quadrilaterals of this subface and store
             ! the intersection edge of the ones that are cut by the
             ! interpolation direction.

             quadJLoop: do j=(jBeg+1),jEnd
               quadILoop: do i=(iBeg+1),iEnd

                 ! Store the 4 radii and polar angles of the quad a
                 ! bit easier.

                 rad(1) = xFace(i-1,j-1,1)
                 rad(2) = xFace(i,  j-1,1)
                 rad(3) = xFace(i,  j  ,1)
                 rad(4) = xFace(i-1,j,  1)

                 theta(1) = xFace(i-1,j-1,3)
                 theta(2) = xFace(i,  j-1,3)
                 theta(3) = xFace(i,  j  ,3)
                 theta(4) = xFace(i-1,j,  3)

                 ! Determine the difference between the minimum and
                 ! maximum radius.

                 rMin = min(rad(1), rad(2), rad(3), rad(4))
                 rMax = max(rad(1), rad(2), rad(3), rad(4))
                 dr   = rMax - rMin

                 ! Determine the quadrants of the nodes.

                 do jj=1,4
                   if(theta(jj) > zero) then
                     if(theta(jj) > half*pi) then
                       nk(jj) = 2
                     else
                       nk(jj) = 1
                     endif
                   else
                     if(theta(jj) > -half*pi) then
                       nk(jj) = 4
                     else
                       nk(jj) = 3
                     endif
                   endif
                 enddo

                 ! Loop over the number of angles to be checked.

                 angleLoop: do kk=1,nAngles

                   ! Initialize nCuts.

                   nCuts = 0

                   ! Loop over the edges of the quadrilateral.

                   edgeLoop: do jj=1,4

                     ! Store the other node of the edge in ii.

                     ii = jj + 1
                     if(ii == 5) ii = 1

                     ! Store the two polar angles of the nodes in
                     ! theta1 and theta2. They may have to be modified
                     ! to account for the singularity at theta == pi.

                     theta1 = theta(jj)
                     theta2 = theta(ii)

                     ! If one of the angles is in the second and the
                     ! other one is in the third quadrant, modify either
                     ! theta1 or theta2, depending on the interpolation
                     ! angle.

                     if(nk(jj) == 2 .and. nk(ii) == 3) then

                       ! Theta1 is in the second and theta2 in the third
                       ! quadrant. Check the interpolation direction and
                       ! either modify either theta1 or theta2.

                       if(quadrants(kk) == 3) then
                         theta1 = theta1 - two*pi
                       else
                         theta2 = theta2 + two*pi
                       endif

                     else if(nk(jj) == 3 .and. nk(ii) == 2) then

                       ! Theta1 is in the third and theta2 in the second
                       ! quadrant. Check the interpolation direction and
                       ! either modify either theta1 or theta2.

                       if(quadrants(kk) == 3) then
                         theta2 = theta2 - two*pi
                       else
                         theta1 = theta1 + two*pi
                       endif

                     endif

                     ! Check whether or not the edge intersects with
                     ! the current interpolation direction.

                     theta1 = theta1 - angles(kk)
                     theta2 = theta2 - angles(kk)

                     if(theta1*theta2 <= zero) then

                       ! Edge intersects. Store it.

                       nCuts = nCuts + 1
                       edgeCut(nCuts) = jj

                       ww1 = abs(theta2)/max(abs(theta2-theta1),eps)
                       ww2 = one - ww1
                       radCut(nCuts) = ww1*rad(jj) + ww2*rad(ii)

                     endif

                   enddo edgeLoop

                   ! Exit the angle loop if intersections were found.

                   if(nCuts > 0) exit angleLoop

                 enddo angleLoop

                 ! Determine the situation we are having here.

                 select case(nCuts)

                   case (4_intType)

                     ! Four intersections found. It is only a
                     ! theoretically possible case for a degenerated
                     ! triangle. It is always a regular edge.

                     if(nLocal == nAlloc) call reallocIntervals
                     nLocal = nLocal + 1

                     rMin = min(radCut(1), radCut(2), &
                                radCut(3), radCut(4))
                     rMax = max(radCut(1), radCut(2), &
                                radCut(3), radCut(4))

                     intervals(nLocal)%rMin         = rMin
                     intervals(nLocal)%rMax         = rMax
                     intervals(nLocal)%regularEdge = .true.

                   !=====================================================

                   case (3_intType)

                     ! 3 intersections. There is a theoretical
                     ! possibility that it the intersection line just
                     ! touches a degenerated quad. In that case the
                     ! minimum and maximum value of radCut is equal
                     ! and the edge should not be stored.

                     rMin = min(radCut(1), radCut(2), radCut(3))
                     rMax = max(radCut(1), radCut(2), radCut(3))

                     if((rMax - rMin) > 1.e-10_realType*dr) then

                       ! Edge must be stored. It is always a
                       ! regular edge.

                       if(nLocal == nAlloc) call reallocIntervals
                       nLocal = nLocal + 1

                       intervals(nLocal)%rMin         = rMin
                       intervals(nLocal)%rMax         = rMax
                       intervals(nLocal)%regularEdge = .true.

                     endif

                   !=====================================================

                   case (2_intType)

                     ! The normal case for an intersection. If the
                     ! integer difference between the cut edges is 2
                     ! then it is a regular face crossing; otherwise
                     ! it is irregular.

                     if(nLocal == nAlloc) call reallocIntervals
                     nLocal = nLocal + 1

                     intervals(nLocal)%rMin = min(radCut(1), radCut(2))
                     intervals(nLocal)%rMax = max(radCut(1), radCut(2))

                     if((edgeCut(2) - edgeCut(1)) == 2) then
                       intervals(nLocal)%regularEdge = .true.
                     else
                       intervals(nLocal)%regularEdge = .false.
                     endif

                 end select

               enddo quadILoop
             enddo quadJLoop

             ! Release the memory of the nodal coordinates again.

             deallocate(xFace, stat=ierr)
             if(ierr /= 0)                       &
               call terminate("mixingIntervals", &
                              "Deallocation failure for xFace.")

           endif testInterface2
         enddo bocoLoop2
       enddo domainLoop2

       ! Release the memory of angles; they are not needed anymore.

       deallocate(angles, stat=ierr)
       if(ierr /= 0)                       &
         call terminate("mixingIntervals", &
                        "Deallocation failure for angles.")
!
!      ******************************************************************
!      *                                                                *
!      * Determine the points of the interpolation intervals by         *
!      * gathering all the local edges and sorting them in increasing   *
!      * order.                                                         *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number of edges which will be received from
       ! every processor in allgatherv later on. Note the use of
       ! MPI_integer instead of sumb_integer for the receiving data
       ! type.

       call mpi_allgather(nLocal, 1, sumb_integer, recvcounts, 1, &
                          mpi_integer, comm, ierr)

       ! Create the array displs, needed in allgatherv later on.

       displs(1) = 0
       do nn=2,nProc
         displs(nn) = displs(nn-1) + recvcounts(nn-1)
       enddo

       ! Determine the global number of edges stored.

       nGlobal = displs(nProc) + recvcounts(nProc)

       ! Allocate the memory for the buffers used in allgatherv.

       allocate(intBuf1(nLocal),  realBuf1(2,nLocal), &
                intBuf2(nGlobal), realBuf2(2,nGlobal), stat=ierr)
       if(ierr /= 0)                       &
         call terminate("mixingIntervals", &
                        "Memory allocation failure for the buffers &
                        &used in allgatherv.")

       ! Copy the local information in intBuf1 and realBuf1.
       ! Release the memory of intervals afterwards.

       do nn=1,nLocal
         if( intervals(nn)%regularEdge ) then
           intBuf1(nn) = 1
         else
           intBuf1(nn) = 0
         endif

         realBuf1(1,nn) = intervals(nn)%rMin
         realBuf1(2,nn) = intervals(nn)%rMax
       enddo

       deallocate(intervals, stat=ierr)
       if(ierr /= 0)                       &
         call terminate("mixingIntervals", &
                        "Deallocation failure for intervals.")

       ! Use allgatherv to gather the data from all processors.
       ! First the integer data.

       nl = nLocal
       call mpi_allgatherv(intBuf1, nl, sumb_integer, intBuf2, &
                           recvcounts, displs, sumb_integer, comm, ierr)

       nl = 2*nLocal
       do nn=1,nProc
         recvcounts(nn) = 2*recvcounts(nn)
         displs(nn)     = 2*displs(nn)
       enddo

       call mpi_allgatherv(realBuf1, nl, sumb_real, realBuf2, &
                           recvcounts, displs, sumb_real, comm, ierr)

       ! Allocate the memory for the global edges and copy the
       ! information from intBuf2 and realBuf2.

       allocate(edges(nGlobal), stat=ierr)
       if(ierr /= 0)                       &
         call terminate("mixingIntervals", &
                        "Memory allocation failure for edges.")

       do nn=1,nGlobal
         if(intBuf2(nn) == 1) then
           edges(nn)%regularEdge = .true.
         else
           edges(nn)%regularEdge = .false.
         endif

         edges(nn)%rMin = realBuf2(1,nn)
         edges(nn)%rMax = realBuf2(2,nn)
       enddo

       ! Release the memory of the communication buffers again.

       deallocate(intBuf1, realBuf1, intBuf2, realBuf2, stat=ierr)
       if(ierr /= 0)                       &
         call terminate("mixingIntervals", &
                        "Deallocation failure for the buffers &
                        &used in allgatherv.")

       ! Sort edges in increasing order.

       call qsortMixingIntervalType(edges, nGlobal)

       ! Determine the true number of intervals by getting rid
       ! of the irregular edges.

       mm = 0
       nn = 0

       sortEdgesLoop: do

         ! Condition to exit the loop.

         if(nn == nGlobal) exit

         ! Update the counter nn and check the status of the edge.

         nn = nn + 1
         testRegular: if( edges(nn)%regularEdge ) then

           ! Regular edge. Update mm and store it.

           mm = mm + 1
           edges(mm) = edges(nn)

         else testRegular

           ! Irregular edge. Determine the situation we are having here.

           if(nn == nGlobal) then

             ! Last cut edge present. Add this edge to the previous
             ! edge computed.

             edges(mm-1)%rMax = edges(nn)%rMax

           else

             ! This cut edge and the next one are merged.
             ! This means that both mm and nn must be updated.
             ! The edge stored is the sum of nn-1 and nn.

             mm = mm + 1
             nn = nn + 1

             edges(mm)%rMin = edges(nn-1)%rMin
             edges(mm)%rMax = edges(nn)%rMax

           endif

         endif testRegular

       enddo sortEdgesLoop

       ! Store the true number of edges in nGlobal.

       nGlobal = mm

       ! Determine the number of points to define the
       ! interpolation intervals. This is the number of edges just
       ! determined + 2 (1 additional point plus halo at the maximum)
       ! + possibly halo at the minimum.

       nMixingPoints = nGlobal + 2

       ! Add a halo node to the minimum boundary, if desired.
       ! Mm will serve as a counter for the points.

       mm = 0
       if(nSingularQuad == 1 .or. addMinHalo) then
         nMixingPoints = nMixingPoints + 1
         mm = 1
       endif

       ! Allocate the memory for mixingPoints and mixingCells.

       allocate(mixingPoints(nMixingPoints), &
                mixingCells(nMixingPoints-1), stat=ierr)
       if(ierr /= 0) &
         call terminate("mixingIntervals", &
                        "Memory allocation failure for mixingPoints &
                        &and mixingCells.")

       ! Store the points.

       do nn=1,nGlobal
         mm = mm + 1
         mixingPoints(mm) = edges(nn)%rMin
       enddo

       ! Add the last point and the halo at the end of the interval.

       mm = mm + 1
       mixingPoints(mm) = edges(nGlobal)%rMax

       mm = nMixingPoints
       mixingPoints(mm) = two*mixingPoints(mm-1) - mixingPoints(mm-2)

       ! Add a halo node to the minimum boundary, if needed.

       if(nSingularQuad == 1) then

         ! Polar singular quad present. Start the interval at r == 0.

         mixingPoints(1) = zero

       else if( addMinHalo ) then

         ! A halo node must be added to the minimum boundary.
         ! Simply extrapolate.

         mixingPoints(1) = two*mixingPoints(2) - mixingPoints(3)

       endif

       ! Release the memory of edges.

       deallocate(edges, stat=ierr)
       if(ierr /= 0)                       &
         call terminate("mixingIntervals", &
                        "Deallocation failure for edges.")

       ! Determine the coordinates of mixingCells, i.e. the cell
       ! centered coordinates of the interpolation intervals.

       do nn=1,(nMixingPoints-1)
         mixingCells(nn) = half*(mixingPoints(nn) &
                          +       mixingPoints(nn+1))
       enddo

       !=================================================================

       contains

         !===============================================================

         subroutine reallocIntervals
!
!        ****************************************************************
!        *                                                              *
!        * ReallocIntervals reallocates the memory for the pointer      *
!        * intervals; the array is increased.                           *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Local variables.
!
         integer :: ierr

         integer(kind=intType) :: nn, nOld

         type(mixingIntervalType), dimension(:), pointer :: tmp
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         tmp => intervals
         nOld = nAlloc

         nAlloc = nAlloc + 1000
         allocate(intervals(nAlloc), stat=ierr)
         if(ierr /= 0)                        &
           call terminate("reallocIntervals", &
                          "Memory allocation failure for intervals.")

         do nn=1,nOld
           intervals(nn) = tmp(nn)
         enddo

         deallocate(tmp, stat=ierr)
         if(ierr /= 0)                        &
           call terminate("reallocIntervals", &
                          "Deallocation failure for tmp.")

         end subroutine reallocIntervals

       end subroutine mixingIntervals
