!
!      ******************************************************************
!      *                                                                *
!      * File:          determineSections.f90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 12-09-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine determineSections
!
!      ******************************************************************
!      *                                                                *
!      * determineSections determines the number of sections, i.e.      *
!      * grid parts between sliding mesh interfaces, present in the     *
!      * entire grid.                                                   *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use block
       use cgnsGrid
       use communication
       use inputTimeSpectral
       use section
       use su_cgns
       implicit none
!
!      Local parameter, threshold for allowed angle difference between
!                       the theoretical and true value, 0.1 degrees.
!
       real(kind=realType), parameter :: threshold = 0.1_realType
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn, mm, ii, jj
       integer(kind=intType) :: nLevel, nSlices, slideID
       integer(kind=intType), dimension(cgnsNDom) :: sectionID, sorted
       integer(kind=intType), dimension(cgnsNsliding,2) :: secSliding

       real(kind=realType) :: cosTheta, cosPhi, cosPsi
       real(kind=realType) :: sinTheta, sinPhi, sinPsi
       real(kind=realType) :: r11, r12, r13, r21, r22, r23
       real(kind=realType) :: r31, r32, r33
       real(kind=realType) :: d1, d2, a1, a0, lamr, lami, angle, dAngle

       logical :: situationChanged
!
!      Function definition.
!
       integer(kind=intType) :: bsearchIntegers
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize sectionID to the cgns block number.

       do nn=1,cgnsNDom
         sectionID(nn) = nn
       enddo

       ! Loop to determine the highest cgns block id in every section.

       loopHighestBlock: do

         ! Initialize situationChanged to .false.

         situationChanged = .false.

         ! Loop over the internal block faces of the blocks.

         do nn=1,cgnsNDom

           ! First the 1 to 1 block connectivities.

           do mm=1,cgnsDoms(nn)%n1to1

             ! Store the neighboring block a bit easier and change the
             ! value of sectionID if the neighboring block has a
             ! higher value. In that case situationChanged is .true.

             ii = cgnsDoms(nn)%conn1to1(mm)%donorBlock
             if(sectionID(ii) > sectionID(nn)) then
               sectionID(nn)    = sectionID(ii)
               situationChanged = .true.
             endif

           enddo

           ! No general connectivities yet.

         enddo

         ! Criterion to exit the loop.

         if(.not. situationChanged) exit

       enddo loopHighestBlock

       ! Copy sectionID in sorted and sort sorted in increasing order.

       do nn=1,cgnsNDom
         sorted(nn) = sectionID(nn)
       enddo

       call qsortIntegers(sorted, cgnsNDom)

       ! Determine the number of sections. Note there is at least one
       ! section, because there is at least one block in the grid.

       nSections = 1
       do nn=2,cgnsNDom
         if(sorted(nn) > sorted(nSections)) then
           nSections = nSections +1
           sorted(nSections) = sorted(nn)
         endif
       enddo

       ! Determine the sections to which the owned blocks belong.

       nLevel = ubound(flowDoms,2)
       do nn=1,nDom

         ! Determine the corresponding cgns block id and search its
         ! section id in sorted.

         mm = flowDoms(nn,1,1)%cgnsBlockID
         mm = bsearchIntegers(sectionID(mm), sorted, nSections)

         if( debug ) then
           if(mm == 0) call terminate("determineSections", &
                                      "Entry not found in sorted.")
         endif

         ! Set the section id for all grid levels for all spectral
         ! time intervals.

         do jj=1,nTimeIntervalsSpectral
           do ii=1,nLevel
             flowDoms(nn,ii,jj)%sectionID = mm
           enddo
         enddo

       enddo

       ! Allocate the memory for sections.

       allocate(sections(nSections), stat=ierr)
       if(ierr /= 0)                         &
         call terminate("determineSections", &
                        "Memory allocation failure for sections.")

       ! Initialize the number of slices for each of the sections to 1,
       ! periodic and rotating to .false. and rotCenter, rotAxis
       ! and rotRate to zero.

       do nn=1,nSections
         sections(nn)%nSlices   = 1
         sections(nn)%periodic  = .false.
         sections(nn)%rotating  = .false.
         sections(nn)%rotCenter = zero
         sections(nn)%rotAxis   = zero
         sections(nn)%rotRate   = zero
       enddo

       ! Determine the number of slices and the periodic transformation
       ! for the sections.

       loopCGNSDom: do nn=1,cgnsNDom

         ! Search for the corresponding section.

         ii = bsearchIntegers(sectionID(nn), sorted, nSections)
         if( debug ) then
           if(ii == 0) call terminate("determineSections", &
                                      "Entry not found in sorted.")
         endif

         ! It is assumed that periodic info is correct. So if this
         ! section has already been treated, there is no need to do
         ! it again.

         if( sections(ii)%periodic ) cycle

         ! Loop over the 1 to 1 subfaces of the cgns block and try to
         ! find a periodic one.

         do mm=1,cgnsDoms(nn)%n1to1
           if( cgnsDoms(nn)%conn1to1(mm)%periodic ) exit
         enddo

         ! Continue with the next block if this block does not have
         ! periodic subfaces.

         if(mm > cgnsDoms(nn)%n1to1) cycle

         ! Subface mm is a periodic one. Set periodic to .true.

         sections(ii)%periodic = .true.

         ! Set the rotation axis of the section to the rotation
         ! angles of the periodic transformation. This may be
         ! overwritten later on using the rotation rate, but for
         ! some cases this is the only rotation information present.

         sections(ii)%rotAxis = cgnsDoms(nn)%conn1to1(mm)%rotationAngles

         ! Construct the rotation matrix, where it is assumed that the
         ! sequence of rotation is first rotation around the x-axis,
         ! followed by rotation around the y-axis and finally rotation
         ! around the z-axis.

         cosTheta = cos(cgnsDoms(nn)%conn1to1(mm)%rotationAngles(1))
         sinTheta = sin(cgnsDoms(nn)%conn1to1(mm)%rotationAngles(1))

         cosPhi = cos(cgnsDoms(nn)%conn1to1(mm)%rotationAngles(2))
         sinPhi = sin(cgnsDoms(nn)%conn1to1(mm)%rotationAngles(2))

         cosPsi = cos(cgnsDoms(nn)%conn1to1(mm)%rotationAngles(3))
         sinPsi = sin(cgnsDoms(nn)%conn1to1(mm)%rotationAngles(3))

         r11 =  cosPhi*cosPsi
         r21 =  cosPhi*sinPsi
         r31 = -sinPhi

         r12 = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi
         r22 = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi
         r32 = sinTheta*cosPhi

         r13 = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi
         r23 = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi
         r33 = cosTheta*cosPhi

         ! Store the rotation matrix, rotation center and translation
         ! vector for this section.

         sections(ii)%rotCenter  = &
                      cgnsDoms(nn)%conn1to1(mm)%rotationCenter
         sections(ii)%translation = &
                      cgnsDoms(nn)%conn1to1(mm)%translation

         sections(ii)%rotMatrix(1,1) = r11
         sections(ii)%rotMatrix(2,1) = r21
         sections(ii)%rotMatrix(3,1) = r31

         sections(ii)%rotMatrix(1,2) = r12
         sections(ii)%rotMatrix(2,2) = r22
         sections(ii)%rotMatrix(3,2) = r32

         sections(ii)%rotMatrix(1,3) = r13
         sections(ii)%rotMatrix(2,3) = r23
         sections(ii)%rotMatrix(3,3) = r33

         ! Determine the coefficients of lambda and lambda^2 of the
         ! characteristic polynomial.

         d2 = -r11 - r22 - r33
         d1 =  r11*r22 + r11*r33 + r22*r33 - r12*r21 - r13*r31 - r23*r32

         ! Make use of the fact that one eigenvalue of the transformation
         ! matrix is 1 and determine the coefficients of the quadratic
         ! equation for the other two eigenvalues.

         a1 = d2 + one
         a0 = a1 + d1

         ! Determine the real and imaginary part of the two eigenvalues.
         ! Neglect the factor 1/2 here.

         lamr = -a1
         lami = sqrt(abs(a1*a1 - four*a0))

         ! Determine the angle in the imaginary plane. Due to the
         ! positive definition of lami, this angle will be between
         ! 0 and pi. Take care of the exceptional case that the angle
         ! is zero.

         angle = atan2(lami, lamr)
         if(angle == zero) angle = two*pi

         ! Determine the number of slices.

         nSlices = nint(two*pi/angle)

         ! Determine the angle difference in degrees between
         ! nSlices*angle and a complete rotation. If this is larger than
         ! the threshold processor 0 will print an error message and exit.

         dAngle = abs(180.0_realType*(two*pi - nSlices*angle)/pi)

         if(dAngle >= threshold) then
           if(myID == 0)                         &
             call terminate("determineSections", &
                            "Periodic angle not a integer divide of &
                            &360 degrees")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         ! Set the number of slices for this section.

         sections(ii)%nSlices = nSlices

       enddo loopCGNSDom

       ! Again loop over the number of block of the original mesh,
       ! but now determine whether or not the section is rotating.

       do nn=1,cgnsNDom

         ! If the block is rotating, copy that information to
         ! the corresponding section. If the section is not
         ! periodic, also set the rotation center.

         if( cgnsDoms(nn)%rotatingFrameSpecified ) then

           ii = bsearchIntegers(sectionID(nn), sorted, nSections)
           sections(ii)%rotating = .true.
           sections(ii)%rotAxis  = cgnsDoms(nn)%rotRate
           sections(ii)%rotRate  = cgnsDoms(nn)%rotRate

           if(.not. sections(ii)%periodic) &
             sections(ii)%rotCenter = cgnsDoms(nn)%rotCenter

         endif
       enddo

       ! Determine the two sections for every sliding mesh
       ! interface.

       secSliding = 0
       do nn=1,cgnsNDom
         do mm=1,cgnsDoms(nn)%nBocos
           if(cgnsDoms(nn)%bocoInfo(mm)%actualFace .and. &
              cgnsDoms(nn)%bocoInfo(mm)%BCType == SlidingInterface) then

             ! Boundary face is part of a sliding mesh interface.
             ! Determine the ID of the interface.

             slideID = abs(cgnsDoms(nn)%bocoInfo(mm)%slidingID)

             ! Determine the section to which this block belongs and
             ! store its id for this sliding interface.

             ii = bsearchIntegers(sectionID(nn), sorted, nSections)
             if(secSliding(slideID,1) == 0) then
               secSliding(slideID,1) = ii
             else if(secSliding(slideID,1) /= ii) then
               secSliding(slideID,2) = ii
             endif

           endif
         enddo
       enddo

       ! Loop over the sliding mesh interfaces to set the rotation axis
       ! for non-rotating sections.

       do ii=1,cgnsNsliding

         ! Store the two sections id's a bit easier.

         mm = secSliding(ii,1)
         nn = secSliding(ii,2)

         ! Print an error message if both sections are not rotating.

         if((.not. sections(mm)%rotating) .and. &
            (.not. sections(nn)%rotating) ) then
           if(myID == 0)                         &
             call terminate("determineSections", &
                            "Encountered sliding interface between &
                            &two non-rotating sections")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         ! Set the rotation axis if section mm is not rotating.
         ! If it is not periodic also set the rotation point.

         if(.not. sections(mm)%rotating) then
           sections(mm)%rotAxis = sections(nn)%rotAxis
           if(.not. sections(mm)%periodic) &
             sections(mm)%rotCenter = sections(nn)%rotCenter
         endif

         ! Idem for section nn.

         if(.not. sections(nn)%rotating) then
           sections(nn)%rotAxis = sections(mm)%rotAxis
           if(.not. sections(nn)%periodic) &
             sections(nn)%rotCenter = sections(mm)%rotCenter
         endif

       enddo

       ! Determine the unit rotation axis for the sections.

       do nn=1,nSections
         d1 = one/max(eps,sqrt(sections(nn)%rotAxis(1)**2 &
            +                  sections(nn)%rotAxis(2)**2 &
            +                  sections(nn)%rotAxis(3)**2))

         sections(nn)%rotAxis(1) = d1*sections(nn)%rotAxis(1)
         sections(nn)%rotAxis(2) = d1*sections(nn)%rotAxis(2)
         sections(nn)%rotAxis(3) = d1*sections(nn)%rotAxis(3)
       enddo

       end subroutine determineSections
