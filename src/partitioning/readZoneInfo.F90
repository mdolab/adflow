!
!      ******************************************************************
!      *                                                                *
!      * File:          readZoneInfo.F90                                *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 12-17-2002                                      *
!      * Last modified: 10-10-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readZoneInfo(cgnsBase, nZone, sortedFamName, &
                               famID, noUnits)
!
!      ******************************************************************
!      *                                                                *
!      * readZoneInfo reads the general information, like zone type     *
!      * and physical dimensions, for the given zone/block.             *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use communication
       use flowVarRefState
       use iteration
       use su_cgns
       use partitionMod
       implicit none
!
!      Subroutine arguments
!
       integer, intent(in) :: cgnsBase, nZone
       character(len=*), dimension(*), intent(in) :: sortedFamName
       integer(kind=intType), dimension(*), intent(in) :: famID

       logical, intent(inout) :: noUnits

#ifdef USE_NO_CGNS

       call terminate("readZoneInfo", &
                      "Routine should not be called if no cgns support &
                      &is selected.")
#else
!
!      Local variables
!
       integer :: cgnsInd
       integer :: i, ierr, nCoords
       integer :: mass, len, time, temp, angle

       integer, dimension(9) :: sizesBlock

       integer(kind=intType) :: ii, nn

       real, dimension(3) :: rotCenter, rotRate

       real(kind=realType) :: mult, trans

       character(len=maxCGNSNameLen) :: familyName
       character(len=maxStringLen)   :: errorMessage

       logical :: overwrite
!
!      Function definitions.
!
       integer(kind=intType) :: bsearchStrings
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set the cgns ID for the "master" file and read the size
       ! of the block as well as the zone name.

       cgnsInd = fileIDs(1)

       call cg_zone_read_f(cgnsInd, cgnsBase, nZone, &
                           cgnsDoms(nZone)%zoneName, sizesBlock, ierr)
       if(ierr /= all_ok)               &
         call terminate("readZoneInfo", &
                        "Something wrong when calling cg_nZones_f")

       ! Check the zone type.

       call cg_zone_type_f(cgnsInd, cgnsBase, nZone, &
                           cgnsDoms(nZone)%zonetype, ierr)
       if(ierr /= all_ok)               &
         call terminate("readZoneInfo", &
                        "Something wrong when calling cg_zone_type_f")

       if(cgnsDoms(nZone)%zonetype /= Structured) then
         write(errorMessage,*) "Zone ",                         &
                                trim(cgnsDoms(nZone)%zoneName), &
                                " of the grid file is not structured"
         if(myID == 0) call terminate("readZoneInfo", errorMessage)
         call mpi_barrier(SUmb_comm_world, ierr)
       endif

       ! Set the values for the number of nodes and cells in i, j and
       ! k-direction.

       cgnsDoms(nZone)%il = sizesBlock(1)
       cgnsDoms(nZone)%jl = sizesBlock(2)
       cgnsDoms(nZone)%kl = sizesBlock(3)

       cgnsDoms(nZone)%nx = sizesBlock(4)
       cgnsDoms(nZone)%ny = sizesBlock(5)
       cgnsDoms(nZone)%nz = sizesBlock(6)

       ! Check the size of this zone for the other grids to be read.
       ! They should be equal. Note that familyName is only used as
       ! a dummy in this call.

       do nn=2,nGridsRead
         call cg_zone_read_f(fileIDs(nn), cgnsBase, nZone, &
                             familyName, sizesBlock, ierr)
         if(ierr /= all_ok)               &
           call terminate("readZoneInfo", &
                          "Something wrong when calling cg_nZones_f")

         if(cgnsDoms(nZone)%il /= sizesBlock(1) .or. &
            cgnsDoms(nZone)%jl /= sizesBlock(2) .or. &
            cgnsDoms(nZone)%kl /= sizesBlock(3)) then

           write(errorMessage,*) "File ", trim(gridFiles(nn)), &
                                 ", zone ", trim(familyName), &
                                 " Zone dimensions are different&
                                 & than in file ", trim(gridFiles(1))
           if(myID == 0) call terminate("readBlockSizes", errorMessage)
           call mpi_barrier(SUmb_comm_world, ierr)
         endif
       enddo

       ! Goto this zone.

       call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", nZone, "end")
       if(ierr /= all_ok)               &
         call terminate("readZoneInfo", &
                        "Something wrong when calling cg_goto_f")
!
!      ******************************************************************
!      *                                                                *
!      * Try to read the family name.                                   *
!      *                                                                *
!      ******************************************************************
!
       call cg_famname_read_f(familyName, ierr)
       if(ierr == error)                &
         call terminate("readZoneInfo", &
                        "Something wrong when calling cg_famname_read_f")

       ! Check if a family name was specified. If so, determine the
       ! corresponding id.

       cgnsDoms(nZone)%familyID = 0
       if(ierr == all_ok) then

         ! Search the family name in the sorted names. For a valid
         ! grid this name must be found.

         nn = cgnsNFamilies
         ii = bsearchStrings(familyName, sortedFamName, nn)
         if(ii == 0) then

           write(errorMessage,100) trim(familyName)
 100       format("Family name",1X,A,1X,"not present in the grid")
           if(myID == 0) call terminate("readZoneInfo", errorMessage)
           call mpi_barrier(SUmb_comm_world, ierr)

         endif

         ! Set the family number.

         ii = famID(ii)
         cgnsDoms(nZone)%familyID = ii

       endif
!
!      ******************************************************************
!      *                                                                *
!      * Try to determine the units of the coordinates.                 *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number of coordinates in this zone.

       call cg_ncoords_f(cgnsInd, cgnsBase, nZone, nCoords, ierr)
       if(ierr /= all_ok)               &
         call terminate("readZoneInfo", &
                        "Something wrong when calling cg_ncoords_f")

       ! Check that 3 coordinates are present. If not, terminate.

       if(nCoords /= 3) then
         write(errorMessage,102) trim(cgnsDoms(nZone)%zoneName), nCoords
 102     format("The number of coordinates of zone ", a, &
                " of base 1 is", i1, ". This should 3.")

         if(myID == 0) call terminate("readZoneInfo", errorMessage)
         call mpi_barrier(SUmb_comm_world, ierr)
       endif

       ! Loop over the three coordinates.

       cgnsDoms(nZone)%gridUnitsSpecified = .false.

       do i=1,3

         ! Go to the correct place in the grid file.

         call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", nZone, &
                        "GridCoordinates_t", 1, "DataArray_t", i, &
                        "end")
         if(ierr /= all_ok)               &
           call terminate("readZoneInfo", &
                          "Something wrong when calling cg_goto_f")

         call cg_units_read_f(mass, len, time, temp, angle, ierr)
         if(ierr == error)                &
           call terminate("readZoneInfo", &
                          "Something wrong when calling cg_units_read_f")

         ! Check if units were specified.

         if(ierr == all_ok .and. len /= Null) then

           ! Copy the units and set gridUnitsSpecified to .true.

           cgnsDoms(nZone)%mass  = mass
           cgnsDoms(nZone)%len   = len
           cgnsDoms(nZone)%time  = time
           cgnsDoms(nZone)%temp  = temp
           cgnsDoms(nZone)%angle = angle

           cgnsDoms(nZone)%gridUnitsSpecified = .true.

           ! Determine the conversion factor to meters.

           call siLen(len, mult, trans)

           cgnsDoms(nZone)%LRef = mult

         endif

       enddo

       ! Check whether units were specified or not.

       if( cgnsDoms(nZone)%gridUnitsSpecified ) then

         ! Units were specified. Check if a global reference length
         ! was specified as well. If so, compare the conversion factor.
         ! If not identical, processor 0 prints a warning.

         if(myID ==0 .and. LRefSpecified .and. &
            cgnsDoms(nZone)%LRef /= LRef)  then

           print "(a)", "#"
           print "(a)", "#                      Warning"
           print 103, trim(cgnsDoms(nZone)%zoneName)
           print 104
           print "(a)", "#"
 103       format("# Zone",1X,A,": Conversion factor to meters &
                  &not identical to global conversion factor.")
 104       format("# Global conversion factor is ignored.")
         endif

         ! In case no global conversion factor was specified,
         ! set LRef to the LRef of this block.

         if(.not. LRefSpecified) LRef = cgnsDoms(nZone)%LRef

       else

         ! No units specified. Set the reference length of the block
         ! to the global conversion factor.

         cgnsDoms(nZone)%LRef = LRef

         ! If the global reference length was not specified, set
         ! noUnits to .true.

         if(.not. LRefSpecified) noUnits = .true.

       endif
!
!      ******************************************************************
!      *                                                                *
!      * Try to determine the rotation rate and center.                 *
!      *                                                                *
!      ******************************************************************
!
       ! Some initializations.

       cgnsDoms(nZone)%rotatingFrameSpecified = .false.
       cgnsDoms(nZone)%rotRate   = zero
       cgnsDoms(nZone)%rotCenter = zero

       mult = one   ! Assuming radians/s for the rotation rate.

       ! Check if a family overwrite is present.

       overwrite = .false.
       nn = cgnsDoms(nZone)%familyID
       if(nn > 0) then
         if( cgnsFamilies(nn)%rotatingFrameSpecified ) &
           overwrite = .true.
       endif

       testOverwrite: if( overwrite ) then

         ! Rotation rate is specified per family.
         ! Set rotatingFrameSpecified to .true. and copy the data
         ! from the corresponding family. Note that the rotation rate
         ! is already in rad/s and therefore mult does not need to
         ! change.

          cgnsDoms(nZone)%rotatingFrameSpecified = .true.
          cgnsDoms(nZone)%rotRate   = cgnsFamilies(nn)%rotRate
          cgnsDoms(nZone)%rotCenter = cgnsFamilies(nn)%rotCenter

       else testOverwrite

         ! Go to the correct location in the cgns file, where
         ! the rotation rate should be specified.

         call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", nZone, &
                        "end")
         if(ierr /= all_ok)               &
           call terminate("readZoneInfo", &
                          "Something wrong when calling cg_goto_f")

         ! No family information specified.
         ! Try to read the rotation rate and center.

         call cg_rotating_read_f(rotRate, rotCenter, ierr)
         if(ierr == error)                &
           call terminate("readZoneInfo", &
                          "Something wrong when calling &
                          &cg_rotating_read_f")

         ! Check if a rotating frame is specified.

         if(ierr == all_ok) then

           ! Set changingGrid to .true.

           changing_grid = .true.

           ! Set rotatingFrameSpecified to .true. and copy
           ! to rotation center and rotation rate.

           cgnsDoms(nZone)%rotatingFrameSpecified = .true.
           cgnsDoms(nZone)%rotRate   = rotRate
           cgnsDoms(nZone)%rotCenter = rotCenter

           ! Determine the conversion factor to rad/s if the
           ! dimensional units are specified. If not, it is assumed
           ! that the dimensions are given in rad/s.

           call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", nZone, &
                          "RotatingCoordinates_t", 1, "DataArray_t", 2, &
                          "end")
           if(ierr /= all_ok)               &
             call terminate("readZoneInfo", &
                            "Something wrong when calling cg_goto_f")

           call cg_units_read_f(mass, len, time, temp, angle, ierr)
           if(ierr == error)                &
             call terminate("readZoneInfo", &
                            "Something wrong when calling &
                            &cg_units_read_f")

           ! Check if units were specified. If not assume radians.

           if(ierr == all_ok .and. angle /= Null) then

             ! Determine the conversion factor to radIans.

             call siAngle(angle, mult, trans)

           else

             ! Angle units not specified. Assume rad/s.
             ! Processor 0 writes a warning to stdout.

             if(myID == 0) then

               print "(a)", "#"
               print "(a)", "#                      Warning"
               print 101, trim(cgnsDoms(nZone)%zoneName)
               print "(a)", "#"
 101           format("# Zone",1X,A,": No unit specified for &
                      &the rotation rate, assuming rad/s.")
             endif

           endif
         endif

       endif testOverwrite

       ! Multiply the rotation center with LRef to obtain the correct
       ! coordinates in meters.

       cgnsDoms(nZone)%rotCenter(1) = cgnsDoms(nZone)%LRef &
                                    * cgnsDoms(nZone)%rotCenter(1)
       cgnsDoms(nZone)%rotCenter(2) = cgnsDoms(nZone)%LRef &
                                    * cgnsDoms(nZone)%rotCenter(2)
       cgnsDoms(nZone)%rotCenter(3) = cgnsDoms(nZone)%LRef &
                                    * cgnsDoms(nZone)%rotCenter(3)

       ! Multiply the rotation rate by the conversion factor to obtain
       ! the correct angular velocity.

       cgnsDoms(nZone)%rotRate(1) = cgnsDoms(nZone)%rotRate(1)*mult
       cgnsDoms(nZone)%rotRate(2) = cgnsDoms(nZone)%rotRate(2)*mult
       cgnsDoms(nZone)%rotRate(3) = cgnsDoms(nZone)%rotRate(3)*mult

#endif

       end subroutine readZoneInfo
