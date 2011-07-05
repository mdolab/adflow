!
!      ******************************************************************
!      *                                                                *
!      * File:          readGrid.F90                                    *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 12-31-2002                                      *
!      * Last modified: 02-26-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readGrid
!
!      ******************************************************************
!      *                                                                *
!      * readGrid reads the coordinates for the blocks or block parts   *
!      * to be stored on this processor.                                *
!      *                                                                *
!      ******************************************************************
!
#ifdef USE_NO_CGNS
       call terminate("readGrid", &
                      "Routine should not be called if no CGNS support &
                      &is selected.")
#else
       use block
       use cgnsGrid
       use cgnsNames
       use communication
       use flowVarRefState
       use inputIO
       use IOModule
       use su_cgns
       use partitionMod
       implicit none
!
!      Local variables.
!
       integer :: cgnsInd, cgnsBase, cgnsZone
       integer :: j, nCoords
       integer :: ierr, realTypeCGNS, datatype

       integer, dimension(3) :: rangeMin, rangeMax

       integer(kind=intType) :: i, ii, jj, kk, ll, nn
       integer(kind=intType) :: il, jl, kl
       integer(kind=intType) :: typeMismatch

       character(len=7)              :: int1String
       character(len=2*maxStringLen) :: errorMessage
       character(len=maxCGNSNameLen) :: coordname

       real(kind=cgnsRealType), allocatable, dimension(:,:,:) :: buffer
!
!      Function definition.
!
       integer :: setCGNSRealType
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set the cgns real type and initialize typeMismatch to 0.
       ! Set cgnsBase to 1, because we will always read from base 1;
       ! possible higher bases are ignored.

       realTypeCGNS = setCGNSRealType()
       typeMismatch = 0
       cgnsBase     = 1

       ! Loop over the number of blocks stored on this processor.

       domainLoop: do i=1,nDom

         ! Abbreviate the nodal block dimensions.

         il = flowDoms(i,1,1)%il
         jl = flowDoms(i,1,1)%jl
         kl = flowDoms(i,1,1)%kl

         ! Store the zone number a bit easier and set the range
         ! for reading the coordinates.

         cgnsZone = flowDoms(i,1,1)%cgnsBlockID

         rangeMin(1) = flowDoms(i,1,1)%iBegor
         rangeMin(2) = flowDoms(i,1,1)%jBegor
         rangeMin(3) = flowDoms(i,1,1)%kBegor

         rangeMax(1) = flowDoms(i,1,1)%iEndor
         rangeMax(2) = flowDoms(i,1,1)%jEndor
         rangeMax(3) = flowDoms(i,1,1)%kEndor

         ! Allocate the memory for the read buffer.

         allocate(buffer(il,jl,kl), stat=ierr)
         if(ierr /= 0)                &
           call terminate("readGrid", &
                          "Memory allocation error for buffer")

         ! Loop over the number of grids to be read.

         nGridLoop: do nn=1,nGridsRead

           ! Store the file index a bit easier.

           cgnsInd = fileIDs(nn)

           ! Determine the number of coordinates in this zone.

           call cg_ncoords_f(cgnsInd, cgnsBase, cgnsZone, &
                             nCoords, ierr)
           if(ierr /= all_ok)           &
             call terminate("readGrid", &
                            "Something wrong when calling cg_ncoords_f")

           ! The coordinates are only read if 3 coordinates are present.

           checkNcoords: if(nCoords == 3) then

             ! Loop over the number of coordinates. Note that the counter j
             ! is an integer. This is for compatibility with cgns.

             coords: do j=1,nCoords

               ! Get the info for this coordinate.

               call cg_coord_info_f(cgnsInd, cgnsBase, cgnsZone, j, &
                                    datatype, coordname, ierr)
               if(ierr /= all_ok)           &
                 call terminate("readGrid", &
                                "Something wrong when calling &
                                &cg_coord_info_f")

               ! Update the value of typeMismatch if the datatype of
               ! the coordinate is not equal to the datatype used in
               ! the solver.

               if(realTypeCGNS /= datatype) &
                 typeMismatch = typeMismatch + 1

               ! Set the value of the counter ll, depending on the name.
               ! Normally the x-coordinate is written first, followed by
               ! the y-coordinate and finally the z-coordinate. But you
               ! never know.

               select case(coordname)
                 case (cgnsCoorX)
                   ll = 1
                 case (cgnsCoorY)
                   ll = 2
                 case (cgnsCoorZ)
                   ll = 3
                 case default
                   write(errorMessage,110)                        &
                               trim(cgnsDoms(cgnsZone)%zoneName), &
                               trim(coordname)
 110               format("Zone ",a," :Unknown coordinate name, ",a, &
                          ",in grid file")
                   call terminate("readGrid", errorMessage)
               end select

               ! Read the coordinates.

               call cg_coord_read_f(cgnsInd, cgnsBase, cgnsZone, &
                                    coordname, realTypeCGNS,     &
                                    rangeMin, rangeMax, buffer, ierr)
               if(ierr /= all_ok)           &
                 call terminate("readGrid", &
                                "Something wrong when calling &
                                &cg_coord_read_f")

               ! Copy the data into IOVar and scale it to meters.

               do kk=1,kl
                 do jj=1,jl
                   do ii=1,il
                     IOVar(i,nn)%w(ii,jj,kk,ll) = buffer(ii,jj,kk) &
                                                * cgnsDoms(cgnsZone)%LRef
                   enddo
                 enddo
               enddo

             enddo coords

           else checkNcoords

             ! There are not three coordinates present in this base.
             ! An error message is printed and an exit is made.

             write(errorMessage,101) trim(gridFiles(nn)),               &
                                     trim(cgnsDoms(cgnsZone)%zoneName), &
                                     nCoords
 101         format("File ", a, ": The number of coordinates of zone ", &
                    a, " should be 3, not ", i1)

             call terminate("readGrid", errorMessage)

           endif checkNcoords

         enddo nGridLoop

         ! Release the memory of buffer.

         deallocate(buffer, stat=ierr)
         if(ierr /= 0) call terminate("readGrid", &
                                      "Deallocation error for buffer")
       enddo domainLoop

       ! Close the cgns files.

       do nn=1,nGridsRead
         call cg_close_f(fileIDs(nn), ierr)
         if(ierr /= all_ok)           &
           call terminate("readGrid", &
                          "Something wrong when calling cg_close_f")
       enddo

       ! Determine the global sum of typeMismatch; the result only
       ! needs to be known on processor 0. Use ii as the global buffer
       ! to store the result. If a type mismatch occured,
       ! print a warning.

       call mpi_reduce(typeMismatch, ii, 1, sumb_integer, &
                       mpi_sum, 0, SUmb_comm_world, ierr)
       if(myID == 0 .and. ii > 0) then

         write(int1String,"(i6)") ii
         int1String = adjustl(int1String)

         print "(a)", "#"
         print "(a)", "#                      Warning"
         print 120, trim(int1String)
         print "(a)", "#"
 120     format("# ",a," type mismatches occured when reading the &
                        &coordinates of the blocks")
       endif

       ! If the coordinates in the solution files must be written in
       ! meters, correct this info for all cgns blocks.

       if( writeCoorMeter ) then
         do i=1,cgnsNDom
           cgnsDoms(i)%mass  = Null
           cgnsDoms(i)%len   = Meter
           cgnsDoms(i)%time  = Null
           cgnsDoms(i)%temp  = Null
           cgnsDoms(i)%angle = Null
           
           cgnsDoms(i)%gridUnitsSpecified = .true.
           cgnsDoms(i)%LRef = one
         enddo

         LRef = one
       endif

#endif

       end subroutine readGrid
