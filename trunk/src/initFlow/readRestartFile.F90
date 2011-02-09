!
!      ******************************************************************
!      *                                                                *
!      * File:          readRestartFile.F90                             *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 03-07-2003                                      *
!      * Last modified: 09-13-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readRestartFile(halosRead)
!
!      ******************************************************************
!      *                                                                *
!      * readRestartFile reads the fine grid solution(s) from the       *
!      * restart file(s). If the restart file(s) do not correspond to   *
!      * the current mesh, the solution(s) are interpolated onto this   *
!      * mesh. It is also allowed to change boundary conditions, e.g.   *
!      * an alpha and/or Mach sweep is possible. Furthermore there is   *
!      * some support when starting from a different turbulence model,  *
!      * although this should be used with care.                        *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use cgnsGrid
       use communication
       use flowVarRefState
       use inputIO
       use inputPhysics
       use inputTimeSpectral
       use iteration
       use monitor
       use su_cgns
       use restartMod
       implicit none
!
!      Subroutine arguments
!
       logical, intent(out) :: halosRead
!
!      Local variables.
!
       integer :: nZones, cellDim, physDim, ierr, nSols

       integer, dimension(9) :: sizes
       integer, dimension(nSolsRead) :: fileIDs

       integer(kind=intType) :: ii, jj, nn
       integer(kind=intType) :: nTypeMismatch
       integer(kind=intType) :: nHiMin, nHjMin, nHkMin
       integer(kind=intType) :: nHiMax, nHjMax, nHkMax

       character(len=7)              :: integerString
       character(len=maxCGNSNameLen) :: cgnsName
       character(len=2*maxStringLen) :: errorMessage
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
#ifdef USE_NO_CGNS

       call terminate("readRestartFile", &
                      "Routine should not be called if no cgns support &
                      &is selected.")

#else
       ! Initialize halosRead to .true. This will be overwritten if
       ! there is at least one block present for which the halo data
       ! cannot be read.

       halosRead = .true.

       ! Initialize nTypeMismatch to 0.

       nTypeMismatch = 0

       ! Loop over the number of files to be read and open them.

       fileOpenLoop: do solID=1,nSolsRead

         ! Open the restart file for reading.

         call cg_open_f(solFiles(solID), mode_read, cgnsInd, ierr)
         if(ierr /= all_ok) then
           write(errorMessage,*) "File ", trim(solFiles(solID)), &
                                 " could not be opened for reading"
           call terminate("readRestartFile", errorMessage)
         endif

         fileIDs(solID) = cgnsInd

         ! Determine the number of bases in the cgns file.
         ! This must be at least 1.

         call cg_nbases_f(cgnsInd, cgnsBase, ierr)
         if(ierr /= all_ok)                  &
           call terminate("readRestartFile", &
                          "Something wrong when calling cg_nbases_f")

         if(CGNSBase < 1) then
           write(errorMessage,*) "CGNS file ", trim(solFiles(solID)), &
                                 " does not contain a base"
           call terminate("readRestartFile", errorMessage)
         endif

         ! Only data from the first base is read. Information from
         ! higher bases is ignored.

         cgnsBase = 1

         ! Read the cell and physical dimensions as well as the name for
         ! this base.

         call cg_base_read_f(cgnsInd, cgnsBase, cgnsName, cellDim, &
                             physDim, ierr)
         if(ierr /= all_ok)                  &
           call terminate("readRestartFile", &
                          "Something wrong when calling cg_base_read_f")

         ! Check the cell and physical dimensions. Both must be 3 for
         ! this code to work.

         if(cellDim /= 3 .or. physDim /= 3) then
           write(errorMessage,100) cellDim, physDim
 100       format("Both the number of cell and physical dimensions &
                  &should be 3, not",1X,I1,1X,"and",1X,I1)
           call terminate("readRestartFile", errorMessage)
         endif

       enddo fileOpenLoop

       ! Read the convergence history. Only processor 0 needs to do that.

       if(myID == 0) call readConvHistory(fileIDs)

       ! Read the time history for an unsteady computation. Again only
       ! done by processor 0. Note that time history only needs to be
       ! present in the first solution file

       if(equationMode == unsteady .and. myID == 0) &
         call readTimeHistory(fileIDs)

       ! Broadcast nTimeStepsRestart and timeUnsteadyRestart to all
       ! processors. These values are needed to perform a consistent
       ! unsteady restart.

       call mpi_bcast(nTimeStepsRestart, 1, sumb_integer, 0, &
                      SUmb_comm_world, ierr)
       call mpi_bcast(timeUnsteadyRestart, 1, sumb_real, 0, &
                      SUmb_comm_world, ierr)

       ! Get the scaling factors for density, pressure and velocity
       ! by reading the reference state.

       call scaleFactors(fileIDs)

       ! Loop over the number of files to be read and read the solution.

       solLoop: do solID=1,nSolsRead

         ! Store the file index a bit easier and set the base to 1.

         cgnsInd  = fileIDs(solID)
         cgnsBase = 1

         ! Determine the number of zones (blocks) in the restart file
         ! and check if this is identical to the number in the grid file.

         call cg_nzones_f(cgnsInd, cgnsBase, nZones, ierr)
         if(ierr /= all_ok)                  &
           call terminate("readRestartFile", &
                          "Something wrong when calling cg_nzones_f")

         if(nZones /= cgnsNdom)              &
           call terminate("readRestartFile", &
                          "Number of blocks in grid file and restart &
                          &file differ")

         ! Create a sorted version of the zone names of the restart file
         ! and store its corresponding zone numbers in zoneNumbers.

         call getSortedZoneNumbers

         ! Loop over the number of blocks stored on this processor.

         domains: do nn=1,nDom

           ! Set the pointers for this block. Make sure that the
           ! correct data is set.

           ii = min(solID,nTimeIntervalsSpectral)
           call setPointers(nn, 1_intType, ii)

           ! Store the zone name of the original grid a bit easier.

           cgnsName = cgnsDoms(nbkGlobal)%zoneName

           ! Search in the sorted zone names of the restart file for
           ! cgnsName. The name must be found; otherwise the restart
           ! is pointless. If found, the zone number is set accordingly.

           jj = bsearchStrings(cgnsname, zoneNames, cgnsNdom)
           if(jj == 0) then
             write(errorMessage,*) "Zone name ", trim(cgnsName),  &
                                   " not found in restart file ", &
                                   trim(solFiles(solID))
             call terminate("readRestartFile", errorMessage)
           else
             jj = zoneNumbers(jj)
           endif

           cgnsZone = jj

           ! Determine the dimensions of the zone and check if these are
           ! identical to the dimensions of the block in the grid file.

           call cg_zone_read_f(cgnsInd, cgnsBase, cgnsZone, &
                               cgnsname, sizes, ierr)
           if(ierr /= all_ok)                  &
             call terminate("readRestartFile", &
                            "Something wrong when calling &
                            &cg_zone_read_f")

           if(cgnsDoms(nbkGlobal)%il /= sizes(1) .or. &
              cgnsDoms(nbkGlobal)%jl /= sizes(2) .or. &
              cgnsDoms(nbkGlobal)%kl /= sizes(3))     &
             call terminate("readRestartFile", &
                            "Corresponding zones in restart file and &
                            &grid file have different dimensions")

           ! Determine the number of flow solutions in this zone and
           ! check if there is a solution stored.

           call cg_nsols_f(cgnsInd, cgnsBase, cgnsZone, nSols, ierr)
           if(ierr /= all_ok)                  &
             call terminate("readRestartFile", &
                            "Something wrong when calling cg_nsols_f")

           if(nSols == 0)                      &
             call terminate("readRestartFile", &
                            "No solution present in restart file")

           ! Check for multiple solutions. A distinction is needed for
           ! overset cases because there will be an extra solution node
           ! for the nodal iblanks.

           if((nSols > 1 .and. .not. oversetPresent) .or. nSols > 2) &
             call terminate("readRestartFile", &
                            "Multiple solutions present in restart file")

           ! Determine the location of the solution variables. A loop is
           ! done over the solution nodes which is either 1 or 2. In the
           ! latter case, pick the node not named "Nodal Blanks".

           do cgnsSol=1,nSols
             call cg_sol_info_f(cgnsInd, cgnsBase, cgnsZone, cgnsSol, &
                                cgnsName, location, ierr)
             if(ierr /= all_ok)                  &
               call terminate("readRestartFile", &
                              "Something wrong when calling &
                              &cg_sol_info_f")

             if (trim(cgnsName) /= "Nodal Blanks") exit
           end do

           ! Determine the rind info.

           call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", &
                          cgnsZone, "FlowSolution_t", cgnsSol, "end")
           if(ierr /= all_ok)                  &
             call terminate("readRestartFile", &
                            "Something wrong when calling cg_goto_f")

           call cg_rind_read_f(sizes, ierr)
           if(ierr /= all_ok)                  &
             call terminate("readRestartFile", &
                            "Something wrong when calling &
                            &cg_rind_read_f")

           ! Check if halo's are present. If not, set halosRead to .false.
           ! This only needs to be done if this is not an older state for
           ! an unsteady computation.

           if(solID == 1 .or. equationMode == timeSpectral) then
             if(sizes(1) == 0 .or. sizes(2) == 0 .or. sizes(3) == 0 .or. &
                sizes(4) == 0 .or. sizes(5) == 0 .or. sizes(6) == 0)     &
               halosRead = .false.
           endif

           ! Initialize the number of halo cells to read to 0.

           nHiMin = 0; nHjMin = 0; nHkMin = 0
           nHiMax = 0; nHjMax = 0; nHkMax = 0

           ! Determine the range which must be read. A few things must be
           ! taken into account: - in iBegor, iEndor, etc. The nodal
           !                       range is stored. As in CGNS the cell
           !                       range start at 1, 1 must be subtracted
           !                       from the upper bound.
           !                     - the rind info must be taken into
           !                       account, because only the upper bound
           !                       is changed in cgns; the lower bound
           !                       remains 1.
           !                     - in case the solution is stored in the
           !                       vertices one extra variable in each
           !                       direction is read. An averaging will
           !                       take place to obtain cell centered
           !                       values.
           ! Also when vertex data is present, set halosRead to .false.,
           ! because it is not possible to determine the halos.

           if(location == CellCenter) then

             ! Correct the number of halo cells to be read.
             ! Only if this is not an older state in time for an
             ! unsteady computation.

             if(solID == 1 .or. equationMode == timeSpectral) then
               if(sizes(1) > 0) nHiMin = 1; if(sizes(2) > 0) nHiMax = 1
               if(sizes(3) > 0) nHjMin = 1; if(sizes(4) > 0) nHjMax = 1
               if(sizes(5) > 0) nHkMin = 1; if(sizes(6) > 0) nHkMax = 1
             endif

             ! Set the cell range to be read from the CGNS file.

             rangeMin(1) = iBegOr + sizes(1) - nHiMin
             rangeMin(2) = jBegOr + sizes(3) - nHjMin
             rangeMin(3) = kBegOr + sizes(5) - nHkMin

             rangeMax(1) = rangeMin(1) + nx-1 + nHiMin + nHiMax
             rangeMax(2) = rangeMin(2) + ny-1 + nHjMin + nHjMax
             rangeMax(3) = rangeMin(3) + nz-1 + nHkMin + nHkMax

           else if(location == Vertex) then

             ! Set the nodal range such that enough info is present
             ! to average the nodal data to create the cell centered
             ! data in the owned cells. No halo cells will be
             ! initialized.

             halosRead   = .false.

             rangeMin(1) = iBegor + sizes(1)
             rangeMin(2) = jBegor + sizes(3)
             rangeMin(3) = kBegor + sizes(5)

             rangeMax(1) = rangeMin(1) + nx
             rangeMax(2) = rangeMin(2) + ny
             rangeMax(3) = rangeMin(3) + nz
           else
             call terminate("readRestartFile", &
                            "Only CellCenter or Vertex data allowed in &
                            &restart file")
           endif

           ! Allocate the memory for buffer, needed to store the variable
           ! to be read, and bufferVertex in case the solution is stored
           ! in the vertices.

           allocate(buffer(2-nHiMin:il+nHiMax, &
                           2-nHjMin:jl+nHjMax, &
                           2-nHkMin:kl+nHkMax), stat=ierr)
           if(ierr /= 0)                       &
             call terminate("readRestartFile", &
                            "Memory allocation failure for buffer")

           if(location == Vertex) then
             allocate(bufferVertex(1:il,1:jl,1:kl), stat=ierr)
             if(ierr /= 0)                       &
               call terminate("readRestartFile", &
                          "Memory allocation failure for bufferVertex")
           endif

           ! Create a sorted version of the variable names and store the
           ! corresponding type in varTypes.

           call getSortedVarNumbers

           ! Read the density and the turbulence variables.

           call readDensity(nTypeMismatch)
           call readTurbvar(nTypeMismatch)

           ! Read the other variables, depending on the situation.

           testPrim: if(solID == 1 .or. equationMode == timeSpectral) then

             ! Either the first solution or time spectral mode. Read
             ! the primitive variables from the restart file.

             call readXvelocity(nTypeMismatch)
             call readYvelocity(nTypeMismatch)
             call readZvelocity(nTypeMismatch)
             call readPressure(nTypeMismatch)

           else testPrim

             ! Old solution in unsteady mode. Read the conservative
             ! variables.

             call readXmomentum(nTypeMismatch)
             call readYmomentum(nTypeMismatch)
             call readZmomentum(nTypeMismatch)
             call readEnergy(nTypeMismatch)

           endif testPrim

           ! Release the memory of buffer, varNames and varTypes.

           deallocate(buffer, varNames, varTypes, stat=ierr)
           if(ierr /= 0)                       &
             call terminate("readRestartFile", &
                            "Deallocation error for buffer, varNames &
                            &and varTypes.")

           ! In case bufferVertex is allocated, release it.

           if(location == Vertex) then
             deallocate(bufferVertex, stat=ierr)
             if(ierr /= 0)                       &
               call terminate("readRestartFile", &
                              "Deallocation error for bufferVertex")
           endif

         enddo domains

         ! Release the memory of zoneNames and zoneNumbers.

         deallocate(zoneNames, zoneNumbers, stat=ierr)
         if(ierr /= 0)                       &
           call terminate("readRestartFile", &
                          "Deallocation failure for zoneNames &
                          &and zoneNumbers.")

         ! Close the cgns solution file.

         call cg_close_f(cgnsInd, ierr)
         if(ierr /= all_ok)                  &
           call terminate("readRestartFile", &
                          "Something wrong when calling cg_close_f")

       enddo solLoop

       ! Write a message about the time step number for which is
       ! restarted.

       if(equationMode == unsteady .and. myID == 0) then
         write(integerString,"(i7)") nTimeStepsRestart+1
         integerString = adjustl(integerString)

         print "(a)", "#"
         print 110, trim(integerString)
         print "(a)", "#"
 110     format("# Restarting at time step",1X,A,".")
       endif

       ! Determine the global sum of nTypeMismatch; the result only
       ! needs to be known on processor 0. Use ii as the global buffer
       ! to store the result. If a type mismatch occured,
       ! print a warning.

       call mpi_reduce(nTypeMismatch, ii, 1, sumb_integer, &
                       mpi_sum, 0, SUmb_comm_world, ierr)
       if(myID == 0 .and. ii > 0) then

         write(integerString,"(i6)") ii
         integerString = adjustl(integerString)

         print "(a)", "#"
         print "(a)", "#                      Warning"
         print 120, trim(integerString)
         print "(a)", "#"
 120     format("# ",a," type mismatches occured when reading the &
                &solution of the blocks")
       endif

#endif

       end subroutine readRestartFile
