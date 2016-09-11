module writeCGNSVolume

contains

  subroutine writeCGNSVolumeSol
    !
    !       writeCGNSVolumeSol and its subroutines write the cell          
    !       centered CGNS solution file(s).                                
    !
    use block
    use cgnsGrid
    use communication
    use inputPhysics
    use IOModule
    use su_cgns
    use outputMod
    use inputIteration
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nn
    integer(kind=intType) :: nVolSolvar, nVolDiscrVar

    character(len=maxCGNSNameLen), &
         dimension(:), allocatable :: solNames

    ! Determine the number and names of the solution files.
    ! Also set the pointers for IOVar needed for the general
    ! treatment of the IO.

    call volSolFileNamesWrite

    ! Return immediately if no volume solution files must be written.

    if(nVolSolToWrite == 0) return

    ! Write a message that the solution file(s) are being written.
    ! Of course only processor 0 does this.

    if(myID == 0 .and. printIterations) then
       print "(a)", "#"
       print "(a,a)", "# Writing volume solution file(s): ",trim(volSolFileNames(1))
    endif

    ! Open the CGNS file(s), the convergence info and if needed the
    ! time accurate data for an unsteady computation. This is only
    ! done by processor 0.

    if(myID == 0) then
       call openCGNSVolumeSol
       call writeCGNSConvInfo

       if(equationMode == unsteady) call writeCGNSTimeHistory
    endif

    ! Determine the number of variables to be written to the volume
    ! solution file(s) as well as the CGNS names.

    call numberOfVolSolVariables(nVolSolvar, nVolDiscrVar)
    allocate(solNames(nVolSolvar+nVolDiscrVar), stat=ierr)
    if(ierr /= 0)                          &
         call terminate("writeCGNSVolumeSol", &
         "Memory allocation failure for solNames")
    call volSolNames(solNames)

    ! Loop over the number of CGNS blocks and write the cell centered
    ! solution(s) of this block.

    do nn=1,cgnsNDom
       call writeSolCGNSZone(nn, nVolSolvar, nVolDiscrVar, solNames)
    enddo

    ! Close the cgns file(s). Only processor 0 does this.

    if(myID == 0) then
       do nn=1,nVolSolToWrite
          call cg_close_f(fileIDs(nn), ierr)
          if(ierr /= CG_OK)                     &
               call terminate("writeCGNSVolumeSol", &
               "Something wrong when calling cg_close_f")
       enddo

    end if

    ! Deallocate the memory of fileIDs and cgnsBases.  These are
    ! allocated ALL PROCESSORS not just processor 0.
    ! Fixed Bug: GKK 

    if (allocated(fileIDs)) then 
       deallocate(fileIDs, stat=ierr)
    end if
    if (allocated(cgnsBases)) then
       deallocate(cgnsBases, stat=ierr)
    end if

    if(ierr /= 0)                          &
         call terminate("writeCGNSVolumeSol", &
         "Deallocation error for fileIDs &
         &and cgnsBases.")

    ! Deallocate the memory of solNames.

    deallocate(solNames, stat=ierr)
    if(ierr /= 0)                          &
         call terminate("writeCGNSVolumeSol", &
         "Deallocation error for solNames.")

    ! Deallocate the memory of IOVar. Note that the first entry
    ! is used as a temporary buffer.

    do nn=1,nDom
       deallocate(IOVar(nn,1)%w, stat=ierr)
       if(ierr /= 0)                          &
            call terminate("writeCGNSVolumeSol", &
            "Deallocation error for IOVar%w")
    enddo

    deallocate(IOVar, stat=ierr)
    if(ierr /= 0)                          &
         call terminate("writeCGNSVolumeSol", &
         "Deallocation error for IOVar")

    ! Wait until all processors (especially processor 0) reach
    ! this point.

    call mpi_barrier(SUmb_comm_world, ierr)

    ! Write a message that the solution file has been written.
    ! Of course only processor 0 does this.

    if(myID == 0 .and. printIterations) then
       print "(a)", "# Volume solution file(s) written"
       print "(a)", "#"
    endif
  end subroutine writeCGNSVolumeSol

  subroutine volSolFileNamesWrite
    !
    !       volSolFileNamesWrite determines the names and number of volume 
    !       solution files to be written. Furthermore it sets the pointers 
    !       and/or allocates the memory for IOVar to make a general        
    !       treatment of the writing possible.                             
    !
    use block
    use inputIO
    use inputPhysics
    use inputTimeSpectral
    use IOModule
    use iteration
    use monitor
    use outputMod
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nn, mm, kk, nAvail
    integer(kind=intType) :: iEnd, jEnd, kEnd

    character(len=7) :: intString

    !       Determine the names and number of volume solution files to be  
    !       written.                                                       
    !
    ! Determine the situation we are having here.

    select case (equationMode)

    case (steady)

       ! Steady state computation. Allocate the memory for the
       ! volume solution file name.

       allocate(volSolFileNames(1), stat=ierr)
       if(ierr /= 0)                            &
            call terminate("volSolFileNamesWrite", &
            "Memory allocation failure for &
            &volSolFileNames")

       ! Volume solution file. Always set the name.
       ! Set nVolSolToWrite to 1 if a solution file must be written;
       ! otherwise set it to 0.

       volSolFileNames(1) = solFile
       if( writeVolume ) then
          nVolSolToWrite = 1
       else
          nVolSolToWrite = 0
       endif

       !===============================================================

    case (unsteady)

       ! Unsteady computation. For a consistent restart nOldLevels
       ! solutions must be written. However, it is possible that not
       ! that many solutions are available or that some of these
       ! solutions have already been written in an earlier time step.

       ! First determine the number of available solutions.

       nAvail = timeStepUnsteady + nTimeStepsRestart + 1
       nAvail = min(nAvail,nOldLevels)

       ! Allocate the memory for the file names. Note that this is
       ! an upper boundary. It is possible that less files need
       ! to be written.

       allocate(volSolFileNames(nAvail), stat=ierr)
       if(ierr /= 0)                            &
            call terminate("volSolFileNamesWrite", &
            "Memory allocation failure for &
            &volSolFileNames")

       ! Set the names of the files.

       do nn=1,nAvail
          write(intString,"(i7)") timeStepUnsteady + &
               nTimeStepsRestart + 1 - nn
          intString = adjustl(intString)

          volSolFileNames(nn)  = trim(solfile)//"&
               &Timestep"//trim(intString)
       enddo

       ! Determine the number of volume solution files to write.

       if( writeVolume ) then

          ! Initialize nVolSolToWrite to 1.

          nVolSolToWrite = 1

          ! Loop over the older levels and check if some of
          ! them must be written as well.

          do nn=1,(nAvail-1)
             if(.not. oldSolWritten(nn) ) then
                nVolSolToWrite = nVolSolToWrite + 1
                volSolFileNames(nVolSolToWrite) = volSolFileNames(nn+1)
             endif
          enddo

       else

          ! No volume solution files need to be written.

          nVolSolToWrite = 0

       endif

       !===============================================================

    case (timeSpectral)

       ! Time spectral computation. Allocate the file names.

       allocate(volSolFileNames(nTimeIntervalsSpectral), stat=ierr)
       if(ierr /= 0)                            &
            call terminate("volSolFileNamesWrite", &
            "Memory allocation failure for &
            &volSolFileNames")

       ! Set the names of the files.

       do nn=1,nTimeIntervalsSpectral
          write(intString,"(i7)") nn
          intString = adjustl(intString)

          volSolFileNames(nn)  = trim(solfile)//"&
               &Spectral"//trim(intString)
       enddo

       ! Set the number of volume solution files to write.
       ! Either they are written or they are not written.

       if( writeVolume ) then
          nVolSolToWrite = nTimeIntervalsSpectral
       else
          nVolSolToWrite = 0
       endif

    end select
    !
    !       Set the pointers for IOVar if volume solution files need to be 
    !       written.                                                       
    !
    testSolsToWrite: if(nVolSolToWrite > 0) then

       ! Allocate the memory for IOVar.

       allocate(IOVar(nDom,nVolSolToWrite), stat=ierr)
       if(ierr /= 0)                            &
            call terminate("volSolFileNamesWrite", &
            "Memory allocation failure for IOVar")

       ! As the writing normally involves other variables than just
       ! the primitive ones, memory for the member variable w must be
       ! allocated to make the general IO treatment possible. This is
       ! a bit of an overhead, but that's a small price to pay for the
       ! general treatment.

       if( storeRindLayer ) then
          do nn=1,nDom
             IOVar(nn,1)%pointerOffset = 0

             iEnd = flowDoms(nn,1,1)%ie
             jEnd = flowDoms(nn,1,1)%je
             kEnd = flowDoms(nn,1,1)%ke

             allocate(IOVar(nn,1)%w(iEnd,jEnd,kEnd,1), stat=ierr)
             if(ierr /= 0)                            &
                  call terminate("volSolFileNamesWrite", &
                  "Memory allocation failure for IOVar%w")
          enddo
       else
          do nn=1,nDom
             IOVar(nn,1)%pointerOffset = 0

             iEnd = flowDoms(nn,1,1)%il
             jEnd = flowDoms(nn,1,1)%jl
             kEnd = flowDoms(nn,1,1)%kl

             allocate(IOVar(nn,1)%w(2:iEnd,2:jEnd,2:kEnd,1), stat=ierr)
             if(ierr /= 0)                            &
                  call terminate("volSolFileNamesWrite", &
                  "Memory allocation failure for IOVar%w")
          enddo
       endif

       ! Set the pointers for the other solutions depending on the
       ! situation.

       select case(equationMode)

       case (steady, timeSpectral)

          ! Actually only time spectral mode, but steady is added to
          ! avoid a compiler warning. Set the pointers for the higher
          ! spectral solution to the first solution.

          do mm=2,nVolSolToWrite
             do nn=1,nDom
                IOVar(nn,mm)%pointerOffset = 0
                IOVar(nn,mm)%w => IOVar(nn,1)%w 
             enddo
          enddo

          !=============================================================

       case (unsteady)

          ! It is possible that for an unsteady computation previous
          ! solutions need to be written. However only the variables
          ! wOld need to be written, so the pointer can be set to the
          ! correct entries. As the starting indices of wOld are 2, 
          ! a pointer shift takes place here. I know this is a pain
          ! in the butt, but that's what we have to live with.

          kk = 1
          do mm=1,(nAvail-1)
             if(.not. oldSolWritten(mm) ) then
                kk = kk + 1
                do nn=1,nDom
                   IOVar(nn,kk)%pointerOffset = -1
                   IOVar(nn,kk)%w => flowDoms(nn,1,1)%wOld(mm,2:,2:,2:,:)
                enddo
             endif
          enddo

       end select

    endif testSolsToWrite

  end subroutine volSolFileNamesWrite

  subroutine openCGNSVolumeSol
    !
    !       openCGNSVolumeSol opens the cgns solution file(s) if needed.   
    !       If opened the files are opened either for writing or for       
    !       modification. When the grid file(s) have been written, these   
    !       files are still open and nothing needs to be done.             
    !       Only processor 0 performs this task.                           
    !
    use cgnsGrid
    use monitor
    use su_cgns
    use outputMod
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nn

    character(len=maxStringLen) :: errorMessage

    ! Check if grid files have been written and if the solution must
    ! be stored in the same file. If so only the CGNS header must
    ! be written. Do this and return immediately.

    if(writeGrid .and. (.not. useLinksInCGNS)) then
       do nn=1,nVolSolToWrite
          call writeCGNSHeader(fileIDs(nn), cgnsBases(nn))
       enddo

       return
    endif

    ! Solution files must be created. Allocate the memory for
    ! the fileIDs and the bases.

    allocate(fileIDs(nVolSolToWrite), cgnsBases(nVolSolToWrite), &
         stat=ierr)
    if(ierr /= 0)                         &
         call terminate("openCGNSVolumeSol", &
         "Memory allocation failure for fileIDs and &
         &cgnsBases")

    ! Determine the situation we are having here.

    testForLinks: if( useLinksInCGNS ) then

       ! Links are used to the coordinates of the zones. This means
       ! that the files must be opened in write mode.

       do nn=1,nVolSolToWrite
          call cg_open_f(volSolFileNames(nn), mode_write, &
               fileIDs(nn), ierr)
          if(ierr /= CG_OK) then
             write(errorMessage,*) "File ", trim(volSolFileNames(nn)), &
                  " could not be opened by CGNS&
                  & for writing"
             call terminate("openCGNSVolumeSol", errorMessage)
          endif

          ! Create the base.

          call cg_base_write_f(fileIDs(nn), cgnsBaseName, cgnsCelldim, &
               cgnsPhysdim, cgnsBases(nn), ierr)
          if(ierr /= CG_OK)                    &
               call terminate("openCGNSVolumeSol", &
               "Something wrong when calling &
               &cg_base_write_f")
       enddo

    else testForLinks

       ! Solutions must be written in the same file(s) as the grid.
       ! As the grid file(s) are not written during the current call
       ! to writeSol, these are old files and must therefore be opened
       ! in modify mode.

       do nn=1,nVolSolToWrite
          call cg_open_f(volSolFileNames(nn), mode_modify, &
               fileIDs(nn), ierr)
          if(ierr /= CG_OK) then
             write(errorMessage,*) "File ", trim(volSolFileNames(nn)), &
                  " could not be opened by CGNS&
                  & for writing"
             call terminate("openCGNSVolumeSol", errorMessage)
          endif

          ! Simply set the base IDs to 1.

          cgnsBases(nn) = 1
       enddo

    endif testForLinks

    ! Write the CGNS header.

    do nn=1,nVolSolToWrite
       call writeCGNSHeader(fileIDs(nn), cgnsBases(nn))
    enddo
  end subroutine openCGNSVolumeSol

  subroutine writeCGNSConvInfo
    !
    !       writeCGNSConvInfo writes the convergence info to the           
    !       cgns file(s).                                                  
    !
    use inputIO
    use inputPhysics
    use monitor
    use su_cgns
    use outputMod
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: nConv, base, cgnsInd, conv
    integer :: i, nn, mm, ierr, realTypeCGNS

    real(kind=8), dimension(:), allocatable :: buf8

    ! Return immediately if the convergence history (of the inner
    ! iterations) does not need to be stored. This logical can
    ! only be .false. for an unsteady computation.

    if(.not. storeConvInnerIter) return

    ! Store the number of iterations to be written in nn.
    ! This is nIterCur + 1, because the array starts at 0.

    nn = nIterCur + 1

    ! Depending on the input option, set the CGNS type and allocate
    ! the memory for either buf4 or buf8.

    ! Always write the cgnsConvergence history to be double since
    ! that's what Tecplot needs

    realTypeCGNS = RealDouble
    allocate(buf8(0:nIterCur), stat=ierr)


    if(ierr /= 0)                         &
         call terminate("writeCGNSConvInfo", &
         "Memory allocation failure for either buf4 &
         &or buf8")

    ! Determine the number of convergence histories to be written.
    ! This depends on the equation mode.

    select case (equationMode)

    case (steady,unsteady)
       nConv = 1

    case (timeSpectral)
       nConv = nVolSolToWrite    ! == number of spectral solutions.

    end select

    ! Loop over the number of convergence histories to be written.

    convLoop: do conv=1,nConv

       ! Abbreviate the corresponding file and base a bit easier.

       cgnsInd = fileIDs(conv)
       base    = cgnsBases(conv)

       ! Go to the correct position in the CGNS file.

       call cg_goto_f(cgnsInd, base, ierr, "end")
       if(ierr /= CG_OK)                    &
            call terminate("writeCGNSConvInfo", &
            "Something wrong when calling cg_goto_f")

       ! Create the convergence history node. Add a small description.

       call cg_convergence_write_f(nn,"L2 norms are computed by taking &
            &the square root of the quotient &
            &the sum of the square of the &
            &residuals and the total number &
            &of cells in the grid.", ierr)
       if(ierr /= CG_OK)                    &
            call terminate("writeCGNSConvInfo", &
            "Something wrong when calling &
            &cg_convergence_write_f")

       ! The convergence history must be written under the node just
       ! created. Go there.

       call cg_goto_f(cgnsInd, base, ierr, &
            "ConvergenceHistory_t", 1, "end")
       if(ierr /= CG_OK)                    &
            call terminate("writeCGNSConvInfo", &
            "Something wrong when calling cg_goto_f")

       ! Loop over the number of monitoring variables.

       monLoop: do i=1,nMon

          ! Copy the convergence info to either buf4 or buf8 and write
          ! it to file.
          do mm=0,nIterCur
             buf8(mm) = convArray(mm,conv,i)
          enddo

          call cg_array_write_f(monNames(i), realTypeCGNS, 1, nn, &
               buf8, ierr)
          if(ierr /= CG_OK)                    &
               call terminate("writeCGNSConvInfo", &
               "Something wrong when calling &
               &cg_array_write_f")
       enddo monLoop
    enddo convLoop

    ! Release the memory of buf8.

    deallocate(buf8, stat=ierr)

    if(ierr /= 0)                         &
         call terminate("writeCGNSConvInfo", &
         "Deallocation failure for either buf4 or buf8")
  end subroutine writeCGNSConvInfo

  subroutine writeCGNSTimeHistory
    !
    !       WriteCGNSTimeHistory writes for unsteady computations          
    !       the time history of the monitoring variables to the            
    !       cgns file.                                                     
    !
    use cgnsNames
    use inputIO
    use monitor
    use su_cgns
    use outputMod
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: i, nn, mm, cgnsInd, base, ierr, realTypeCGNS

    real(kind=4), dimension(:), allocatable :: buf4
    real(kind=8), dimension(:), allocatable :: buf8

    character(len=maxCGNSNameLen) :: cgnsName

    ! Store the file index and base a bit easier. Note that this info
    ! is only written to the 1st file, because this is an unsteady
    ! computation.

    cgnsInd = fileIDs(1)
    base    = cgnsBases(1)

    ! Store the total number of time steps in nn.

    nn = timeStepUnsteady + nTimeStepsRestart

    ! Depending on the input option, set the CGNS type and allocate
    ! the memory for either buf4 or buf8.

    ! Set the cgns real type depending on the input option.

    select case (precisionSol)
    case (precisionSingle)
       realTypeCGNS = RealSingle
       allocate(buf4(nn), stat=ierr)

       !===============================================================

    case (precisionDouble)
       realTypeCGNS = RealDouble
       allocate(buf8(nn), stat=ierr)
    end select

    if(ierr /= 0)                            &
         call terminate("writeCGNSTimeHistory", &
         "Memory allocation failure for either buf4 &
         &or buf8")

    ! Go to the correct position in the cgns file.

    call cg_goto_f(cgnsInd, base, ierr, "end")
    if(ierr /= CG_OK)                       &
         call terminate("writeCGNSTimeHistory", &
         "Something wrong when calling cg_goto_f")

    ! Create the name of the base iterative data node.

    cgnsName = "TimeHistory"

    ! Create the base iterative node and check if everything
    ! went okay.

    call cg_biter_write_f(cgnsInd, base, cgnsName, nn, ierr)
    if(ierr /= CG_OK)                       &
         call terminate("writeCGNSTimeHistory", &
         "Something wrong when calling cg_biter_write_f")

    ! The time history must be written under the node just created.
    ! Go there.

    call cg_goto_f(cgnsInd, base, ierr, &
         "BaseIterativeData_t" , 1, "end")
    if(ierr /= CG_OK)                          &
         call terminate("writeCGNSTimeHistory", &
         "Something wrong when calling cg_goto_f")

    ! Write the time values.

    cgnsName = cgnsTimeValue
    call cg_array_write_f(cgnsName, realTypeCGNS, 1, nn, &
         timeArray, ierr)
    if(ierr /= CG_OK)                       &
         call terminate("writeCGNSTimeHistory", &
         "Something wrong when calling cg_array_write_f")

    ! Loop over the number of monitoring variables and write
    ! their time history.

    monLoop: do i=1,nMon

       ! Copy the time history to either buf4 or buf8 and write it
       ! to file.

       select case (precisionSol)
       case (precisionSingle)
          do mm=1,nn
             buf4(mm) = timeDataArray(mm,i)
          enddo

          call cg_array_write_f(monNames(i), realTypeCGNS, 1, nn, &
               buf4, ierr)

          !=============================================================

       case (precisionDouble)
          do mm=1,nn
             buf8(mm) = timeDataArray(mm,i)
          enddo

          call cg_array_write_f(monNames(i), realTypeCGNS, 1, nn, &
               buf8, ierr)
       end select

       if(ierr /= CG_OK)                       &
            call terminate("writeCGNSTimeHistory", &
            "Something wrong when calling &
            &cg_array_write_f")
    enddo monLoop

    ! Release the memory of buf4 or buf8.

    select case (precisionSol)
    case (precisionSingle)
       deallocate(buf4, stat=ierr)

    case (precisionDouble)
       deallocate(buf8, stat=ierr)
    end select

    if(ierr /= 0)                            &
         call terminate("writeCGNSTimeHistory", &
         "Deallocation failure for either buf4 or buf8") 
  end subroutine writeCGNSTimeHistory

  subroutine writeSolCGNSZone(zone, nSolVar, nDiscrVar, solNames)
    !
    !       writeSolCGNSZone writes a volume solution of the given zone    
    !       to the cgns file(s). In case the solution must be written to a 
    !       separate file, useLinksInCGNS == .true., a link to the zone of 
    !       the grid file is created.                                      
    !
    use blockPointers
    use cgnsGrid
    use cgnsNames
    use communication
    use flowVarRefState
    use inputIO
    use inputPhysics
    use iteration
    use su_cgns
    use outputMod
    use utils, only : terminate, setPointers
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: zone
    integer(kind=intType), intent(in) :: nSolVar, nDiscrVar

    character(len=*), dimension(*), intent(in) :: solNames
    !
    !      Local variables.
    !
    integer :: ierr
    integer :: source, bufSize, size, nnVar
    integer :: cgnsInd, cgnsBase, cgnsZone, cgnsSol, realTypeCGNS

    integer, dimension(mpi_status_size) :: status
    integer, dimension(9)               :: sizes

    integer(kind=intType) :: i, j, nn, mm, ll, ind, nVarWritten
    integer(kind=intType) :: nBlocks, nSubblocks, offset
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd
    integer(kind=intType) :: iBegCGNS, jBegCGNS, kBegCGNS
    integer(kind=intType) :: iEndCGNS, jEndCGNS, kEndCGNS

    integer(kind=intType), dimension(nProc) :: nMessages

    integer(kind=intType), dimension(:), allocatable :: proc

    integer(kind=intType), dimension(:,:,:), allocatable :: subRanges

    real(kind=realType), dimension(:), allocatable :: buffer

    real(kind=4), dimension(:), allocatable :: sol4
    real(kind=8), dimension(:), allocatable :: sol8

    logical :: rindLayerThisSol, unsteadyHigherSol, writeLink

    character(len=maxStringLen) :: linkName, solName

    ! Set the cgns real type depending on the input option.

    ! Store the number of local blocks and the offset in
    ! blocksCGNSblock for this zone a bit easier.

    offset  = nBlocksCGNSblock(zone-1)
    nBlocks = nBlocksCGNSblock(zone) - offset

    ! Determine the amount of block parts each processor will send to
    ! processor 0.

    call mpi_gather(nBlocks, 1, sumb_integer, nMessages, 1, &
         sumb_integer, 0, SUmb_comm_world, ierr)

    ! At the moment the writing of the cgns file is sequential and done
    ! by processor 0. This means that this processor gathers all info
    ! from the other processors and writes it to file.

    rootproc: if(myID == 0) then
       !
       !         I am processor 0 and poor me has to do all the work.         
       !
       ! Allocate the memory for the array used to write the solution
       ! to file. The size depends whether or not rind layers are to
       ! be written and the precision of the floating point type.

       if( storeRindLayer ) then
          ll = (cgnsDoms(zone)%kl+1) * (cgnsDoms(zone)%jl+1) &
               * (cgnsDoms(zone)%il+1) 
       else
          ll = (cgnsDoms(zone)%kl-1) * (cgnsDoms(zone)%jl-1) &
               * (cgnsDoms(zone)%il-1) 
       endif

       select case(precisionSol)
       case (precisionSingle)
          allocate(sol4(ll), sol8(0), stat=ierr)
       case (precisionDouble)
          allocate(sol4(0), sol8(ll), stat=ierr)
       end select

       if(ierr /= 0)                        &
            call terminate("writeSolCGNSZone", &
            "Memory allocation failure for sol")

       ! First determine the number of subblocks into the original cgns
       ! block is split.

       nSubblocks = 0
       do i=1,nProc
          nSubblocks = nSubblocks + nMessages(i)
       enddo

       ! Allocate the memory for the ranges and the processor
       ! where the subblock is stored.

       allocate(subRanges(3,2,nSubblocks), proc(nSubblocks), stat=ierr)
       if(ierr /= 0)                        &
            call terminate("writeSolCGNSZone", &
            "Memory allocation failure for subRanges &
            &and proc")

       ! Determine the processor ID's where the subRanges are stored.
       ! Note that 1 must be substracted, because the processor numbering
       ! starts at 0.

       nSubblocks = 0
       do i=1,nProc
          do j=1,nMessages(i)
             nSubblocks = nSubblocks + 1
             proc(nSubblocks) = i - 1
          enddo
       enddo

       ! Determine the subranges for the 1st solution.

       rindLayerThisSol = storeRindLayer
       call getSubRangesSol

       ! Allocate the memory for buffer.

       allocate(buffer(bufSize), stat=ierr)
       if(ierr /= 0)                        &
            call terminate("writeSolCGNSZone", &
            "Memory allocation failure for buffer")

       ! Loop over the number of solutions.

       solLoopRoot: do ind=1,nVolSolToWrite

          ! Determine whether or not we are dealing with an unsteady
          ! higher solution here.

          unsteadyHigherSol = .false.
          if(ind > 1 .and. equationMode == unsteady) &
               unsteadyHigherSol = .true.

          ! For unsteady mode on rigid meshes it is possible that
          ! only the solution must be written; the coordinates are
          ! not written to a file. In case links are used for that
          ! case the link should not be created. The logical
          ! writeLink takes care of that.

          writeLink = useLinksInCGNS
          if(unsteadyHigherSol .and. (.not. deforming_Grid)) &
               writeLink = .false.

          ! Store the file and base ID a bit easier and set
          ! rindLayerThisSol. A rind layer is not written for the
          ! higher solution in unsteady mode.

          cgnsInd          = fileIDs(ind)
          cgnsBase         = cgnsBases(ind)
          rindLayerThisSol = storeRindLayer
          if( unsteadyHigherSol ) rindLayerThisSol = .false.

          ! Check if the subranges must be recomputed. Some dirty stuff
          ! must be done, because Fortran 90/95 does not allow a
          ! comparison between logicals.

          nn = 0; if( storeRindLayer )   nn = 1
          mm = 0; if( rindLayerThisSol ) mm = 1
          if(nn /= mm) call getSubRangesSol

          ! Check whether a zone must be created, which is
          ! true if links are used. Create the zone, if needed.

          createZoneTest: if( useLinksInCGNS ) then

             ! A zone must be created. Use the same name as
             ! in the original grid file.

             sizes(1) = cgnsDoms(zone)%il
             sizes(2) = cgnsDoms(zone)%jl
             sizes(3) = cgnsDoms(zone)%kl
             sizes(4) = cgnsDoms(zone)%nx
             sizes(5) = cgnsDoms(zone)%ny
             sizes(6) = cgnsDoms(zone)%nz
             sizes(7) = 0
             sizes(8) = 0
             sizes(9) = 0

             call cg_zone_write_f(cgnsInd, cgnsBase,              &
                  cgnsDoms(zone)%zonename, sizes, &
                  Structured, cgnsZone, ierr)
             if(ierr /= CG_OK)                 then

                call terminate("writeSolCGNSZone", &
                     "Something wrong when calling &
                     &cg_zone_write_f")
             end if
             ! Check if a link must actually be created.

             writeLinkTest: if( writeLink ) then

                ! Create the link of the coordinates to the zone in the
                ! original grid.  First move to the correct location.

                call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", &
                     cgnsZone, "end")
                if(ierr /= CG_OK)                   &
                     call terminate("writeSolCGNSZone", &
                     "Something wrong when calling cg_goto_f")

                ! Determine the link name and write the link to the grid
                ! coordinates.

                linkName = cgnsBasename//"&
                     &/"//cgnsDoms(zone)%zonename//"&
                     &/"//"GridCoordinates"

                call cg_link_write_f("GridCoordinates", &
                     gridFileNames(ind), linkName, ierr)
                if(ierr /= CG_OK)                   &
                     call terminate("writeSolCGNSZone", &
                     "Something wrong when calling &
                     &cg_link_write_f")
             endif writeLinkTest

          else createZoneTest

             ! The zone already exists. Simply set cgnsZone to zone.

             cgnsZone = zone

          endif createZoneTest

          ! Create the flow solution node.

          call cg_sol_write_f(cgnsInd, cgnsBase, cgnsZone,          &
               "Flow solution", CellCenter, cgnsSol, &
               ierr)
          if(ierr /= CG_OK)                   &
               call terminate("writeSolCGNSZone", &
               "Something wrong when calling &
               &cg_sol_write_f")

          ! Create the rind layers. If rind layers must be stored put
          ! 1 layer on every face of the block; otherwise put 0 layers.
          ! Use sizes as a buffer to store the rind data. The rind data
          ! must be created under the just created solution node.

          call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", &
               cgnsZone, "FlowSolution_t", cgnsSol, "end")
          if(ierr /= CG_OK)                   &
               call terminate("writeSolCGNSZone", &
               "Something wrong when calling cg_goto_f")

          if( rindLayerThisSol ) then
             sizes(1) = 1; sizes(2) = 1; sizes(3) = 1
             sizes(4) = 1; sizes(5) = 1; sizes(6) = 1
          else
             sizes(1) = 0; sizes(2) = 0; sizes(3) = 0
             sizes(4) = 0; sizes(5) = 0; sizes(6) = 0
          endif

          call cg_rind_write_f(sizes, ierr)
          if(ierr /= CG_OK)                   &
               call terminate("writeSolCGNSZone", &
               "Something wrong when calling &
               &cg_rind_write_f")

          ! Determine the index range of the solution of the zone
          ! to be written.

          iBegCGNS = 2 - sizes(1)
          jBegCGNS = 2 - sizes(3)
          kBegCGNS = 2 - sizes(5)

          iEndCGNS = cgnsDoms(zone)%il + sizes(2)
          jEndCGNS = cgnsDoms(zone)%jl + sizes(4)
          kEndCGNS = cgnsDoms(zone)%kl + sizes(6)

          ! Determine the number of variables to be written.
          ! For the unsteady higher solutions only the variables
          ! needed for a restart are written.

          nVarWritten = nSolVar+nDiscrVar
          if( unsteadyHigherSol ) nVarWritten = nw

          ! Loop over the number of variables to be written.

          varWriteLoop: do nn=1,nVarWritten

             ! Copy solNames(nn) in solName for later purposes.
             ! Correct this value if the conservative variables
             ! for a consistent unsteady restart must be written.
             ! No need to correct the turbulent variables, because
             ! these names are okay.

             solName = solNames(nn)

             if( unsteadyHigherSol ) then
                nnVar = nn

                select case(nnVar)
                case (irho)
                   solName = cgnsDensity

                case (imx)
                   solName = cgnsMomx

                case (imy)
                   solName = cgnsMomy

                case (imz)
                   solName = cgnsMomz

                case (irhoE)
                   solName = cgnsEnergy
                end select
             endif

             ! Loop over the number of subblocks stored on
             ! this processor.

             do mm=1,nBlocks

                ! Set the pointers to the local domain.

                if( unsteadyHigherSol ) then
                   call setPointers(blocksCGNSblock(mm+offset), &
                        1_intType, 1_intType)
                else
                   call setPointers(blocksCGNSblock(mm+offset), &
                        1_intType, ind)
                endif

                ! Determine the cell range I have to write. This depends
                ! whether or not halo's must be written.

                iBeg = 2; iEnd = il
                jBeg = 2; jEnd = jl
                kBeg = 2; kEnd = kl

                if(storeRindLayer .and. (.not. unsteadyHigherSol)) then
                   if(iBegor == 1) iBeg = 1
                   if(jBegor == 1) jBeg = 1
                   if(kBegor == 1) kBeg = 1

                   if(iEndor == cgnsDoms(zone)%il) iEnd = ie
                   if(jEndor == cgnsDoms(zone)%jl) jEnd = je
                   if(kEndor == cgnsDoms(zone)%kl) kEnd = ke
                endif

                ! Fill the buffer with the correct solution variable.
                ! The routine called depends on the situation.

                if( unsteadyHigherSol ) then

                   ! Unsteady higher solution. Write the old solution.

                   call storeOldSolInBuffer(buffer, ind, nn, iBeg, iEnd, &
                        jBeg, jEnd, kBeg, kEnd)
                else

                   ! Standard solution must be written.

                   call storeSolInBuffer(buffer, .true., solName, &
                        iBeg, iEnd, jBeg, jEnd,  &
                        kBeg, kEnd)
                endif

                ! And store it in sol. The routine called depends on
                ! the desired precision.

                select case (precisionSol)
                case (precisionSingle)
                   call copyDataBufSinglePrecision(sol4, buffer,        &
                        iBegCGNS, jBegCGNS, &
                        kBegCGNS, iEndCGNS, &
                        jEndCGNS, kEndCGNS, &
                        subRanges(1,1,mm))
                case (precisionDouble)
                   call copyDataBufDoublePrecision(sol8, buffer,        &
                        iBegCGNS, jBegCGNS, &
                        kBegCGNS, iEndCGNS, &
                        jEndCGNS, kEndCGNS, &
                        subRanges(1,1,mm))
                end select

             enddo

             ! Loop over the number of subblocks stored on
             ! other processors.

             do mm=(nBlocks+1),nSubblocks

                ! Receive the range of subblock mm.

                source = proc(mm)
                call mpi_recv(buffer, bufSize, sumb_real, source, &
                     source+1, SUmb_comm_world, status, ierr)

                ! And store it in sol.

                select case (precisionSol)
                case (precisionSingle)
                   call copyDataBufSinglePrecision(sol4, buffer,        &
                        iBegCGNS, jBegCGNS, &
                        kBegCGNS, iEndCGNS, &
                        jEndCGNS, kEndCGNS, &
                        subRanges(1,1,mm))
                case (precisionDouble)
                   call copyDataBufDoublePrecision(sol8, buffer,        &
                        iBegCGNS, jBegCGNS, &
                        kBegCGNS, iEndCGNS, &
                        jEndCGNS, kEndCGNS, &
                        subRanges(1,1,mm))
                end select

             enddo

             ! Write the solution variable to file. Source is just used
             ! as a dummy variable and does not have a meaning.
             select case(precisionSol)
             case (precisionSingle)
                call cg_field_write_f(cgnsInd, cgnsBase, cgnsZone, &
                     cgnsSol, realSingle, solName, sol4, source, ierr)
             case (precisionDouble)
                call cg_field_write_f(cgnsInd, cgnsBase, cgnsZone, &
                     cgnsSol, realDouble, solName, sol8, source, ierr)
             end select
             if(ierr /= CG_OK)    &
                  call terminate("writeSolCGNSZone", &
                  "Something wrong when calling &
                  &cg_field_write_f")
          enddo varWriteLoop

       enddo solLoopRoot

       ! Release some memory only allocated on the root processor.
       deallocate(sol4, sol8, subRanges, proc, stat=ierr)
       if(ierr /= 0) call terminate("writeSolCGNSZone", &
            "Deallocation error on root proc")

    else rootproc
       !
       !         I am not the root processor and may have to send some data   
       !         to the root processor.                                       
       ! 
       ! Determine the subranges for the 1st solution.

       rindLayerThisSol = storeRindLayer
       call getSubRangesSol

       ! Allocate the memory for buffer.

       allocate(buffer(bufSize), stat=ierr)
       if(ierr /= 0)                        &
            call terminate("writeSolCGNSZone", &
            "Memory allocation failure for buffer")

       ! Loop over the number of solutions.

       solLoopOthers: do ind=1,nVolSolToWrite

          ! Determine whether or not we are dealing with an unsteady
          ! higher solution here.

          unsteadyHigherSol = .false.
          if(ind > 1 .and. equationMode == unsteady) &
               unsteadyHigherSol = .true.

          ! Set rindLayerThisSol. A rind layer is not written for
          ! the higher solutions in unsteady mode.

          rindLayerThisSol = storeRindLayer
          if( unsteadyHigherSol ) rindLayerThisSol = .false.

          ! Check if the subranges must be recomputed. Some dirty stuff
          ! must be done, because Fortran 90/95 does not allow a
          ! comparison between logicals.

          nn = 0; if( storeRindLayer )   nn = 1
          mm = 0; if( rindLayerThisSol ) mm = 1
          if(nn /= mm) call getSubRangesSol

          ! Determine the number of variables to be written.
          ! For the unsteady higher solutions only the variables
          ! needed for a restart are written.

          nVarWritten = nSolVar+nDiscrVar
          if( unsteadyHigherSol ) nVarWritten = nw

          ! Loop over the number of variables to be written.

          do nn=1,nVarWritten

             ! Loop over the number of subblocks stored on
             ! this processor.

             do mm=1,nBlocks

                ! Set the pointers to the local domain.

                if( unsteadyHigherSol ) then
                   call setPointers(blocksCGNSblock(mm+offset), &
                        1_intType, 1_intType)
                else
                   call setPointers(blocksCGNSblock(mm+offset), &
                        1_intType, ind)
                endif

                ! Determine the cell range I have to write. This depends
                ! whether or not halo's must be written.

                iBeg = 2; iEnd = il
                jBeg = 2; jEnd = jl
                kBeg = 2; kEnd = kl

                if(storeRindLayer .and. (.not. unsteadyHigherSol)) then
                   if(iBegor == 1) iBeg = 1
                   if(jBegor == 1) jBeg = 1
                   if(kBegor == 1) kBeg = 1

                   if(iEndor == cgnsDoms(zone)%il) iEnd = ie
                   if(jEndor == cgnsDoms(zone)%jl) jEnd = je
                   if(kEndor == cgnsDoms(zone)%kl) kEnd = ke
                endif

                ! Fill the buffer with the correct solution variable.
                ! The routine called depends on the situation.

                if( unsteadyHigherSol ) then

                   ! Unsteady higher solution. Write the old solution.

                   call storeOldSolInBuffer(buffer, ind, nn, iBeg, iEnd, &
                        jBeg, jEnd, kBeg, kEnd)
                else

                   ! Standard solution must be written.

                   call storeSolInBuffer(buffer, .true., solNames(nn), &
                        iBeg, iEnd, jBeg, jEnd,       &
                        kBeg, kEnd)
                endif

                ! And send it to processor 0.

                ll   = (iEnd-iBeg+1)*(jEnd-jBeg+1)*(kEnd-kBeg+1)
                size = ll

                call mpi_send(buffer, size, sumb_real, 0, myID+1, &
                     SUmb_comm_world, ierr)
             enddo
          enddo

       enddo solLoopOthers
    endif rootproc

    ! Release some memory.

    deallocate(buffer, stat=ierr)
    if(ierr /= 0)                        &
         call terminate("writeSolCGNSZone", &
         "Deallocation error for buffer")

    !=================================================================

  contains

    !===============================================================

    subroutine getSubRangesSol
      !
      !         getSubRangesSol determines the subranges of the              
      !         computational blocks that contribute to the CGNS block which 
      !         is currently written. Also the size of the largest subblock  
      !         is determined.                                               
      !
      implicit none
      !
      !        Local variables.
      !
      integer :: source

      integer(kind=intType) :: i, j, ll
      integer(kind=intType), dimension(6) :: ii

      ! Initialize bufSize.

      bufSize = 0

      ! Test if I'm the root processor or not.

      testRoot: if(myID == 0) then

         ! I'm the root processor. Determine the subRanges of the 
         ! subblocks stored on locally. Note that nBlocks can be 0.

         do i=1,nBlocks

            ! Store the local block ID a bit easier in j.

            j = blocksCGNSblock(i+offset)

            ! Determine the range; this is the same for all spectral
            ! solutions, so the first one can be used.

            subRanges(1,1,i) = flowDoms(j,1,1)%iBegor + 1
            subRanges(1,2,i) = flowDoms(j,1,1)%iEndor

            subRanges(2,1,i) = flowDoms(j,1,1)%jBegor + 1
            subRanges(2,2,i) = flowDoms(j,1,1)%jEndor

            subRanges(3,1,i) = flowDoms(j,1,1)%kBegor + 1
            subRanges(3,2,i) = flowDoms(j,1,1)%kEndor

            ! Correct in case rind layers must be stored.

            if( rindLayerThisSol ) then

               if(subRanges(1,1,i) == 2) subRanges(1,1,i) = 1
               if(subRanges(2,1,i) == 2) subRanges(2,1,i) = 1
               if(subRanges(3,1,i) == 2) subRanges(3,1,i) = 1

               if(subRanges(1,2,i) == cgnsDoms(zone)%il) &
                    subRanges(1,2,i) =  cgnsDoms(zone)%il + 1
               if(subRanges(2,2,i) == cgnsDoms(zone)%jl) &
                    subRanges(2,2,i) =  cgnsDoms(zone)%jl + 1
               if(subRanges(3,2,i) == cgnsDoms(zone)%kl) &
                    subRanges(3,2,i) =  cgnsDoms(zone)%kl + 1

            endif

         enddo

         ! The rest of the block ranges must be obtained by
         ! communication.

         do i=(nBlocks+1),nSubblocks

            ! Receive the range of subblock i.

            source = proc(i)
            call mpi_recv(ii, 6, sumb_integer, source, source, &
                 SUmb_comm_world, status, ierr)

            subRanges(1,1,i) = ii(1)
            subRanges(1,2,i) = ii(2)
            subRanges(2,1,i) = ii(3)
            subRanges(2,2,i) = ii(4)
            subRanges(3,1,i) = ii(5)
            subRanges(3,2,i) = ii(6)
         enddo

         ! Determine the size of the largest subblock.

         do i=1,nSubBlocks
            ll = (subRanges(1,2,i) - subRanges(1,1,i) + 1) &
                 * (subRanges(2,2,i) - subRanges(2,1,i) + 1) &
                 * (subRanges(3,2,i) - subRanges(3,1,i) + 1)
            bufSize = max(bufSize, ll)
         enddo

      else testRoot

         ! Loop over the number of subblocks stored on this processor.

         do i=1,nBlocks

            ! Store the local block id a bit easier in j.

            j = blocksCGNSblock(i+offset)

            ! Copy the range of this subblock into the buffer ii.
            ! This is the same for all spectral solutions, so the
            ! first one can be used.

            ii(1) = flowDoms(j,1,1)%iBegor + 1
            ii(2) = flowDoms(j,1,1)%iEndor
            ii(3) = flowDoms(j,1,1)%jBegor + 1
            ii(4) = flowDoms(j,1,1)%jEndor
            ii(5) = flowDoms(j,1,1)%kBegor + 1
            ii(6) = flowDoms(j,1,1)%kEndor

            ! Correct in case rind layers must be stored.

            if( rindLayerThisSol ) then

               if(ii(1) == 2) ii(1) = 1
               if(ii(2) == cgnsDoms(zone)%il) ii(2) = ii(2) + 1
               if(ii(3) == 2) ii(3) = 1
               if(ii(4) == cgnsDoms(zone)%jl) ii(4) = ii(4) + 1
               if(ii(5) == 2) ii(5) = 1
               if(ii(6) == cgnsDoms(zone)%kl) ii(6) = ii(6) + 1

            endif

            ! Send the buffer to processor 0.

            call mpi_send(ii, 6, sumb_integer, 0, myID, &
                 SUmb_comm_world, ierr)

            ! Check the size of this subblock and update bufSize
            ! if needed.

            ll = (ii(2) - ii(1) + 1) * (ii(4) - ii(3) + 1) &
                 * (ii(6) - ii(5) + 1)
            bufSize = max(bufSize, ll)

         enddo

      endif testRoot

    end subroutine getSubRangesSol
  end subroutine writeSolCGNSZone
end module writeCGNSVolume
