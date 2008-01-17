!
!      ******************************************************************
!      *                                                                *
!      * File:          releaseMemory.f90                               *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 08-16-2004                                      *
!      * Last modified: 11-27-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine releaseMemoryPart1
!
!      ******************************************************************
!      *                                                                *
!      * releaseMemoryPart1 releases all the memory on the coarser      *
!      * grids of flowDoms and the fine grid memory which is not needed *
!      * for the possible interpolation of the spectral solution.       *
!      *                                                                *
!      ******************************************************************
!
       use block
       use inputIteration
       use inputTimeSpectral
       use inputPhysics
       use inputUnsteady
       implicit none
!
!      Local variables
!
       integer :: ierr

       integer(kind=intType) :: sps, nLevels, level, nn
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number of grid levels present in flowDoms.

       nLevels = ubound(flowDoms,2)

       ! Loop over the number of spectral solutions.

       spectralLoop: do sps=1,nTimeIntervalsSpectral

         ! Loop over the coarser grid levels and local blocks and
         ! deallocate all the memory.

         do level=2,nLevels
           do nn=1,nDom
             call deallocateBlock(nn, level, sps)
           enddo
         enddo

         ! Release some memory of the fine grid, which is not needed
         ! anymore.

         do nn=1,nDom
           deallocate(flowDoms(nn,1,sps)%dw,   flowDoms(nn,1,sps)%fw,   &
                      flowDoms(nn,1,sps)%dtl,  flowDoms(nn,1,sps)%radI, &
                      flowDoms(nn,1,sps)%radJ, flowDoms(nn,1,sps)%radK, &
                      stat=ierr)
           if(ierr /= 0)                          &
             call terminate("releaseMemoryPart1", &
                            "Deallocation error for dw, fw, dtl and &
                            &spectral radii.")

           ! Nullify the pointers, such that no attempt is made to
           ! release the memory again.

           nullify(flowDoms(nn,1,sps)%dw)
           nullify(flowDoms(nn,1,sps)%fw)
           nullify(flowDoms(nn,1,sps)%dtl)
           nullify(flowDoms(nn,1,sps)%radI)
           nullify(flowDoms(nn,1,sps)%radJ)
           nullify(flowDoms(nn,1,sps)%radK)

           ! Check if the zeroth stage runge kutta memory has been
           ! allocated. If so deallocate it and nullify the pointers.

           if(smoother == RungeKutta) then

             deallocate(flowDoms(nn,1,sps)%wn, flowDoms(nn,1,sps)%pn, &
                        stat=ierr)
             if(ierr /= 0)                          &
               call terminate("releaseMemoryPart1", &
                              "Deallocation error for wn and pn")

             nullify(flowDoms(nn,1,sps)%wn)
             nullify(flowDoms(nn,1,sps)%pn)

           endif

           ! Release the memory of the old residuals for the time
           ! accurate Runge-Kutta schemes.

           if(equationMode          == unsteady .and. &
              timeIntegrationScheme == explicitRK) then

             deallocate(flowDoms(nn,1,sps)%dwOldRK, stat=ierr)
             if(ierr /= 0)                          &
               call terminate("releaseMemoryPart1", &
                              "Deallocation error for dwOldRK,")

             nullify(flowDoms(nn,1,sps)%dwOldRK)
           endif

         enddo

       enddo spectralLoop

       end subroutine releaseMemoryPart1

!      ==================================================================

       subroutine releaseMemoryPart2
!
!      ******************************************************************
!      *                                                                *
!      * releaseMemoryPart2 releases all the memory of flowDoms on the  *
!      * finest grid as well as the memory allocated in the other       *
!      * modules.                                                       *
!      *                                                                *
!      ******************************************************************
!
       use block
       use inputTimeSpectral
       implicit none
!
!      Local variables
!
       integer :: ierr

       integer(kind=intType) :: nn, sps
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Release the memory of flowDoms of the finest grid and of the
       ! array flowDoms afterwards.

       do sps=1,nTimeIntervalsSpectral
         do nn=1,nDom
           call deallocateBlock(nn, 1_intType, sps)
         enddo
       enddo

       deallocate(flowDoms, stat=ierr)
       if(ierr /= 0)                          &
         call terminate("releaseMemoryPart2", &
                        "Deallocation failure for flowDoms")

       ! Some more memory should be deallocated if this code is to
       ! be used in combination with adaptation.

       end subroutine releaseMemoryPart2

!      ==================================================================

       subroutine deallocateBlock(nn, level, sps)
!
!      ******************************************************************
!      *                                                                *
!      * deallocateBlock deallocates all the allocated memory of the    *
!      * given block.                                                   *
!      *                                                                *
!      ******************************************************************
!
       use block
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nn, level, sps
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: i

       type(viscSubfaceType), dimension(:), pointer :: viscSubface
       type(BCDataType),      dimension(:), pointer :: BCData

       logical :: deallocationFailure
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize deallocationFailure to .false.

       deallocationFailure = .false.

       ! Set the pointer for viscSubface and deallocate the memory
       ! stored in there. Initialize ierr to 0, such that the terminate
       ! routine is only called at the end if a memory deallocation
       ! failure occurs.

       ierr = 0
       viscSubface => flowDoms(nn,level,sps)%viscSubface
       do i=1,flowDoms(nn,level,sps)%nViscBocos
         deallocate(viscSubface(i)%tau,  viscSubface(i)%q, &
                    viscSubface(i)%utau, stat=ierr)
         if(ierr /= 0) deallocationFailure = .true.

         nullify(viscSubface(i)%tau)
         nullify(viscSubface(i)%q)
         nullify(viscSubface(i)%utau)
       enddo

       ! Set the pointer for BCData and deallocate the memory
       ! stored in there.

       BCData => flowDoms(nn,level,sps)%BCData
       do i=1,flowDoms(nn,level,sps)%nBocos

         if( associated(BCData(i)%norm) ) &
            deallocate(BCData(i)%norm, stat=ierr)
         if(ierr /= 0) deallocationFailure = .true.

         if( associated(BCData(i)%rface) ) &
            deallocate(BCData(i)%rface, stat=ierr)
         if(ierr /= 0) deallocationFailure = .true.

         if( associated(BCData(i)%uSlip) ) &
            deallocate(BCData(i)%uSlip, stat=ierr)
         if(ierr /= 0) deallocationFailure = .true.

         if( associated(BCData(i)%TNS_Wall) ) &
            deallocate(BCData(i)%TNS_Wall, stat=ierr)
         if(ierr /= 0) deallocationFailure = .true.

         if( associated(BCData(i)%ptInlet) ) &
            deallocate(BCData(i)%ptInlet, stat=ierr)
         if(ierr /= 0) deallocationFailure = .true.

         if( associated(BCData(i)%ttInlet) ) &
            deallocate(BCData(i)%ttInlet, stat=ierr)
         if(ierr /= 0) deallocationFailure = .true.

         if( associated(BCData(i)%htInlet) ) &
            deallocate(BCData(i)%htInlet, stat=ierr)
         if(ierr /= 0) deallocationFailure = .true.

         if( associated(BCData(i)%flowXdirInlet) ) &
            deallocate(BCData(i)%flowXdirInlet, stat=ierr)
         if(ierr /= 0) deallocationFailure = .true.

         if( associated(BCData(i)%flowYdirInlet) ) &
            deallocate(BCData(i)%flowYdirInlet, stat=ierr)
         if(ierr /= 0) deallocationFailure = .true.

         if( associated(BCData(i)%flowZdirInlet) ) &
            deallocate(BCData(i)%flowZdirInlet, stat=ierr)
         if(ierr /= 0) deallocationFailure = .true.

         if( associated(BCData(i)%rho) ) &
            deallocate(BCData(i)%rho, stat=ierr)
         if(ierr /= 0) deallocationFailure = .true.

         if( associated(BCData(i)%velx) ) &
            deallocate(BCData(i)%velx, stat=ierr)
         if(ierr /= 0) deallocationFailure = .true.

         if( associated(BCData(i)%vely) ) &
            deallocate(BCData(i)%vely, stat=ierr)
         if(ierr /= 0) deallocationFailure = .true.

         if( associated(BCData(i)%velz) ) &
            deallocate(BCData(i)%velz, stat=ierr)
         if(ierr /= 0) deallocationFailure = .true.

         if( associated(BCData(i)%ps) ) &
            deallocate(BCData(i)%ps, stat=ierr)
         if(ierr /= 0) deallocationFailure = .true.

         if( associated(BCData(i)%turbInlet) ) &
           deallocate(BCData(i)%turbInlet, stat=ierr)
         if(ierr /= 0) deallocationFailure = .true.

         nullify(BCData(i)%norm)
         nullify(BCData(i)%rface)
         nullify(BCData(i)%uSlip)
         nullify(BCData(i)%TNS_Wall)
         nullify(BCData(i)%ptInlet)
         nullify(BCData(i)%ttInlet)
         nullify(BCData(i)%htInlet)
         nullify(BCData(i)%flowXdirInlet)
         nullify(BCData(i)%flowYdirInlet)
         nullify(BCData(i)%flowZdirInlet)
         nullify(BCData(i)%rho)
         nullify(BCData(i)%velx)
         nullify(BCData(i)%vely)
         nullify(BCData(i)%velz)
         nullify(BCData(i)%ps)
         nullify(BCData(i)%turbInlet)

       enddo

       ! Deallocate all pointers of this block.

       if( associated(flowDoms(nn,level,sps)%BCType) ) &
         deallocate(flowDoms(nn,level,sps)%BCType, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%BCFaceID) ) &
         deallocate(flowDoms(nn,level,sps)%BCFaceID, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%cgnsSubface) ) &
         deallocate(flowDoms(nn,level,sps)%cgnsSubface, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%inBeg) ) &
         deallocate(flowDoms(nn,level,sps)%inBeg, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%inEnd) ) &
         deallocate(flowDoms(nn,level,sps)%inEnd, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%jnBeg) ) &
         deallocate(flowDoms(nn,level,sps)%jnBeg, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%jnEnd) ) &
         deallocate(flowDoms(nn,level,sps)%jnEnd, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%knBeg) ) &
         deallocate(flowDoms(nn,level,sps)%knBeg, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%knEnd) ) &
         deallocate(flowDoms(nn,level,sps)%knEnd, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%dinBeg) ) &
         deallocate(flowDoms(nn,level,sps)%dinBeg, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%dinEnd) ) &
         deallocate(flowDoms(nn,level,sps)%dinEnd, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%djnBeg) ) &
         deallocate(flowDoms(nn,level,sps)%djnBeg, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%djnEnd) ) &
         deallocate(flowDoms(nn,level,sps)%djnEnd, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%dknBeg) ) &
         deallocate(flowDoms(nn,level,sps)%dknBeg, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%dknEnd) ) &
         deallocate(flowDoms(nn,level,sps)%dknEnd, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%icBeg) ) &
         deallocate(flowDoms(nn,level,sps)%icBeg, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%icEnd) ) &
         deallocate(flowDoms(nn,level,sps)%icEnd, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%jcBeg) ) &
         deallocate(flowDoms(nn,level,sps)%jcBeg, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%jcEnd) ) &
         deallocate(flowDoms(nn,level,sps)%jcEnd, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%kcBeg) ) &
         deallocate(flowDoms(nn,level,sps)%kcBeg, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%kcEnd) ) &
         deallocate(flowDoms(nn,level,sps)%kcEnd, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%neighBlock) ) &
         deallocate(flowDoms(nn,level,sps)%neighBlock, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%neighProc) ) &
         deallocate(flowDoms(nn,level,sps)%neighProc, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%l1) ) &
         deallocate(flowDoms(nn,level,sps)%l1, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%l2) ) &
         deallocate(flowDoms(nn,level,sps)%l2, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%l3) ) &
         deallocate(flowDoms(nn,level,sps)%l3, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%groupNum) ) &
         deallocate(flowDoms(nn,level,sps)%groupNum, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.


       if( associated(flowDoms(nn,level,sps)%iblank) ) &
         deallocate(flowDoms(nn,level,sps)%iblank, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%ibndry) ) &
         deallocate(flowDoms(nn,level,sps)%ibndry, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%idonor) ) &
         deallocate(flowDoms(nn,level,sps)%idonor, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%overint) ) &
         deallocate(flowDoms(nn,level,sps)%overint, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%neighBlockOver) ) &
         deallocate(flowDoms(nn,level,sps)%neighBlockOver, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%neighProcOver) ) &
         deallocate(flowDoms(nn,level,sps)%neighProcOver, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.


       if( associated(flowDoms(nn,level,sps)%BCData) ) &
         deallocate(flowDoms(nn,level,sps)%BCData, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%viscSubface) ) &
         deallocate(flowDoms(nn,level,sps)%viscSubface, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.


       if( associated(flowDoms(nn,level,sps)%viscIminPointer) ) &
         deallocate(flowDoms(nn,level,sps)%viscIminPointer, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%viscImaxPointer) ) &
         deallocate(flowDoms(nn,level,sps)%viscImaxPointer, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%viscJminPointer) ) &
         deallocate(flowDoms(nn,level,sps)%viscJminPointer, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%viscJmaxPointer) ) &
         deallocate(flowDoms(nn,level,sps)%viscJmaxPointer, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%viscKminPointer) ) &
         deallocate(flowDoms(nn,level,sps)%viscKminPointer, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%viscKmaxPointer) ) &
         deallocate(flowDoms(nn,level,sps)%viscKmaxPointer, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.


       if( associated(flowDoms(nn,level,sps)%x) ) &
         deallocate(flowDoms(nn,level,sps)%x, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%xOld) ) &
         deallocate(flowDoms(nn,level,sps)%xOld, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%si) ) &
         deallocate(flowDoms(nn,level,sps)%si, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%sj) ) &
         deallocate(flowDoms(nn,level,sps)%sj, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%sk) ) &
         deallocate(flowDoms(nn,level,sps)%sk, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%vol) ) &
         deallocate(flowDoms(nn,level,sps)%vol, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%volOld) ) &
         deallocate(flowDoms(nn,level,sps)%volOld, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%pori) ) &
         deallocate(flowDoms(nn,level,sps)%pori, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%porj) ) &
         deallocate(flowDoms(nn,level,sps)%porj, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%pork) ) &
         deallocate(flowDoms(nn,level,sps)%pork, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%indFamilyI) ) &
         deallocate(flowDoms(nn,level,sps)%indFamilyI, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%indFamilyJ) ) &
         deallocate(flowDoms(nn,level,sps)%indFamilyJ, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%indFamilyK) ) &
         deallocate(flowDoms(nn,level,sps)%indFamilyK, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%factFamilyI) ) &
         deallocate(flowDoms(nn,level,sps)%factFamilyI, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%factFamilyJ) ) &
         deallocate(flowDoms(nn,level,sps)%factFamilyJ, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%factFamilyK) ) &
         deallocate(flowDoms(nn,level,sps)%factFamilyK, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%rotMatrixI) ) &
         deallocate(flowDoms(nn,level,sps)%rotMatrixI, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%rotMatrixJ) ) &
         deallocate(flowDoms(nn,level,sps)%rotMatrixJ, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%rotMatrixK) ) &
         deallocate(flowDoms(nn,level,sps)%rotMatrixK, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%sFaceI) ) &
         deallocate(flowDoms(nn,level,sps)%sFaceI, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%sFaceJ) ) &
         deallocate(flowDoms(nn,level,sps)%sFaceJ, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%sFaceK) ) &
         deallocate(flowDoms(nn,level,sps)%sFaceK, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.


       if( associated(flowDoms(nn,level,sps)%w) ) &
         deallocate(flowDoms(nn,level,sps)%w, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%wOld) ) &
         deallocate(flowDoms(nn,level,sps)%wOld, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%p) ) &
         deallocate(flowDoms(nn,level,sps)%p, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%gamma) ) &
         deallocate(flowDoms(nn,level,sps)%gamma, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%rlv) ) &
         deallocate(flowDoms(nn,level,sps)%rlv, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%rev) ) &
         deallocate(flowDoms(nn,level,sps)%rev, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%s) ) &
         deallocate(flowDoms(nn,level,sps)%s, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.


       if( associated(flowDoms(nn,level,sps)%p1) ) &
         deallocate(flowDoms(nn,level,sps)%p1, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%dw) ) &
         deallocate(flowDoms(nn,level,sps)%dw, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%fw) ) &
         deallocate(flowDoms(nn,level,sps)%fw, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%dwOldRK) ) &
         deallocate(flowDoms(nn,level,sps)%dwOldRK, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%w1) ) &
         deallocate(flowDoms(nn,level,sps)%w1, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%wr) ) &
         deallocate(flowDoms(nn,level,sps)%wr, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%mgIFine) ) &
         deallocate(flowDoms(nn,level,sps)%mgIFine, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%mgJFine) ) &
         deallocate(flowDoms(nn,level,sps)%mgJFine, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%mgKFine) ) &
         deallocate(flowDoms(nn,level,sps)%mgKFine, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%mgIWeight) ) &
         deallocate(flowDoms(nn,level,sps)%mgIWeight, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%mgJWeight) ) &
         deallocate(flowDoms(nn,level,sps)%mgJWeight, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%mgKWeight) ) &
         deallocate(flowDoms(nn,level,sps)%mgKWeight, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%mgICoarse) ) &
         deallocate(flowDoms(nn,level,sps)%mgICoarse, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%mgJCoarse) ) &
         deallocate(flowDoms(nn,level,sps)%mgJCoarse, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%mgKCoarse) ) &
         deallocate(flowDoms(nn,level,sps)%mgKCoarse, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%iCo) ) &
         deallocate(flowDoms(nn,level,sps)%iCo, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%jCo) ) &
         deallocate(flowDoms(nn,level,sps)%jCo, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%kCo) ) &
         deallocate(flowDoms(nn,level,sps)%kCo, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.


       if( associated(flowDoms(nn,level,sps)%wn) ) &
         deallocate(flowDoms(nn,level,sps)%wn, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%pn) ) &
         deallocate(flowDoms(nn,level,sps)%pn, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%dtl) ) &
         deallocate(flowDoms(nn,level,sps)%dtl, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%radI) ) &
         deallocate(flowDoms(nn,level,sps)%radI, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%radJ) ) &
         deallocate(flowDoms(nn,level,sps)%radJ, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%radK) ) &
         deallocate(flowDoms(nn,level,sps)%radK, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.


       if( associated(flowDoms(nn,level,sps)%d2Wall) ) &
         deallocate(flowDoms(nn,level,sps)%d2Wall, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.


       if( associated(flowDoms(nn,level,sps)%bmti1) ) &
         deallocate(flowDoms(nn,level,sps)%bmti1, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%bmti2) ) &
         deallocate(flowDoms(nn,level,sps)%bmti2, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%bmtj1) ) &
         deallocate(flowDoms(nn,level,sps)%bmtj1, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%bmtj2) ) &
         deallocate(flowDoms(nn,level,sps)%bmtj2, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%bmtk1) ) &
         deallocate(flowDoms(nn,level,sps)%bmtk1, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%bmtk2) ) &
         deallocate(flowDoms(nn,level,sps)%bmtk2, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.


       if( associated(flowDoms(nn,level,sps)%bvti1) ) &
         deallocate(flowDoms(nn,level,sps)%bvti1, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%bvti2) ) &
         deallocate(flowDoms(nn,level,sps)%bvti2, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%bvtj1) ) &
         deallocate(flowDoms(nn,level,sps)%bvtj1, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%bvtj2) ) &
         deallocate(flowDoms(nn,level,sps)%bvtj2, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%bvtk1) ) &
         deallocate(flowDoms(nn,level,sps)%bvtk1, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       if( associated(flowDoms(nn,level,sps)%bvtk2) ) &
         deallocate(flowDoms(nn,level,sps)%bvtk2, stat=ierr)
       if(ierr /= 0) deallocationFailure = .true.

       ! Check for errors in the deallocation.

       if( deallocationFailure ) &
         call terminate("deallocateBlock", &
                        "Something went wrong when deallocating memory")

       ! Nullify the pointers of this block.

       call nullifyFlowDomPointers(nn,level,sps)
 
       end subroutine deallocateBlock
