!
!      ******************************************************************
!      *                                                                *
!      * File:          tempMemory.f90                                  *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-10-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine deallocateTempMemory(resNeeded)
!
!      ******************************************************************
!      *                                                                *
!      * deallocateTempMemory deallocates memory used in the solver,    *
!      * but which is not needed to store the actual solution. In this  *
!      * way the memory can be used differently, e.g. when writing the  *
!      * solution or computing the wall distances.                      *
!      *                                                                *
!      ******************************************************************
!
       use block
       use communication
       use inputIteration
       use inputTimeSpectral
       implicit none
!
!      Subroutine arguments.
!
       logical, intent(in) :: resNeeded
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn, mm
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Deallocate the communication buffers

       deallocate(sendBuffer, recvBuffer, stat=ierr)
       if(ierr /= 0)                              &
         call terminate("deallocateTempMemory", &
                        "Deallocation error for communication buffers")

       ! Loop over the spectral modes and domains. Note that only memory
       ! on the finest grid is released, because a) most of these
       ! variables are only allocated on the fine grid and b) the coarser
       ! grids do not contribute that much in the memory usage anyway.

       spectralModes: do mm=1,nTimeIntervalsSpectral
         domains: do nn=1,nDom

           ! Check if the residual, time step, etc. Is needed.

           if(.not. resNeeded) then

             ! Residual, etc. Not needed.
             ! Deallocate residual, the time step and the spectral radii
             ! of the fine level.

             deallocate(flowDoms(nn,1,mm)%dw,   flowDoms(nn,1,mm)%fw,   &
                        flowDoms(nn,1,mm)%dtl,  flowDoms(nn,1,mm)%radI, &
                        flowDoms(nn,1,mm)%radJ, flowDoms(nn,1,mm)%radK, &
                        stat=ierr)
             if(ierr /= 0)                            &
               call terminate("deallocateTempMemory", &
                              "Deallocation error for dw, fw, dtl and &
                              &spectral radii.")
           endif

           ! The memory for the zeroth Runge Kutta stage
           ! if a Runge Kutta scheme is used.

           if(smoother == RungeKutta) then

             deallocate(flowDoms(nn,1,mm)%wn, flowDoms(nn,1,mm)%pn, &
                        stat=ierr)
             if(ierr /= 0)                            &
               call terminate("deallocateTempMemory", &
                              "Deallocation error for wn and pn")
           endif

         enddo domains
       enddo spectralModes

       end subroutine deallocateTempMemory

!      ==================================================================

       subroutine allocateTempMemory(resNeeded)
!
!      ******************************************************************
!      *                                                                *
!      * AllocateTempMemory allocates the memory again that was         *
!      * temporarily deallocted by deallocateTempMemory.                *
!      *                                                                *
!      ******************************************************************
!
       use block
       use communication
       use constants
       use flowVarRefState
       use inputIteration
       use inputTimeSpectral
       use iteration
       implicit none
!
!      Subroutine arguments.
!
       logical, intent(in) :: resNeeded
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn,mm
       integer(kind=intType) :: il, jl, kl, ie, je, ke, ib, jb, kb
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! The memory for the receive buffers.

       allocate(sendBuffer(sendBufferSize), &
                recvBuffer(recvBufferSize), stat=ierr)
       if(ierr /= 0)                          &
         call terminate("allocateTempMemory", &
                        "Memory allocation failure for comm buffers")

       ! Loop over the spectral modes and domains. Note that only memory
       ! on the finest mesh level needs to be reallocated, because the
       ! memory on the coarser levels has not been released or is not
       ! needed .

       spectralModes: do mm=1,nTimeIntervalsSpectral
         domains: do nn=1,nDom

           ! Store some dimensions a bit easier.

           il = flowDoms(nn,1,mm)%il
           jl = flowDoms(nn,1,mm)%jl
           kl = flowDoms(nn,1,mm)%kl

           ie = flowDoms(nn,1,mm)%ie
           je = flowDoms(nn,1,mm)%je
           ke = flowDoms(nn,1,mm)%ke

           ib = flowDoms(nn,1,mm)%ib
           jb = flowDoms(nn,1,mm)%jb
           kb = flowDoms(nn,1,mm)%kb

           ! Check if the residual, time step, etc. was deallocated.

           if(.not. resNeeded) then

             ! Allocate the residual, the time step and
             ! the spectral radii.

             allocate(flowDoms(nn,1,mm)%dw(0:ib,0:jb,0:kb,1:nw),  &
                      flowDoms(nn,1,mm)%fw(0:ib,0:jb,0:kb,1:nwf), &
                      flowDoms(nn,1,mm)%dtl(1:ie,1:je,1:ke),      &
                      flowDoms(nn,1,mm)%radI(1:ie,1:je,1:ke),     &
                      flowDoms(nn,1,mm)%radJ(1:ie,1:je,1:ke),     &
                      flowDoms(nn,1,mm)%radK(1:ie,1:je,1:ke), stat=ierr)
             if(ierr /= 0)                            &
               call terminate("allocateTempMemory", &
                              "Memory allocation failure for dw, fw, &
                              &dtl and the spectral radii.")

             ! Initialize dw and fw to zero to avoid possible overflows
             ! of the halo's.

             flowDoms(nn,1,mm)%dw = zero
             flowDoms(nn,1,mm)%fw = zero

           endif

           ! The memory for the zeroth runge kutta stage
           ! if a runge kutta scheme is used.

           if(smoother == RungeKutta) then

             allocate(flowDoms(nn,1,mm)%wn(2:il,2:jl,2:kl,1:nMGVar), &
                      flowDoms(nn,1,mm)%pn(2:il,2:jl,2:kl), stat=ierr)
             if(ierr /= 0)                            &
               call terminate("allocateTempMemory", &
                              "Memory allocation failure for wn and pn")
           endif

         enddo domains
       enddo spectralModes

       end subroutine allocateTempMemory
