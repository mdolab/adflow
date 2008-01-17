!
!      ******************************************************************
!      *                                                                *
!      * File:          volSolFileNamesWrite.f90                        *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 10-11-2005                                      *
!      * Last modified: 10-13-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine volSolFileNamesWrite
!
!      ******************************************************************
!      *                                                                *
!      * volSolFileNamesWrite determines the names and number of volume *
!      * solution files to be written. Furthermore it sets the pointers *
!      * and/or allocates the memory for IOVar to make a general        *
!      * treatment of the writing possible.                             *
!      *                                                                *
!      ******************************************************************
!
       use block
       use inputIO
       use inputPhysics
       use inputTimeSpectral
       use IOModule
       use iteration
       use monitor
       use outputMod
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn, mm, kk, nAvail
       integer(kind=intType) :: iEnd, jEnd, kEnd

       character(len=7) :: intString
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
!      ******************************************************************
!      *                                                                *
!      * Determine the names and number of volume solution files to be  *
!      * written.                                                       *
!      *                                                                *
!      ******************************************************************
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
!      ******************************************************************
!      *                                                                *
!      * Set the pointers for IOVar if volume solution files need to be *
!      * written.                                                       *
!      *                                                                *
!      ******************************************************************
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
