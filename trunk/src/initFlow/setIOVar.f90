!
!      ******************************************************************
!      *                                                                *
!      * File:          setIOVar.f90                                    *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 07-05-2005                                      *
!      * Last modified: 10-10-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setIOVar
!
!      ******************************************************************
!      *                                                                *
!      * setIOVar allocates the memory for the derived data type IOVar, *
!      * which is simplifies the reading. If an interpolation must be   *
!      * performed for the time spectral method also the solution of    *
!      * this IO type is allocated. For all other cases the pointers of *
!      * IOVar are set the the appropriate entries of flowDoms, with    *
!      * possible offset due to the usage of pointers.                  *
!      *                                                                *
!      ******************************************************************
!
       use block
       use flowVarRefState
       use inputPhysics
       use IOModule
       use restartMod
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn, mm, il, jl, kl
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Allocate the memory for IOVar.

       allocate(IOVar(nDom,nSolsRead), stat=ierr)
       if(ierr /= 0)                &
         call terminate("setIOVar", &
                        "Memory allocation failure for solRead")

       ! Determine the equation mode we are solving and set the pointers
       ! of coorRead accordingly, or even allocate the memory, if needed.

       select case(equationMode)

         case (steady)

           ! Steady computation. Only one solution needs to be read.
           ! Loop over the number of blocks and set the pointers.
           ! No pointer offset is needed.

           do nn=1,nDom
             IOVar(nn,1)%pointerOffset = 0
             IOVar(nn,1)%w => flowDoms(nn,1,1)%w(1:,1:,1:,:)
           enddo

         !===============================================================

         case (unsteady)

           ! Unsteady computation. The first solution should be stored in
           ! w. For the others the pointers point to wOld. As the
           ! starting indices of wOld are 2, a pointer shift takes place
           ! here. I know this is a pain in the butt, but that's what
           ! we have to live with.

           do nn=1,nDom
             IOVar(nn,1)%pointerOffset = 0
             IOVar(nn,1)%w => flowDoms(nn,1,1)%w(1:,1:,1:,:)

             do mm=2,nSolsRead
               IOVar(nn,mm)%pointerOffset = -1
               IOVar(nn,mm)%w => flowDoms(nn,1,1)%wOld(mm-1,2:,2:,2:,:)
             enddo
           enddo

         !===============================================================

         case (timeSpectral)

           ! Time spectral mode. A further check is required.

           testAllocSolRead: if( interpolSpectral ) then

             ! A restart is performed using a different number of time
             ! instances than the previous computation. Consequently the
             ! solutions will be interpolated later on. Hence some
             ! additional storage is required for the solutions and thus
             ! the w variables of IOVar are allocated. No halo data is
             ! needed here. No pointer offset either due to the explicit
             ! allocation.

             do nn=1,nDom
               il = flowDoms(nn,1,1)%il
               jl = flowDoms(nn,1,1)%jl
               kl = flowDoms(nn,1,1)%kl

               do mm=1,nSolsRead
                 IOVar(nn,mm)%pointerOffset = 0

                 allocate(IOVar(nn,mm)%w(2:il,2:jl,2:kl,nw), stat=ierr)
                 if(ierr /= 0)                &
                   call terminate("setIOVar", &
                                  "Memory allocation failure for w")
               enddo

             enddo

           else testAllocSolRead

             ! A restart is made using either 1 solution or the correct
             ! number of solution instances. In both cases simply set
             ! the pointers of IOVar. No pointer offset is needed.

             do nn=1,nDom
               do mm=1,nSolsRead
                 IOVar(nn,mm)%pointerOffset = 0
                 IOVar(nn,mm)%w => flowDoms(nn,1,mm)%w(1:,1:,1:,:)
               enddo
             enddo

           endif testAllocSolRead

       end select

       end subroutine setIOVar
