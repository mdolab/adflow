!
!      ******************************************************************
!      *                                                                *
!      * File:          allocCoorFineGrid.f90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 06-23-2005                                      *
!      * Last modified: 10-10-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine allocCoorFineGrid
!
!      ******************************************************************
!      *                                                                *
!      * allocCoorFineGrid allocates the memory for all the coordinates *
!      * of all local blocks. Also the memory for the derived data type *
!      * used for the reading is allocated. If an interpolation must be *
!      * performed for the time spectral method the variables of this   *
!      * IO type are allocated as well. For all other cases the pointer *
!      * of the variables are set to the appropriate entry in flowDoms. *
!      *                                                                *
!      ******************************************************************
!
       use block
       use inputPhysics
       use inputTimeSpectral
       use IOModule
       use iteration
       use partitionMod
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn, mm
       integer(kind=intType) :: il, jl, kl, ie, je, ke
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the local blocks and allocate the memory for the
       ! coordinates.

       blockLoop: do nn=1,nDom

         ! Some abbreviations of the block dimensions.

         il = flowDoms(nn,1,1)%il
         jl = flowDoms(nn,1,1)%jl
         kl = flowDoms(nn,1,1)%kl

         ie = flowDoms(nn,1,1)%ie
         je = flowDoms(nn,1,1)%je
         ke = flowDoms(nn,1,1)%ke

         ! Loop over the number of spectral modes and allocate the
         ! memory of the coordinates, single halo's included. The halo
         ! values will be initialized in the preprocessing part.

         do mm=1,nTimeIntervalsSpectral
           allocate(flowDoms(nn,1,mm)%x(0:ie,0:je,0:ke,3), stat=ierr)
           if(ierr /= 0)                         &
             call terminate("allocCoorFineGrid", &
                            "Memory allocation failure for flowDoms%x")
         enddo

         ! For a time accurate computation on deforming meshes, allocate
         ! the memory for the old coordinates. As this is not the time
         ! spectral mode, the third index of flowDoms is always 1.

         if(deforming_Grid .and. equationMode == unsteady) then

           allocate(flowDoms(nn,1,1)%xOld(nOldLevels,0:ie,0:je,0:ke,3), &
                    stat=ierr)
           if(ierr /= 0)                         &
             call terminate("allocCoorFineGrid", &
                            "Memory allocation failure for xOld")
         endif

       enddo blockLoop

       ! Allocate the memory for IOVar.

       allocate(IOVar(nDom,nGridsRead), stat=ierr)
       if(ierr /= 0)                         &
         call terminate("allocCoorFineGrid", &
                        "Memory allocation failure for IOVar")

       ! Determine the equation mode we are solving and set the pointers
       ! of IOVar accordingly, or even allocate the memory, if needed.

       select case(equationMode)

         case (steady)

           ! Steady computation. Only one grid needs to be read.
           ! Loop over the number of blocks and set the pointer.
           ! No pointer offsets are needed for the coordinates.

           do nn=1,nDom
             IOVar(nn,1)%pointerOffset = 0
             IOVar(nn,1)%w => flowDoms(nn,1,1)%x(1:,1:,1:,:)
           enddo

         !===============================================================

         case (unsteady)

           ! Unsteady computation. The first set of coordinates should
           ! be stored in x, other sets (if present) in xOld.
           ! No pointer offsets are needed for the coordinates.

           do nn=1,nDom
             IOVar(nn,1)%pointerOffset = 0
             IOVar(nn,1)%w => flowDoms(nn,1,1)%x(1:,1:,1:,:)

             do mm=2,nGridsRead
               IOVar(nn,mm)%pointerOffset = 0
               IOVar(nn,mm)%w => flowDoms(nn,1,1)%xOld(mm-1,1:,1:,1:,:)
             enddo
           enddo

         !===============================================================

         case (timeSpectral)

           ! Time spectral mode. A further check is required.

           testAllocIOVar: if(interpolSpectral .and. &
                              nGridsRead > 1) then

             ! A restart is performed for a deforming mesh using a
             ! different number of time instances than the previous
             ! computation. Consequently the coordinates, or better
             ! the deformations, will be interpolated later on. Hence
             ! some additional storage is required for the coordinates
             ! to be read and thus the memory for the variables w of
             ! IOVar is allocated. No halo data is needed here.
             ! Note that for the 1st time instance the pointer is set
             ! to the coordinates of flowDoms, because these are not
             ! interpolated. Only the higher time instances are
             ! interpolated. No pointer offsets are needed for the
             ! coordinates.

             do nn=1,nDom
               il = flowDoms(nn,1,1)%il
               jl = flowDoms(nn,1,1)%jl
               kl = flowDoms(nn,1,1)%kl

               IOVar(nn,1)%pointerOffset = 0
               IOVar(nn,1)%w => flowDoms(nn,1,1)%x(1:,1:,1:,:)

               do mm=2,nGridsRead
                 IOVar(nn,mm)%pointerOffset = 0
                 allocate(IOVar(nn,mm)%w(il,jl,kl,3), stat=ierr)
                 if(ierr /= 0)                         &
                   call terminate("allocCoorFineGrid", &
                                  "Memory allocation failure for &
                                  &IOVar%w")
               enddo
             enddo

           else testAllocIOVar

             ! One of the following options is true.
             ! - The computation starts from scratch.
             ! - A restart is performed using a rigid grid, possibly
             !   moving. The number of time instances does not have
             !   to be the same.
             ! - A restart with a deforming mesh is performed with the
             !   same number of time instances.
             !
             ! In all these situations the pointers of IOVar are
             ! simply set to the coordinates of flowDoms.

             do nn=1,nDom
               do mm=1,nGridsRead
                 IOVar(nn,mm)%pointerOffset = 0
                 IOVar(nn,mm)%w => flowDoms(nn,1,mm)%x(1:,1:,1:,:)
               enddo
             enddo

           endif testAllocIOVar

       end select

       end subroutine allocCoorFineGrid
