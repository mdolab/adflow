!
!      ******************************************************************
!      *                                                                *
!      * File:          setBCDataFineGrid.f90                           *
!      * Author:        Edwin van der Weide, Seonghyeon Hahn            *
!      * Starting date: 07-07-2004                                      *
!      * Last modified: 11-28-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setBCDataFineGrid(initializationPart)
!
!      ******************************************************************
!      *                                                                *
!      * setBCDataFineGrid extracts the boundary condition data from    *
!      * the cgnsGrid and stores it in useable form in the BCData       *
!      * arrays of the currently finest grid, i.e. groundLevel.         *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use cgnsGrid
       use communication
       use inputTimeSpectral
       use iteration
       use BCDataMod
       implicit none
!
!      Subroutine arguments.
!
       logical, intent(in) :: initializationPart
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: i, j, sps

       logical :: allTurbMassBleedInflow,  allTurbSubsonicInflow
       logical :: allFlowSupersonicInflow, allTurbSupersonicInflow
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize axAssumed and massflowPrescribed to .false.,
       ! indicating that no assumption is made about the axial direction
       ! and no subsonic inflow boundaries with prescribed mass flow
       ! are present.

       axAssumed          = .false.
       massflowPrescribed = .false.

       ! Initialize all the prescribed turbulence as well as flow
       ! variables for supersonic inlet to .true.

       allTurbMassBleedInflow  = .true.
       allTurbSubsonicInflow   = .true.
       allTurbSupersonicInflow = .true.

       allFlowSupersonicInflow = .true.

       ! Loop over the number of spectral solutions and local blocks.

       spectralLoop: do sps=1,nTimeIntervalsSpectral
         domainsLoop: do i=1,nDom

           ! Set the pointers to this block on groundLevel to make
           ! the code readable.

           call setPointers(i,groundLevel,sps)

           ! Loop over the number of boundary condition subfaces.

           bocoLoop: do j=1,nBocos

             ! Store the cgns boundary subface number, the number of
             ! boundary condition data sets and the data sets a bit easier.

             cgnsBoco = cgnsSubface(j)
             nDataSet = &
                   cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%nDataSet
             dataSet => &
                   cgnsDoms(nbkGlobal)%bocoInfo(cgnsBoco)%dataSet

             ! Store the range of the boundary subface a bit easier.

             iBeg = BCData(j)%icBeg; iEnd = BCData(j)%icEnd
             jBeg = BCData(j)%jcBeg; jEnd = BCData(j)%jcEnd

             ! Determine the boundary condition we are having here and
             ! call the appropriate routine.

             select case (BCType(j))

               case (MassBleedInflow)
                 call BCDataMassBleedInflow(j, allTurbMassBleedInflow)

               !=========================================================

               case (NSWallIsothermal)
                 call BCDataIsothermalWall(j)

               !=========================================================

               case (SupersonicInflow)
                 call BCDataSupersonicInflow(j, allFlowSupersonicInflow, &
                                             allTurbSupersonicInflow)

               !=========================================================

               case (SubsonicInflow)
                 call BCDataSubsonicInflow(j, allTurbSubsonicInflow)

               !=========================================================

               case (SubsonicOutflow)
                 call BCDataSubsonicOutflow(j)

               !=========================================================

               case (DomainInterfaceAll, DomainInterfaceRhoUVW, &
                     DomainInterfaceP,   DomainInterfaceRho,    &
                     DomainInterfaceTotal)
                 if( initializationPart ) call BCDataDomainInterface(j)

             end select

           enddo bocoLoop
         enddo domainsLoop
       enddo spectralLoop

       ! If this is the initialization part perform some checks
       ! to see if certain assumptions were made.

       checkInit: if( initializationPart ) then

         ! Check whether or not an assumption was made on the axial
         ! direction. If so, processor 0 prints a warning.

         i = 0
         if( axAssumed ) i = 1
         call mpi_reduce(i, j, 1, sumb_integer, mpi_max, 0, &
                         SUmb_comm_world, ierr)

         if(myID == 0 .and. j == 1) then

           print "(a)", "#"
           print "(a)", "#*==================== !!! Warning !!! &
                        &======================"
           print "(a)", "# Radial boundary data given while no &
                         &rotation axis is present."
           print "(a)", "# It is assumed that the X-axis is the axial &
                        &direction."
           print "(a)", "#*=====================================&
                        &======================"
           print "(a)", "#"

         endif

         ! Check whether or not subsonic inflow boundaries are present with
         ! a prescribed mass flow. If so print a warning that the flow
         ! problem should not be a choked one.

         i = 0
         if( massflowPrescribed ) i = 1
         call mpi_reduce(i, j, 1, sumb_integer, mpi_max, 0, &
                         SUmb_comm_world, ierr)

         if(myID == 0 .and. j == 1) then

           print "(a)", "#"
           print "(a)", "#*==================== !!! Warning !!! &
                        &======================"
           print "(a)", "# Subsonic inflow boundaries present with &
                        &prescribed mass flow."
           print "(a)", "# This is only a well posed problem if the &
                        &flow is not choked."
           print "(a)", "#*=====================================&
                        &======================"
           print "(a)", "#"

         endif

         ! Check whether or not mass bleed inflow regions are present
         ! for which the free stream turbulence is used.

         i = 0
         if(.not. allTurbMassBleedInflow) i = 1
         call mpi_reduce(i, j, 1, sumb_integer, mpi_max, 0, &
                         SUmb_comm_world, ierr)

         if(myID == 0 .and. j == 1) then

           print "(a)", "#"
           print "(a)", "#*==================== !!! Warning !!! &
                        &======================"
           print "(a)", "# Inflow bleed regions present for which the &
                        &turbulence"
           print "(a)", "# quantities are not or insufficiently &
                        &prescribed."
           print "(a)", "# Using free stream values instead."
           print "(a)", "#*=====================================&
                        &======================"
           print "(a)", "#"

         endif

         ! Check whether or not subsonic inflow regions are present
         ! for which the free stream turbulence is used.

         i = 0
         if(.not. allTurbSubsonicInflow) i = 1
         call mpi_reduce(i, j, 1, sumb_integer, mpi_max, 0, &
                         SUmb_comm_world, ierr)

         if(myID == 0 .and. j == 1) then

           print "(a)", "#"
           print "(a)", "#*==================== !!! Warning !!! &
                        &======================"
           print "(a)", "# Subsonic inflow regions present for which &
                        &the turbulence"
           print "(a)", "# quantities are not or insufficiently &
                        &prescribed."
           print "(a)", "# Using free stream values instead."
           print "(a)", "#*=====================================&
                        &======================"
           print "(a)", "#"

         endif

         ! Check whether or not supersonic inflow regions are present
         ! for which the free stream variables is used.

         i = 0
         if(.not. allFlowSupersonicInflow) i = 1
         call mpi_reduce(i, j, 1, sumb_integer, mpi_max, 0, &
                         SUmb_comm_world, ierr)

         if(myID == 0 .and. j == 1) then

           print "(a)", "#"
           print "(a)", "#*==================== !!! Warning !!! &
                        &======================"
           print "(a)", "# Supersonic inflow regions present for which &
                        &the flow variables"
           print "(a)", "# are not or insufficiently prescribed."
           print "(a)", "# Using free stream values instead."
           print "(a)", "#*=====================================&
                        &======================"
           print "(a)", "#"

         endif

         ! Check whether or not supersonic inflow regions are present
         ! for which the free stream turbulence is used.

         i = 0
         if(.not. allTurbSupersonicInflow) i = 1
         call mpi_reduce(i, j, 1, sumb_integer, mpi_max, 0, &
                         SUmb_comm_world, ierr)

         if(myID == 0 .and. j == 1) then

           print "(a)", "#"
           print "(a)", "#*==================== !!! Warning !!! &
                        &======================"
           print "(a)", "# Supersonic inflow regions present for which &
                        &the turbulence"
           print "(a)", "# quantities are not or insufficiently &
                        &prescribed."
           print "(a)", "# Using free stream values instead."
           print "(a)", "#*=====================================&
                        &======================"
           print "(a)", "#"

         endif

       endif checkInit

       end subroutine setBCDataFineGrid

