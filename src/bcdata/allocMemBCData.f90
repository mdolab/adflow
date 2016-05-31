!
!      ******************************************************************
!      *                                                                *
!      * File:          allocMemBCData.f90                              *
!      * Author:        Edwin van der Weide, Seonghyeon Hahn            *
!      * Starting date: 09-13-2004                                      *
!      * Last modified: 09-13-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine allocMemBCData
!
!      ******************************************************************
!      *                                                                *
!      * allocMemBCData allocates the memory for the prescribed         *
!      * boundary data for all multigrid levels and all spectral        *
!      * solutions for all blocks.                                      *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use flowVarRefState
       use inputTimeSpectral
       use iteration
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: mm, nn, sps, level, nLevels
       integer(kind=intType) :: iBeg, jBeg, iEnd, jEnd
       integer(kind=intType) :: inodeBeg, jnodeBeg, inodeEnd, jnodeEnd
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number of multigrid levels.

       nLevels = ubound(flowDoms,2)

       ! Loop over the number of multigrid level, spectral solutions
       ! and local blocks.

       levelLoop: do level=1,nLevels
         spectralLoop: do sps=1,nTimeIntervalsSpectral
           domainsLoop: do nn=1,nDom

             ! Have the pointers in blockPointers point to the
             ! current block to make everything more readable.

             call setPointers(nn, level, sps)

             ! Loop over the number of boundary subfaces for this block.

             bocoLoop: do mm=1,nBocos

               ! Store the cell range of the boundary subface
               ! a bit easier.

               iBeg = BCData(mm)%icbeg; iEnd = BCData(mm)%icend
               jBeg = BCData(mm)%jcbeg; jEnd = BCData(mm)%jcend

               inodeBeg = BCData(mm)%inbeg; inodeEnd = BCData(mm)%inend
               jnodeBeg = BCData(mm)%jnbeg; jnodeEnd = BCData(mm)%jnend

               ! Initialize all mask values to 1
               bcData(mm)%mask = 1

               ! Determine the boundary condition we are having here
               ! and allocate the memory accordingly.

               select case (BCType(mm))

                 case (NSWallAdiabatic)

                   allocate(BCData(mm)%uSlip(iBeg:iEnd,jBeg:jEnd,3), &
                            BCData(mm)%uSlipALE(0:nALEsteps,iBeg:iEnd,jBeg:jEnd,3), &
                            BCData(mm)%TNS_Wall(iBeg:iEnd,jBeg:jEnd), &
                            BCData(mm)%F(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd,3), &
                            BCData(mm)%dualArea(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd), &
                            BCData(mm)%sHeatFlux(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd), &
                            BCData(mm)%Fp(iBeg:iEnd,jBeg:jEnd,3), &
                            BCData(mm)%Fv(iBeg:iEnd,jBeg:jEnd,3), &
                            ! Note: iBlank/delta are cell based, but uses the
                            ! node number to guarantee a halo exists. 
                            BCData(mm)%iBlank(iNodeBeg:iNodeEnd+1, jNodeBeg:jnodeEnd+1), &
                            BCData(mm)%delta(iNodeBeg:iNodeEnd+1, jNodeBeg:jnodeEnd+1), &
                            BCData(mm)%deltaNode(iNodeBeg:iNodeEnd, jNodeBeg:jNodeEnd), &
                            stat=ierr)
                   if(ierr /= 0)                      &
                     call returnFail("allocMemBCData", &
                                    "Memory allocation failure for &
                                    &an adiabatic wall")
                   BCData(mm)%iBlank = 1

                 !=======================================================

                 case (NSWallIsothermal)

                   allocate(BCData(mm)%uSlip(iBeg:iEnd,jBeg:jEnd,3),  &
                            BCData(mm)%uSlipALE(0:nALEsteps,iBeg:iEnd,jBeg:jEnd,3), &
                            BCData(mm)%TNS_Wall(iBeg:iEnd,jBeg:jEnd), &
                            BCData(mm)%F(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd,3), &
                            BCData(mm)%dualArea(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd), &
                            BCData(mm)%sHeatFlux(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd), &
                            BCData(mm)%Fp(iBeg:iEnd,jBeg:jEnd,3), &
                            BCData(mm)%Fv(iBeg:iEnd,jBeg:jEnd,3), &
                            ! Note: iBlank/delta are cell based, but uses the
                            ! node number to guarantee a halo exists. 
                            BCData(mm)%iBlank(iNodeBeg:iNodeEnd+1, jNodeBeg:jnodeEnd+1), &
                            BCData(mm)%delta(iNodeBeg:iNodeEnd+1, jNodeBeg:jnodeEnd+1), &
                            BCData(mm)%deltaNode(iNodeBeg:iNodeEnd, jNodeBeg:jNodeEnd), &
                            stat=ierr)
                   if(ierr /= 0)                      &
                     call returnFail("allocMemBCData", &
                                    "Memory allocation failure for &
                                    &an isothermal wall")
                   BCData(mm)%iBlank = 1

                 !=======================================================

                 case (EulerWall)

                   allocate(BCData(mm)%rface(iBeg:iEnd,jBeg:jEnd), &
                            BCData(mm)%rFaceALE(0:nALEsteps,iBeg:iEnd,jBeg:jEnd), &
                            BCData(mm)%TNS_Wall(iBeg:iEnd,jBeg:jEnd), &
                            BCData(mm)%F(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd,3), &
                            BCData(mm)%dualArea(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd), &
                            BCData(mm)%sHeatFlux(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd), &
                            BCData(mm)%Fp(iBeg:iEnd,jBeg:jEnd,3), &
                            BCData(mm)%Fv(iBeg:iEnd,jBeg:jEnd,3), &
                            ! Note: iBlank/delta are cell based, but uses the
                            ! node number to guarantee a halo exists. 
                            BCData(mm)%iBlank(iNodeBeg:iNodeEnd+1, jNodeBeg:jnodeEnd+1), &
                            BCData(mm)%delta(iNodeBeg:iNodeEnd+1, jNodeBeg:jnodeEnd+1), &
                            BCData(mm)%deltaNode(iNodeBeg:iNodeEnd, jNodeBeg:jNodeEnd), &
                            stat=ierr)
                   if(ierr /= 0)                      &
                     call returnFail("allocMemBCData", &
                                    "Memory allocation failure for &
                                    &an Euler wall")
                   BCData(mm)%iBlank = 1
                   !=======================================================

                 case (farField)

                   ! Just allocate the memory for the normal mesh
                   ! velocity.

                   allocate(BCData(mm)%rface(iBeg:iEnd,jBeg:jEnd), &
                            BCData(mm)%rFaceALE(0:nALEsteps,iBeg:iEnd,jBeg:jEnd), &
                            BCData(mm)%F(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd,3), &
                            BCData(mm)%sHeatFlux(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd), &
                            BCData(mm)%dualArea(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd), &

                        stat=ierr)
                   if(ierr /= 0)                      &
                     call returnFail("allocMemBCData", &
                                    "Memory allocation failure for &
                                    &a farfield")
                   !=======================================================

                 case (symm, symmPolar)

                   ! Allocate for symm as well. This is not necessary
                   ! but we need it for the reverse AD.

     ! *******************************
     ! Modified by HDN
     ! *******************************
                   allocate(BCData(mm)%rface(iBeg:iEnd,jBeg:jEnd), &
                            BCData(mm)%rFaceALE(0:nALEsteps,iBeg:iEnd,jBeg:jEnd), &
                            BCData(mm)%F(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd,3), &
                            BCData(mm)%sHeatFlux(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd), &
                            BCData(mm)%dualArea(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd), &
                            stat=ierr)
                   if(ierr /= 0)                      &
                     call returnFail("allocMemBCData", &
                                    "Memory allocation failure for &
                                    &a symm")

                 !=======================================================

                 case (SupersonicInflow, DomainInterfaceAll)

                   ! Supersonic inflow or a domain interface with
                   ! all the data prescribed. Allocate the memory for
                   ! the entire state vector to be prescribed.

                   allocate(BCData(mm)%rho(iBeg:iEnd,jBeg:jEnd),   &
                            BCData(mm)%velx(iBeg:iEnd,jBeg:jEnd),  &
                            BCData(mm)%vely(iBeg:iEnd,jBeg:jEnd),  &
                            BCData(mm)%velz(iBeg:iEnd,jBeg:jEnd),  &
                            BCData(mm)%ps(iBeg:iEnd,jBeg:jEnd),    &
                            BCData(mm)%F(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd,3), &
                            BCData(mm)%sHeatFlux(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd), &
                            stat=ierr)
                   if(ierr /= 0)                      &
                     call returnFail("allocMemBCData", &
                                    "Memory allocation failure for &
                                    &a supersonic inflow")

                   ! Check if memory for the turbulent variables must
                   ! be allocated. If so, do so.

                   if(nt2 >= nt1) then
                     allocate(&
                      BCData(mm)%turbInlet(iBeg:iEnd,jBeg:jEnd,nt1:nt2), &
                      stat=ierr)
                     if(ierr /= 0)                      &
                       call returnFail("allocMemBCData", &
                                      "Memory allocation failure for &
                                      &turbInlet for a supersonic &
                                      &inflow")
                   endif

                 !=======================================================

     ! *******************************
     ! Added by HDN
     ! *******************************
                 case (SupersonicOutflow)
                   ! No state is needed for this boco
                   ! To avoid seg fault during force and heat flux extraction,
                   ! F and sHeatFlux are allocated

                   allocate(BCData(mm)%F(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd,3), &
                            BCData(mm)%sHeatFlux(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd), &
                            stat=ierr)
                   if(ierr /= 0)                      &
                     call terminate("allocMemBCData", &
                                    "Memory allocation failure for &
                                    &a supersonic outflow")

                 !=======================================================

                 case (SubsonicInflow)

                   ! Subsonic inflow. Allocate the memory for the
                   ! variables needed. Note the there are two ways to
                   ! specify boundary conditions for a subsonic inflow.

                   allocate(BCData(mm)%flowXdirInlet(iBeg:iEnd,jBeg:jEnd), &
                            BCData(mm)%flowYdirInlet(iBeg:iEnd,jBeg:jEnd), &
                            BCData(mm)%flowZdirInlet(iBeg:iEnd,jBeg:jEnd), &
                            BCData(mm)%ptInlet(iBeg:iEnd,jBeg:jEnd),       &
                            BCData(mm)%ttInlet(iBeg:iEnd,jBeg:jEnd),       &
                            BCData(mm)%htInlet(iBeg:iEnd,jBeg:jEnd),       &
                            BCData(mm)%rho(iBeg:iEnd,jBeg:jEnd),           &
                            BCData(mm)%velx(iBeg:iEnd,jBeg:jEnd),          &
                            BCData(mm)%vely(iBeg:iEnd,jBeg:jEnd),          &
                            BCData(mm)%velz(iBeg:iEnd,jBeg:jEnd),          &
                            BCData(mm)%F(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd,3), &
                            BCData(mm)%sHeatFlux(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd), &
                            stat=ierr)
                   if(ierr /= 0)                      &
                     call returnFail("allocMemBCData", &
                                    "Memory allocation failure for &
                                    &a subsonic inflow")

                   ! Check if memory for the turbulent variables must
                   ! be allocated. If so, do so.

                   if(nt2 >= nt1) then
                     allocate(&
                      BCData(mm)%turbInlet(iBeg:iEnd,jBeg:jEnd,nt1:nt2), &
                      stat=ierr)
                     if(ierr /= 0)                      &
                       call returnFail("allocMemBCData", &
                                      "Memory allocation failure for &
                                      &turbInlet for a subsonic inflow")
                   endif

                 !=======================================================

                 case (SubsonicOutflow, MassBleedOutflow, &
                       DomainInterfaceP)

                   ! Subsonic outflow, outflow mass bleed or domain
                   ! interface with prescribed pressure. Allocate the
                   ! memory for the static pressure.

                   allocate(BCData(mm)%ps(iBeg:iEnd,jBeg:jEnd), &
                            BCData(mm)%F(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd,3), &
                            BCData(mm)%sHeatFlux(iNodeBeg:iNodeEnd,jNodeBeg:jNodeEnd), &
                            stat=ierr)
                   if(ierr /= 0)                      &
                     call returnFail("allocMemBCData", &
                                    "Memory allocation failure for &
                                    &a subsonic outflow, outflow mass &
                                    &bleed or domain interface with &
                                    &prescribed pressure.")

                   ! Initialize the pressure to avoid problems for
                   ! the bleed flows.

                   BCData(mm)%ps = zero

                 !=======================================================

                 case (DomainInterfaceRhoUVW)

                   ! Domain interface with prescribed density and 
                   ! velocities, i.e. mass flow is prescribed. Allocate
                   ! the memory for the variables needed.

                   allocate(BCData(mm)%rho(iBeg:iEnd,jBeg:jEnd),  &
                            BCData(mm)%velx(iBeg:iEnd,jBeg:jEnd), &
                            BCData(mm)%vely(iBeg:iEnd,jBeg:jEnd), &
                            BCData(mm)%velz(iBeg:iEnd,jBeg:jEnd), &
                            stat=ierr)
                   if(ierr /= 0)                      &
                     call returnFail("allocMemBCData", &
                                    "Memory allocation failure for &
                                    &a domain interface with a &
                                    &prescribed mass flow")

                   ! Check if memory for the turbulent variables must
                   ! be allocated. If so, do so.

                   if(nt2 >= nt1) then
                     allocate(&
                      BCData(mm)%turbInlet(iBeg:iEnd,jBeg:jEnd,nt1:nt2), &
                      stat=ierr)
                     if(ierr /= 0)                      &
                       call returnFail("allocMemBCData", &
                                      "Memory allocation failure for &
                                      &turbInlet for a domain interface &
                                      &with a prescribed mass flow")
                   endif

                 !=======================================================

                 case (DomainInterfaceTotal)

                   ! Domain interface with prescribed total conditions.
                   ! Allocate the memory for the variables needed.

                   allocate(BCData(mm)%flowXdirInlet(iBeg:iEnd,jBeg:jEnd), &
                            BCData(mm)%flowYdirInlet(iBeg:iEnd,jBeg:jEnd), &
                            BCData(mm)%flowZdirInlet(iBeg:iEnd,jBeg:jEnd), &
                            BCData(mm)%ptInlet(iBeg:iEnd,jBeg:jEnd),       &
                            BCData(mm)%ttInlet(iBeg:iEnd,jBeg:jEnd),       &
                            BCData(mm)%htInlet(iBeg:iEnd,jBeg:jEnd),       &
                            stat=ierr)
                   if(ierr /= 0)                      &
                     call returnFail("allocMemBCData", &
                                    "Memory allocation failure for &
                                    &a domain interface with total &
                                    &conditions")

                   ! Check if memory for the turbulent variables must
                   ! be allocated. If so, do so.

                   if(nt2 >= nt1) then
                     allocate(&
                      BCData(mm)%turbInlet(iBeg:iEnd,jBeg:jEnd,nt1:nt2), &
                      stat=ierr)
                     if(ierr /= 0)                      &
                       call returnFail("allocMemBCData", &
                                      "Memory allocation failure for &
                                      &turbInlet for a domain interface &
                                      &with a prescribed mass flow")
                   endif

                 !=======================================================

                 case (domainInterfaceRho)

                   ! Domain interface with prescribed density. 
                   ! Allocate the memory for the density.

                   allocate(BCData(mm)%rho(iBeg:iEnd,jBeg:jEnd), &
                            stat=ierr)
                   if(ierr /= 0)                      &
                     call returnFail("allocMemBCData", &
                                    "Memory allocation failure for &
                                    &a domain interface")

               end select

             enddo bocoLoop

           enddo domainsLoop
         enddo spectralLoop
       enddo levelLoop

       end subroutine allocMemBCData
