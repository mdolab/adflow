!
!      ******************************************************************
!      *                                                                *
!      * File:          releaseExtraMemBCs.f90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-23-2004                                      *
!      * Last modified: 09-27-2004                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine releaseExtraMemBCs
!
!      ******************************************************************
!      *                                                                *
!      * releaseExtraMemBCs releases the extra memory allocated in      *
!      * allocMemBcdata. This additional memory was allocated, such     *
!      * that alternative boundary condition treatments can be handled  *
!      * in setBCDataFineGrid.                                          *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use inputTimeSpectral
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: mm, nn, sps, level, nLevels, ii
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

               ! Determine the boundary condition we are having here.

               inflowType: select case (BCType(mm))

                 case (SubsonicInflow)

                   ! Subsonic inflow. Determine the boundary condition
                   ! treatment and release the accordingly. Note that
                   ! the boundary condition treatment of the finest mesh
                   ! must be used, because this data is not available yet
                   ! on the coarse grid.

                   ii = &
                   flowDoms(nn,1,sps)%BCData(mm)%subsonicInletTreatment

                   select case (ii)

                     case (totalConditions)

                       ! Total conditions used. Release the memory for
                       ! the density and velocities and nullify their
                       ! pointers.

                       deallocate(BCData(mm)%rho,  BCData(mm)%velx, &
                                  BCData(mm)%vely, BCData(mm)%velz, &
                                  stat=ierr)
                       if(ierr /= 0) &
                          call terminate("releaseExtraMemBCs", &
                                         "Deallocation failure for rho, &
                                         &velx, vely and velz")

                       nullify(BCData(mm)%rho)
                       nullify(BCData(mm)%velx)
                       nullify(BCData(mm)%vely)
                       nullify(BCData(mm)%velz)

                     !===================================================

                     case (massFlow)

                       ! Full velocity vector and density prescribed at
                       ! inlet boundaries. Release the memory for the
                       ! total conditions and flow directions and nullify
                       ! their pointers.

                       deallocate(BCData(mm)%ptInlet,        &
                                  BCData(mm)%ttInlet,        &
                                  BCData(mm)%htInlet,        &
                                  BCData(mm)%flowXdirInlet, &
                                  BCData(mm)%flowYdirInlet, &
                                  BCData(mm)%flowZdirInlet, stat=ierr)
                       if(ierr /= 0) &
                          call terminate("releaseExtraMemBCs", &
                                         "Deallocation failure for the &
                                         &total conditions.")

                       nullify(BCData(mm)%ptInlet)
                       nullify(BCData(mm)%ttInlet)
                       nullify(BCData(mm)%htInlet)
                       nullify(BCData(mm)%flowXdirInlet)
                       nullify(BCData(mm)%flowYdirInlet)
                       nullify(BCData(mm)%flowZdirInlet)

                   end select

               end select inflowType

             enddo bocoLoop
           enddo domainsLoop
         enddo spectralLoop
       enddo levelLoop

       end subroutine releaseExtraMemBCs
