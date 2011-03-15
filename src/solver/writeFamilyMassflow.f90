!
!      ******************************************************************
!      *                                                                *
!      * File:          writeFamilyMassflow.f90                         *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-18-2007                                      *
!      * Last modified: 10-30-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine writeFamilyMassflow
!
!      ******************************************************************
!      *                                                                *
!      * writeFamilyMassflow writes the information of the mass flow    *
!      * through families and sliding mesh interfaces to stdout.        *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use cgnsGrid
       use communication
       use inputPhysics
       use inputTimeSpectral
       use monitor
       implicit none
!
!      Local variables.
!
       integer :: nSize, ierr

       integer(kind=intType) :: nn, ii, jj, kk, sps
       integer(kind=intType) :: i1, i2, lenMax, nCharWrite

       real(kind=realType), dimension(2*cgnsNSliding+cgnsNFamilies) :: tmp
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if the mass flow should not be monitored.

       if((.not. monMassSliding) .and. (.not. monMassFamilies)) return

       ! Determine the maximum length of the family names involved.

       lenMax = 0
       if( monMassSliding ) then
         do nn=1,cgnsNSliding

           ! Store the family ID's a bit easier.

           i1 = famIDsSliding(nn,1)
           i2 = famIDsSliding(nn,2)

           ii = len_trim(cgnsFamilies(i1)%familyName)
           lenMax = max(lenMax, ii)

           ii = len_trim(cgnsFamilies(i2)%familyName)
           lenMax = max(lenMax, ii)
         enddo
       endif

       do nn=1,cgnsNFamilies
         if(cgnsFamilies(nn)%monitorMassflow            .and. &
            cgnsFamilies(nn)%BCType /= MassBleedInflow  .and. &
            cgnsFamilies(nn)%BCType /= MassBleedOutflow .and. &
            cgnsFamilies(nn)%BCType /= SlidingInterface) then
           ii = len_trim(cgnsFamilies(nn)%familyName)
           lenMax = max(lenMax, ii)
         endif
       enddo

       lenMax = max(lenMax, 4_intType)

       ! Determine the global values of the mass flow through the
       ! families. These values are only needed on processor 0.

       nSize = ubound(massFlowFamilyInv,1)

       do nn=1,nTimeIntervalsSpectral
         call mpi_reduce(massFlowFamilyInv(1,nn), tmp, nSize, &
                         sumb_real, mpi_sum, 0, SUmb_comm_world, ierr)
         do ii=1,nSize
           massFlowFamilyInv(ii,nn) = tmp(ii)
         enddo

         call mpi_reduce(massFlowFamilyDiss(1,nn), tmp, nSize, &
                         sumb_real, mpi_sum, 0, SUmb_comm_world, ierr)
         do ii=1,nSize
           massFlowFamilyDiss(ii,nn) = tmp(ii)
         enddo
       enddo

       ! Write the data to screen, which is done by processor 0.

       testRootProc: if(myID == 0) then

         ! Determine the number of characters to write on a line.

         nCharWrite = lenMax + 3 + 3*(fieldWidth+1)
         if( equationMode == timeSpectral) &
           nCharWrite = nCharWrite + 11

         ! Check if sliding mesh mass flows must be monitored.

         checkSlidingMesh: if( monMassSliding ) then

           ! Write a nice header.

           print "(a)", "#"
           print "(a)", "# Mass flow through the sliding mesh &
                        &interfaces"

           call headerMassFlow

           ! Loop over the number of number of time instances and the
           ! number of sliding meshes.

           do sps=1,nTimeIntervalsSpectral
             do nn=1,cgnsNSliding

               ! Store the family ID's a bit easier.

               i1 = famIDsSliding(nn,1)
               i2 = famIDsSliding(nn,2)

               ! Determine the first entry in the monitoring arrays.

               ii = 2*abs(cgnsFamilies(i1)%slidingID) - 1

               ! Write the information for the first family.

               write(*,"(a)",advance="no") "#"
               if(equationMode == timeSpectral) &
                 write(*,"(i8,3x)",advance="no") sps

               jj = len_trim(cgnsFamilies(i1)%familyName)
               write(*,"(a)",advance="no") trim(cgnsFamilies(i1)%familyName)
               do kk=jj,lenMax
                 write(*,"(1x)",advance="no")
               enddo

               print 101, massFlowFamilyInv(ii,sps), &
                          massFlowFamilyDiss(ii,sps),&
                          massFlowFamilyInv(ii,sps)  &
                        + massFlowFamilyDiss(ii,sps)

               ! Write the info for the second family.

               ii = ii + 1
               write(*,"(a)",advance="no") "#"
               if(equationMode == timeSpectral) &
                 write(*,"(i8,3x)",advance="no") sps

               jj = len_trim(cgnsFamilies(i2)%familyName)
               write(*,"(a)",advance="no") trim(cgnsFamilies(i2)%familyName)
               do kk=jj,lenMax
                 write(*,"(1x)",advance="no")
               enddo

               print 101, massFlowFamilyInv(ii,sps), &
                          massFlowFamilyDiss(ii,sps),&
                          massFlowFamilyInv(ii,sps)  &
                        + massFlowFamilyDiss(ii,sps)

             enddo
           enddo

           ! Set the value for offset to 2*cgnsNSliding.

           ii = 2*cgnsNSliding

         else checkSlidingMesh

           ! Sliding mesh mass flows are not monitored. Set offset to 0.

           ii = 0

         endif checkSlidingMesh

         ! Check if mass flows for families must be monitored.

         checkFamilies: if( monMassFamilies ) then

           ! Write a nice header.

           print "(a)", "#"
           print "(a)", "# Mass flow through boundary families"

           call headerMassFlow

           ! Loop over the number of number of time instances and the
           ! number of families

           do sps=1,nTimeIntervalsSpectral
             do nn=1,cgnsNFamilies

               ! Check for a family for which the mass flow must be
               ! monitored and which is not a bleed region or a
               ! sliding mesh interface.

               if(cgnsFamilies(nn)%monitorMassflow            .and. &
                  cgnsFamilies(nn)%BCType /= MassBleedInflow  .and. &
                  cgnsFamilies(nn)%BCType /= MassBleedOutflow .and. &
                  cgnsFamilies(nn)%BCType /= SlidingInterface) then

                 ! Update the counter for the entry in the arrays
                 ! massFlowFamilyInv and massFlowFamilyDiss.

                 ii = ii + 1

                 ! Write the information for this family.

                 write(*,"(a)",advance="no") "#"
                 if(equationMode == timeSpectral) &
                   write(*,"(i8,3x)",advance="no") sps

                 jj = len_trim(cgnsFamilies(nn)%familyName)
                 write(*,"(a)",advance="no") trim(cgnsFamilies(nn)%familyName)
                 do kk=jj,lenMax
                   write(*,"(1x)",advance="no")
                 enddo

                 print 101, massFlowFamilyInv(ii,sps), &
                            massFlowFamilyDiss(ii,sps),&
                            massFlowFamilyInv(ii,sps)  &
                          + massFlowFamilyDiss(ii,sps)
               endif
             enddo
           enddo
 
         endif checkFamilies

         ! Print a comment sign for a nice output.

         print "(a)", "#"

 101     format(3(1x,e12.5))

       endif testRootProc

       !=================================================================

       contains

         !===============================================================

         subroutine headerMassFlow
!
!        ****************************************************************
!        *                                                              *
!        * This subroutine writes a nice header for the mass flow       *
!        * output to stdout.                                            *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         write(*,"(a)",advance="no") "#"
         do nn=2,nCharWrite
           write(*,"(a)",advance="no") "-"
         enddo 
         print "(1x)"

         write(*,"(a)",advance="no") "#"
         if(equationMode == timeSpectral) &
           write(*,"(a)",advance="no") " Spectral |"
                      
         write(*,"(a)",advance="no") " Name"
         do nn=5,lenMax
           write(*,"(1x)",advance="no")
         enddo
         write(*,"(a)",advance="no") "|"

         write(*,"(a)",advance="no") "   Central  |"
         write(*,"(a)",advance="no") "Dissipation |"
         write(*,"(a)",advance="no") "    Total   |"
         print "(1x)"

         write(*,"(a)",advance="no") "#"
         do nn=2,nCharWrite
           write(*,"(a)",advance="no") "-"
         enddo
         print "(1x)"

         end subroutine headerMassFlow

       end subroutine writeFamilyMassflow
