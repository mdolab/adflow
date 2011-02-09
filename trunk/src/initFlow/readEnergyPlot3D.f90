!
!      ******************************************************************
!      *                                                                *
!      * File:          readEnergyPlot3D.f90                            *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 07-22-2005                                      *
!      * Last modified: 10-18-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readEnergyPlot3D
!
!      ******************************************************************
!      *                                                                *
!      * readEnergyPlot3D reads the total energy from the restart file  *
!      * when Plot3D format is used. It checks if the total energy is   *
!      * present, determines the position and reads it. If not present  *
!      * it is attempted to read the pressure and compute the total     *
!      * energy with the already known values of density, velocity and  *
!      * possibly turbulent kinetic energy.                             *
!      *                                                                *
!      ******************************************************************
!
       use block
       use cgnsNames
       use communication
       use flowVarRefState
       use inputTimeSpectral
       use IOModule
       use restartMod
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn, i, j, k, p, ip, jp, kp
       integer(kind=intType) :: il, jl, kl

       real(kind=realType) :: vx, vy, vz, pres, rhoInv, dummyK
       real(kind=realType), dimension(:,:,:,:), pointer :: w
!
!      Function definitions.
!
       integer(kind=intType) :: bsearchStrings
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Find the index of the total energy in the solution file and
       ! check if it is present.

       p  = nVar
       nn = bsearchStrings(cgnsEnergy, varNames, p)

       testEnergyPresent: if(nn > 0) then

         ! Energy is present. Determine the index in the solution file,
         ! the corresponding offset and read it. Store it at the correct
         ! location in the  w-avriables.

         nn = sorted2Or(nn)
         p  = irhoE
         P3D_Offset = sizeHeader + (nn-1)*sizeVolumeSol

         call readPlot3DVar(solID, pScale, p)

       else testEnergyPresent

         ! Energy is not present. Check for the pressure.

         p  = nVar
         nn = bsearchStrings(cgnsPressure, varNames, p)

         if(nn == 0) then
           if(myID == 0)                        &
             call terminate("readEnergyPlot3D", &
                            "Neither total energy nor pressure present &
                            &in the restart file")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         ! Pressure is present. Determine the index in the solution file,
         ! the corresponding offset and read it. Store at the position
         ! of rhoE.

         nn = sorted2Or(nn)
         p  = irhoE
         P3D_Offset = sizeHeader + (nn-1)*sizeVolumeSol

         call readPlot3DVar(solID, pScale, p)

         ! Loop over the blocks and compute the total energy from the
         ! pressure, density and momentum. This means that it is assumed
         ! that the conservative variables are stored in w and not the
         ! primitive ones. I.e. this routine is only called to construct
         ! states in the past for an unsteady computation.

         domainLoop: do nn=1,nDom

           ! Store the end of the owned cells and the pointer offset
           ! a bit easier. Also abbreviate the pointer for w.

           il = flowDoms(nn,1,1)%il
           jl = flowDoms(nn,1,1)%jl
           kl = flowDoms(nn,1,1)%kl

           p  = IOVar(nn,solID)%pointerOffset
           w => IOVar(nn,solID)%w

           ! Check whether or not a turbulent kinetic energy is present
           ! and compute the total energy for the owned cells.

           testKPresent: if( kPresent ) then

             do k=2,kl
               kp = k+p
               do j=2,jl
                 jp = j+p
                 do i=2,il
                   ip = i+p

                   rhoInv = one/w(ip,jp,kp,irho)
                   vx   = w(ip,jp,kp,imx)*rhoInv
                   vy   = w(ip,jp,kp,imy)*rhoInv
                   vz   = w(ip,jp,kp,imz)*rhoInv
                   pres = w(ip,jp,kp,irhoE)

                   call etotArray(w(ip,jp,kp,irho), vx, vy, vz, pres,  &
                                  w(ip,jp,kp,itu1), w(ip,jp,kp,irhoE), &
                                  kPresent, 1_intType)
                 enddo
               enddo
             enddo

           else testKPresent

             dummyK = zero

             do k=2,kl
               kp = k+p
               do j=2,jl
                 jp = j+p
                 do i=2,il
                   ip = i+p

                   rhoInv = one/w(ip,jp,kp,irho)
                   vx   = w(ip,jp,kp,imx)*rhoInv
                   vy   = w(ip,jp,kp,imy)*rhoInv
                   vz   = w(ip,jp,kp,imz)*rhoInv
                   pres = w(ip,jp,kp,irhoE)

                   call etotArray(w(ip,jp,kp,irho), vx, vy, vz, pres,  &
                                  dummyK,           w(ip,jp,kp,irhoE), &
                                  kPresent, 1_intType)
                 enddo
               enddo
             enddo

           endif testKPresent

         enddo domainLoop

       endif testEnergyPresent

       end subroutine readEnergyPlot3D
