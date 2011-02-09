!
!      ******************************************************************
!      *                                                                *
!      * File:          readPressurePlot3D.f90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 07-22-2005                                      *
!      * Last modified: 11-27-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readPressurePlot3D(includeHalos)
!
!      ******************************************************************
!      *                                                                *
!      * readPressurePlot3D reads the pressure from the restart file    *
!      * when Plot3D format is used. It checks if the pressure is       *
!      * present, determines the position and reads it. If not present  *
!      * it is attempted to read the total energy and compute the       *
!      * pressure with the already known values of density, velocity    *
!      * and possibly turbulent kinetic energy.                         *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use cgnsNames
       use communication
       use inputTimeSpectral
       use IOModule
       use restartMod
       implicit none
!
!      Subroutine arguments
!
       logical, intent(in) :: includeHalos
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn, pp
       integer(kind=intType) :: nHalo, iBeg, iEnd, jBeg, jEnd, kBeg, kEnd
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
       ! Set the value of nHalo, depending on whether or not the
       ! halos must be included.

       if( includeHalos ) then
         nHalo = 1
       else
         nHalo = 0
       endif

       ! Find the index of the pressure in the solution file and check
       ! if it is present.

       pp = nVar
       nn = bsearchStrings(cgnsPressure, varNames, pp)

       testPresPresent: if(nn > 0) then

         ! Pressure is present. Determine the index in the solution file,
         ! the corresponding offset and read it. As the general reading
         ! routine is used store it in the irhoE position of w.

         nn = sorted2Or(nn)
         pp = irhoE
         P3D_Offset = sizeHeader + (nn-1)*sizeVolumeSol

         call readPlot3DVar(solID, pScale, pp)

       else testPresPresent

         ! Pressure is not present. Check for the total energy.

         pp = nVar
         nn = bsearchStrings(cgnsEnergy, varNames, pp)

         if(nn == 0) then
           if(myID == 0)                          &
             call terminate("readPressurePlot3D", &
                            "Neither pressure nor total energy present &
                            &in the restart file")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         ! Total energy is present. Determine the index in the solution
         ! file, the corresponding offset and read it.

         nn = sorted2Or(nn)
         pp = irhoE
         P3D_Offset = sizeHeader + (nn-1)*sizeVolumeSol

         call readPlot3DVar(solID, pScale, pp)

         ! Loop over the blocks and compute the pressure
         ! from the total energy, density and velocity.

         domainLoopCompute: do nn=1,nDom

           ! Set the pointers for this block. Make sure that the
           ! correct data is set and array boundaries are not
           ! overwritten.

           pp = min(solID,nTimeIntervalsSpectral)
           call setPointers(nn, 1_intType, pp)

           ! Reset the pointer for w to IOVar.

           w => IOVar(nn,solID)%w

           ! Call the routine computePressure.

           iBeg = 2-nHalo; iEnd = il + nHalo
           jBeg = 2-nHalo; jEnd = jl + nHalo
           kBeg = 2-nHalo; kEnd = kl + nHalo

           call computePressure(iBeg, iEnd, jBeg, jEnd, kBeg, kEnd, &
                                IOVar(nbkLocal,solID)%pointerOffset)

         enddo domainLoopCompute

       endif testPresPresent

       end subroutine readPressurePlot3D
