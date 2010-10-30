!
!      ******************************************************************
!      *                                                                *
!      * File:          readZmomentumPlot3D.f90                         *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 07-22-2005                                      *
!      * Last modified: 10-18-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readZmomentumPlot3D
!
!      ******************************************************************
!      *                                                                *
!      * readZmomentumPlot3D reads the momentum in z-direction from the *
!      * restart file when Plot3D format is used. It checks if          *
!      * z-momentum if present and reads it. If not present it checks   *
!      * if the z-velocity is present and constructs the momentum from  *
!      * that. The density is already known.                            *
!      *                                                                *
!      ******************************************************************
!
       use block
       use cgnsNames
       use communication
       use IOModule
       use restartMod
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn, i, j, k, p, ip, jp, kp
       integer(kind=intType) :: il, jl, kl

       real(kind=realType) :: momScale
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
       ! Find the index of the z-momentum in the solution file and check
       ! if it is present.

       p  = nVar
       nn = bsearchStrings(cgnsMomZ, varNames, p)

       testZMomPresent: if(nn > 0) then

         ! z-Momentum present. Determine the index in the solution file,
         ! the corresponding offset and read it. Store it at the correct
         ! location.

         nn = sorted2Or(nn)
         p  = imz
         P3D_Offset = sizeHeader + (nn-1)*sizeVolumeSol

         momScale = rhoScale*velScale
         call readPlot3DVar(solID, momScale, p)

       else testZMomPresent

         ! No z-momentum. Try z-velocity. If that is not present,
         ! print an error message and exit.

         p  = nVar
         nn = bsearchStrings(cgnsVelZ, varNames, p)

         if(nn == 0) then
           if(myID == 0)                           &
             call terminate("readZmomentumPlot3D", &
                            "Neither z-velocity nor z-momentum present &
                            &in the restart file")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         ! z-velocity is present. Determine the offset in the solution
         ! file and read it. Store it at the position of the z-momentum.

         nn = sorted2Or(nn)
         p  = imz
         P3D_Offset = sizeHeader + (nn-1)*sizeVolumeSol

         call readPlot3DVar(solID, velScale, p)

         ! Loop over the owned cells of the blocks and multiply the
         ! current value by the density to obtain the momentum.

         domainLoop: do nn=1,nDom

           ! Store the end of the owned cells and the pointer offset
           ! a bit easier. Also abbreviate the pointer for w.

           il = flowDoms(nn,1,1)%il
           jl = flowDoms(nn,1,1)%jl
           kl = flowDoms(nn,1,1)%kl

           p  = IOVar(nn,solID)%pointerOffset
           w => IOVar(nn,solID)%w

           do k=2,kl
             kp = k+p
             do j=2,jl
               jp = j+p
               do i=2,il
                 ip = i+p
                 w(ip,jp,kp,imz) = w(ip,jp,kp,imz)*w(ip,jp,kp,irho)
               enddo
             enddo
           enddo

         enddo domainLoop

       endif testZMomPresent

       end subroutine readZmomentumPlot3D
