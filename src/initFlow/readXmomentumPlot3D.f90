!
!      ******************************************************************
!      *                                                                *
!      * File:          readXmomentumPlot3D.f90                         *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 07-22-2005                                      *
!      * Last modified: 10-18-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readXmomentumPlot3D
!
!      ******************************************************************
!      *                                                                *
!      * readXmomentumPlot3D reads the momentum in x-direction from the *
!      * restart file when Plot3D format is used. It checks if          *
!      * x-momentum if present and reads it. If not present it checks   *
!      * if the x-velocity is present and constructs the momentum from  *
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
       ! Find the index of the x-momentum in the solution file and check
       ! if it is present.

       p  = nVar
       nn = bsearchStrings(cgnsMomX, varNames, p)

       testXMomPresent: if(nn > 0) then

         ! x-Momentum present. Determine the index in the solution file,
         ! the corresponding offset and read it. Store it at the correct
         ! location.

         nn = sorted2Or(nn)
         p  = imx
         P3D_Offset = sizeHeader + (nn-1)*sizeVolumeSol

         momScale = rhoScale*velScale
         call readPlot3DVar(solID, momScale, p)

       else testXMomPresent

         ! No x-momentum. Try x-velocity. If that is not present,
         ! print an error message and exit.

         p  = nVar
         nn = bsearchStrings(cgnsVelX, varNames, p)

         if(nn == 0) then
           if(myID == 0)                           &
             call terminate("readXmomentumPlot3D", &
                            "Neither x-velocity nor x-momentum present &
                            &in the restart file")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         ! x-velocity is present. Determine the offset in the solution
         ! file and read it. Store it at the position of the x-momentum.

         nn = sorted2Or(nn)
         p  = imx
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
                 w(ip,jp,kp,imx) = w(ip,jp,kp,imx)*w(ip,jp,kp,irho)
               enddo
             enddo
           enddo

         enddo domainLoop

       endif testXMomPresent

       end subroutine readXmomentumPlot3D
