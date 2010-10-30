!
!      ******************************************************************
!      *                                                                *
!      * File:          readYvelocityPlot3D.f90                         *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 07-21-2005                                      *
!      * Last modified: 10-18-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readYvelocityPlot3D
!
!      ******************************************************************
!      *                                                                *
!      * readYvelocityPlot3D reads the velocity in y-direction from the *
!      * restart file when Plot3D format is used. It checks if          *
!      * y-velocity if present and reads it. If not present it checks   *
!      * if the y-Momentum is present and constructs the velocity from  *
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
       ! Find the index of the y-velocity in the solution file and check
       ! if it is present.

       p  = nVar
       nn = bsearchStrings(cgnsVelY, varNames, p)

       testYVelPresent: if(nn > 0) then

         ! y-Velocity present. Determine the index in the solution file,
         ! the corresponding offset and read it. Store it at the correct
         ! location.

         nn = sorted2Or(nn)
         p  = ivy
         P3D_Offset = sizeHeader + (nn-1)*sizeVolumeSol

         call readPlot3DVar(solID, velScale, p)

       else testYVelPresent

         ! No y-velocity. Try y-momentum. If that is not present,
         ! print an error message and exit.

         p  = nVar
         nn = bsearchStrings(cgnsMomY, varNames, p)

         if(nn == 0) then
           if(myID == 0)                           &
             call terminate("readYvelocityPlot3D", &
                            "Neither y-velocity nor y-momentum present &
                            &in the restart file")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         ! y-Momentum is present. Determine the offset in the solution
         ! file and read it. Store it at the position of the y-velocity.

         nn = sorted2Or(nn)
         p  = ivy
         P3D_Offset = sizeHeader + (nn-1)*sizeVolumeSol

         momScale = rhoScale*velScale
         call readPlot3DVar(solID, momScale, p)

         ! Loop over the owned cells of the blocks and divide the
         ! current value by the density to obtain the velocity.

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
                 w(ip,jp,kp,ivy) = w(ip,jp,kp,ivy)/w(ip,jp,kp,irho)
               enddo
             enddo
           enddo

         enddo domainLoop

       endif testYVelPresent

       end subroutine readYvelocityPlot3D
