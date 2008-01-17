!
!      ******************************************************************
!      *                                                                *
!      * File:          readTurbEddyVisPlot3D.f90                       *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 07-22-2005                                      *
!      * Last modified: 10-18-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readTurbEddyVisPlot3D(indW, eddyVisPresent)
!
!      ******************************************************************
!      *                                                                *
!      * readTurbEddyVisPlot3D reads or constructs the eddy viscosity   *
!      * from the restart file when Plot3D format is used. If it is not *
!      * possible to create it eddyVisPresent will be set to .false.    *
!      *                                                                *
!      ******************************************************************
!
       use block
       use cgnsNames
       use communication
       use flowVarRefState
       use inputPhysics
       use IOModule
       use restartMod
       implicit none
!
!      Subroutine arguments.
!
       integer, intent(in)  :: indW
       logical, intent(out) :: eddyVisPresent
!
!      Local variables.
!
       integer(kind=intType) :: nn, i, j, k, p, ip, jp, kp
       integer(kind=intType) :: il, jl, kl

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
       ! Find the index of the eddy-viscosity in the solution file and
       ! check if it is present.

       p  = nVar
       nn = bsearchStrings(cgnsEddy, varNames, p)

       testEddyPresent: if(nn > 0) then

         ! Eddy viscosity is present. Determine the index in the solution
         ! and the corresponding offset sol and read it. Store it in
         ! position indW and set eddyVisPresent to .true.

         nn = sorted2Or(nn)
         p  = indW
         P3D_Offset = sizeHeader + (nn-1)*sizeVolumeSol

         call readPlot3DVar(solID, muScale, p)

         eddyVisPresent = .true.

       else testEddyPresent

         ! Eddy viscosity is not present. Check if the eddy viscosity
         ! ratio is present.

         p  = nVar
         nn = bsearchStrings(cgnsEddyRatio, varNames, p)

         checkEddyRatioPresent: if(nn > 0) then

           ! Eddy viscosity ratio is present. Determine the index in the
           ! solution file, the corresponding offset and read it. Store
           ! it in position indW and set eddyVisPresent to .true.

           nn = sorted2Or(nn)
           p  = indW
           P3D_Offset = sizeHeader + (nn-1)*sizeVolumeSol

           call readPlot3DVar(solID, one, p)

           eddyVisPresent = .true.

           ! Multiply the eddy viscosity by the laminar viscosity such
           ! that it contains the correct nonDimensional value.
           ! As the laminar viscosity is not yet know, use the free
           ! stream viscosity.

           domainLoop: do nn=1,nDom

             ! Store the end of the owned cells and the pointer offset
             ! a bit easier. Also abbreviate the pointer for w.

             il = flowDoms(nn,1,1)%il
             jl = flowDoms(nn,1,1)%jl
             kl = flowDoms(nn,1,1)%kl

             p  = IOVar(nn,solID)%pointerOffset
             w => IOVar(nn,solID)%w

             ! Loop over the owned cells and compute the eddy viscosity.

             do k=2,kl
               kp = k+p
               do j=2,jl
                 jp = j+p
                 do i=2,il
                   ip = i+p

                   w(ip,jp,kp,indW) = w(ip,jp,kp,indW)*muInf
                 enddo
               enddo
             enddo

           enddo domainLoop

         else checkEddyRatioPresent

           ! Eddy viscosity could not be constructed.
           ! Set eddyVisPresent to .false. and initialize it to
           ! the free stream value.

           eddyVisPresent = .false.

           domainLoopSet: do nn=1,nDom

             ! Store the end of the owned cells and the pointer offset
             ! a bit easier. Also abbreviate the pointer for w.

             il = flowDoms(nn,1,1)%il
             jl = flowDoms(nn,1,1)%jl
             kl = flowDoms(nn,1,1)%kl

             p  = IOVar(nn,solID)%pointerOffset
             w => IOVar(nn,solID)%w

             ! Loop over the owned cells and compute the eddy viscosity.

             do k=2,kl
               kp = k+p
               do j=2,jl
                 jp = j+p
                 do i=2,il
                   ip = i+p

                   w(ip,jp,kp,indW) = muInf*eddyVisInfRatio
                 enddo
               enddo
             enddo

           enddo domainLoopSet

         endif checkEddyRatioPresent

       endif testEddyPresent

      end subroutine readTurbEddyVisPlot3D
