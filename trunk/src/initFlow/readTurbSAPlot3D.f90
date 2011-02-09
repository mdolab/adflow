!
!      ******************************************************************
!      *                                                                *
!      * File:          readTurbSAPlot3D.f90                            *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 07-22-2005                                      *
!      * Last modified: 10-18-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readTurbSAPlot3D
!
!      ******************************************************************
!      *                                                                *
!      * readTurbSAPlot3D reads or constructs the nu tilde transport    *
!      * variable for the Spalart-Allmaras type of turbulence model     *
!      * from the restart file when Plot3D format is used. If it is not *
!      * possible to create it nu tilde will be set to free stream and  *
!      * a warning will be printed.                                     *
!      *                                                                *
!      ******************************************************************
!
       use block
       use cgnsNames
       use communication
       use flowVarRefState
       use IOModule
       use restartMod
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: nn, i, j, k, p, ip, jp, kp
       integer(kind=intType) :: il, jl, kl

       real(kind=realType) :: ratio, nu, nuScale

       real(kind=realType), dimension(:,:,:,:), pointer :: w

       logical :: eddyVisPresent
!
!      Function definitions.
!
       integer(kind=intType) :: bsearchStrings
       real(kind=realType)   :: saNuKnownEddyRatio
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Compute the scale for the kinematic viscosity.

       nuScale = muScale/rhoScale

       ! Find the index of the nu tilde in the solution file and check
       ! if it is present.

       p  = nVar
       nn = bsearchStrings(cgnsTurbSANu, varNames, p)

       testNuPresent: if(nn > 0) then

         ! Transport variable is present. Determine the index in the
         ! solution file, the offset from the beginning of the file
         ! where this solution record starts and call the general
         ! reading routine.

         nn = sorted2Or(nn)
         p  = itu1
         P3D_Offset = sizeHeader + (nn-1)*sizeVolumeSol

         call readPlot3DVar(solID, nuScale, p)

       else testNuPresent

         ! Nu tilde not present. Try to construct the eddy viscosity.

         call readTurbEddyVisPlot3D(itu1, eddyVisPresent)

         ! Check if the eddy viscosity could be constructed.

         checkEddyPresent: if( eddyVisPresent ) then

           ! Eddy viscosity is present. Compute nu tilde from this
           ! value. Loop over the number of domains.

           domainLoopCompute: do nn=1,nDom

             ! Store the end of the owned cells and the pointer offset
             ! a bit easier. Also abbreviate the pointer for w.

             il = flowDoms(nn,1,1)%il
             jl = flowDoms(nn,1,1)%jl
             kl = flowDoms(nn,1,1)%kl

             p  = IOVar(nn,solID)%pointerOffset
             w => IOVar(nn,solID)%w

             ! Loop over the owned cells and take the possible pointer
             ! offset into account. At the moment the eddy viscosity is
             ! stored in position itu1 of w. Furthermore the laminar
             ! viscosity is not known yet and therefore assume it to be
             ! the free stream value.

             do k=2,kl
               kp = k+p
               do j=2,jl
                 jp = j+p
                 do i=2,il
                   ip = i+p

                   ratio = w(ip,jp,kp,itu1)/muInf
                   nu    = muInf/w(ip,jp,kp,irho)

                   w(ip,jp,kp,itu1) = saNuKnownEddyRatio(ratio, nu)
                 enddo
               enddo
             enddo

           enddo domainLoopCompute

         else checkEddyPresent

           ! No turbulence info is present in the restart file.
           ! Initialize nu tilde to the free stream value.

           domainLoopSet: do nn=1,nDom

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
                   w(ip,jp,kp,itu1) = wInf(itu1)
                 enddo
               enddo
             enddo

           enddo domainLoopSet

           ! Print a warning that nu tilde has been set to the
           ! free stream values. Only processor 0 does this.

           if(myID == 0) then

             print "(a)", "#"
             print "(a)", "#                 Warning"
             print "(a)", "# No turbulent info present in the &
                           &restart file."
             print "(a)", "# Nu tilde for Spalart-Allmaras model has &
                          &been set to the free stream value."
             print "(a)", "#"

           endif

         endif checkEddyPresent

       endif testNuPresent

       end subroutine readTurbSAPlot3D
