!
!      ******************************************************************
!      *                                                                *
!      * File:          readTurbKwTypePlot3D.f90                        *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 07-23-2005                                      *
!      * Last modified: 10-18-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readTurbKwTypePlot3D
!
!      ******************************************************************
!      *                                                                *
!      * readTurbKwTypePlot3D reads or constructs the k and omega       *
!      * values for two-equations turbulence models of the k-omega      *
!      * type. If not all the information could be retrieved some       *
!      * engineering guess of the turbulent variables is made.          *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use cgnsNames
       use communication
       use flowVarRefState
       use inputPhysics
       use inputTimeSpectral
       use IOModule
       use restartMod
       implicit none
!
!      Local variables.
!
       integer :: indW

       integer(kind=intType) :: nn, i, j, k, pp, ip, jp, kp

       real(kind=realType) :: nuScale, kScale, omegaScale, tauScale

       logical :: eddyVisPresent, turbKPresent, omegaPresent
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
       ! Compute the scales for nu, k, omega and tau.

       nuScale    = muScale/rhoScale
       kScale     = velScale**2
       omegaScale = kScale/nuScale
       tauScale   = one/omegaScale

       ! Find the index of turbulent kinetic energy in the solution file
       ! and check if it is present.

       pp = nVar
       nn = bsearchStrings(cgnsTurbK, varNames, pp)

       turbKPresent = .false.
       if(nn > 0) then

         ! Turbulent kinetic energy is present. Set turbKPresent to
         ! .true. for later purposes, determine the index in the solution
         ! file, compute the corresponding offset and read it. Store it
         ! at the correct location, i.e. itu1.

         turbKPresent = .true.
         nn = sorted2Or(nn)
         pp = itu1
         P3D_Offset = sizeHeader + (nn-1)*sizeVolumeSol

         call readPlot3DVar(solID, kScale, pp)

       endif

       ! Find the index of turbulent dissipation rate, i.e. omega, in
       ! the solution file and check if it is present.

       pp = nVar
       nn = bsearchStrings(cgnsTurbOmega, varNames, pp)

       omegaPresent = .false.
       
       if(nn > 0) then

         ! Omega is present. Set omegaPresent to .true. for later 
         ! purposes, compute the corresponding offset and read it. Store
         ! it at the correct location, i.e. itu2.

         omegaPresent = .true.
         nn = sorted2Or(nn)
         pp = itu2
         P3D_Offset = sizeHeader + (nn-1)*sizeVolumeSol

         call readPlot3DVar(solID, omegaScale, pp)
       endif

       ! If omega is not present, check if tau is present and
       ! initialize omega accordingly.

       testForTau: if(.not. omegaPresent) then

         pp = nVar
         nn = bsearchStrings(cgnsTurbTau, varNames, pp)

         testTauPresent: if(nn > 0) then

           ! Tau is present. Set omegaPresent to .true. for later 
           ! purposes, compute the corresponding offset and read it.
           ! Store it at the correct location, i.e. itu2.

           omegaPresent = .true.
           nn = sorted2Or(nn)
           pp = itu2
           P3D_Offset = sizeHeader + (nn-1)*sizeVolumeSol

           call readPlot3DVar(solID, tauScale, pp)

           ! Transform tau to omega.

           domainLoopTransform: do nn=1,nDom

             ! Store the end of the owned cells and the pointer offset
             ! a bit easier. Also abbreviate the pointer for w.

             il = flowDoms(nn,1,1)%il
             jl = flowDoms(nn,1,1)%jl
             kl = flowDoms(nn,1,1)%kl

             pp = IOVar(nn,solID)%pointerOffset
             w => IOVar(nn,solID)%w

             do k=2,kl
               kp = k+pp
               do j=2,jl
                 jp = j+pp
                 do i=2,il
                   ip = i+pp
                   w(ip,jp,kp,itu2) = one/w(ip,jp,kp,itu2)
                 enddo
               enddo
             enddo

           enddo domainLoopTransform

         endif testTauPresent

       endif testForTau

       ! Check if both variables were present.
       ! If so go to the check to transform omega to tau.

       if(turbKPresent .and. omegaPresent) goto 1001

       ! K and omega are not both present. It is tried to construct
       ! their values with the information that is present.

       ! Try to construct the eddy viscosity. The place to store these
       ! values depend whether k or omega is present.

       indW = itu1
       if( turbKPresent ) indW = itu2

       call readTurbEddyVisPlot3D(indW, eddyVisPresent)

       ! The eddy viscosity is either known or still initialized
       ! to the free stream value. In any case determine the
       ! situation we are dealing with and try to initialize k and
       ! omega accordingly.

       checkKPresent: if( turbKPresent ) then

         domainLoopKPresent: do nn=1,nDom

           ! K is present. Store the end of the owned cells and the
           ! pointer offset a bit easier. Also abbreviate the pointer
           ! for w.

           il = flowDoms(nn,1,1)%il
           jl = flowDoms(nn,1,1)%jl
           kl = flowDoms(nn,1,1)%kl

           pp = IOVar(nn,solID)%pointerOffset
           w => IOVar(nn,solID)%w

           ! Compute omega using the eddy viscosity.
           ! Assume that the standard k-omega formula is also valid
           ! for the SST-model.

           do k=2,kl
             kp = k+pp
             do j=2,jl
               jp = j+pp
               do i=2,il
                 ip = i+pp
                 w(ip,jp,kp,itu2) = w(ip,jp,kp,irho)*w(ip,jp,kp,itu1) &
                                  / w(ip,jp,kp,itu2)
               enddo
             enddo
           enddo

         enddo domainLoopKPresent

         ! Print a warning that omega was not present and has been
         ! constructed. Only processor 0 does this.

         if(myID == 0) then

           print "(a)", "#"
           print "(a)", "#                 Warning"
           print "(a)", "# Omega is not present in the restart file."
           if( eddyVisPresent ) then
             print "(a)", "# It is initialized using the turbulent &
                          &kinetic energy and eddy viscosity."
           else
             print "(a)", "# It is initialized using the turbulent &
                          &kinetic energy and free stream eddy &
                          &viscosity."
           endif
           print "(a)", "#"

         endif

         ! K and omega are initialized.
         ! Go to the check to transform omega to tau.

         goto 1001

       endif checkKPresent

       ! K is not present. Check if omega is present and compute
       ! k accordingly.

       checkOmegaPresent: if( omegaPresent ) then

         domainLoopOmegaPresent: do nn=1,nDom

           ! Omega is present. Store the end of the owned cells and the
           ! pointer offset a bit easier. Also abbreviate the pointer
           ! for w.

           il = flowDoms(nn,1,1)%il
           jl = flowDoms(nn,1,1)%jl
           kl = flowDoms(nn,1,1)%kl

           pp = IOVar(nn,solID)%pointerOffset
           w => IOVar(nn,solID)%w

           ! Compute k using the eddy viscosity.
           ! Assume that the standard k-omega formula is also valid
           ! for the SST-model.

           do k=2,kl
             kp = k+pp
             do j=2,jl
               jp = j+pp
               do i=2,il
                 ip = i+pp
                 w(ip,jp,kp,itu1) = w(ip,jp,kp,itu1)*w(ip,jp,kp,itu2) &
                                  / w(ip,jp,kp,irho)
               enddo
             enddo
           enddo

         enddo domainLoopOmegaPresent

         ! Print a warning that k was not present and has been
         ! constructed. Only processor 0 does this.

         if(myID == 0) then

           print "(a)", "#"
           print "(a)", "#                 Warning"
           print "(a)", "# Turbulent kinetic energy is not present &
                        &in the restart file."
           if( eddyVisPresent ) then
             print "(a)", "# It is initialized using omega and &
                          &the eddy viscosity."
           else
             print "(a)", "# It is initialized using omega and &
                          &the free stream eddy viscosity."
           endif
           print "(a)", "#"

         endif

         ! K and omega are initialized.
         ! Go to the check to transform omega to tau.

         goto 1001

       endif checkOmegaPresent

       ! Both k and omega are not present. Use a guess for omega
       ! and compute k using the known value of the eddy viscosity.

       domainLoopGuess: do nn=1,nDom

         ! Set the pointers for this block. Make sure that no bounds
         ! are exceeded in case a spectral interpolation must be done.

         pp = min(solID, nTimeIntervalsSpectral)
         call setPointers(nn, 1_intType, pp)

         ! Reset the pointer for w to IOVar.

         w => IOVar(nn,solID)%w

         ! The laminar viscosity is not yet known, so set it to the
         ! free stream value.

         rlv = muInf
         call initKOmega(IOVar(nn,solID)%pointerOffset)

       enddo domainLoopGuess

       ! Print a warning that both k and omega are not present in
       ! the restart file. Only processor 0 does this.

       if(myID == 0) then

         print "(a)", "#"
         print "(a)", "#                 Warning"
         print "(a)", "# The turbulent kinetic energy and omega are &
                      &not present in the restart file."
         if( eddyVisPresent ) then
           print "(a)", "# They have been initialized using the &
                        &eddy viscosity."
         else
           print "(a)", "# The default initialization has been used."
         endif

         print "(a)", "#"

       endif

       ! For the k-tau model omega must be transformed to tau.

 1001  select case (turbModel)

         case (ktau)

           domainLoopTau: do nn=1,nDom

             ! Store the end of the owned cells and the pointer offset a
             ! bit easier. Also abbreviate the pointer for w.

             il = flowDoms(nn,1,1)%il
             jl = flowDoms(nn,1,1)%jl
             kl = flowDoms(nn,1,1)%kl

             pp = IOVar(nn,solID)%pointerOffset
             w => IOVar(nn,solID)%w

             ! Invert omega to tau.

             do k=2,kl
               kp = k+pp
               do j=2,jl
                 jp = j+pp
                 do i=2,il
                   ip = i+pp
                   w(ip,jp,kp,itu2) = one/w(ip,jp,kp,itu2)
                 enddo
               enddo
             enddo

           enddo domainLoopTau

       end select

       end subroutine readTurbKwTypePlot3D
