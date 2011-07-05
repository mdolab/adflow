!
!      ******************************************************************
!      *                                                                *
!      * File:          readTurbV2fPlot3D.f90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 07-23-2005                                      *
!      * Last modified: 10-18-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readTurbV2fPlot3D
!
!      ******************************************************************
!      *                                                                *
!      * readTurbV2fPlot3D reads or constructs the four transport       *
!      * variables for the v2f model. At the moment it is not possible  *
!      * to restart from a different turbulence model for v2f, i.e. all *
!      * 4 variables need to be present. Otherwise they will be set to  *
!      * free stream values.                                            *
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
       integer, dimension(4) :: indW

       integer(kind=intType) :: ii, nn, i, j, k, p, ip, jp, kp
       integer(kind=intType) :: il, jl, kl

       real(kind=realType) :: nuScale, kScale, fScale, epsScale
       real(kind=realType), dimension(4) :: turbScale
       real(kind=realType), dimension(:,:,:,:), pointer :: w

       character(len=maxCGNSNameLen), dimension(4) :: namesVar
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
       ! Set the names and indices for the four variables.

       indW(1) = itu1; namesVar(1) = cgnsTurbK
       indW(2) = itu2; namesVar(2) = cgnsTurbEpsilon
       indW(3) = itu3; namesVar(3) = cgnsTurbV2
       indW(4) = itu4; namesVar(4) = cgnsTurbF

       ! Also set the scaling for the 4 variables; v2 has the same
       ! scaling as k.

       nuScale  = muScale/rhoScale
       kScale   = velScale**2
       fScale   = kScale/nuScale
       epsScale = kScale*fScale

       turbScale(1) = kScale
       turbScale(2) = epsScale
       turbScale(3) = kScale
       turbScale(4) = fScale

       ! Loop over the 4 variables and read their values.

       varLoop: do ii=1,4

         ! Find the index of the variable in the solution file and check
         ! if it is present. If not exit the loop.

         p  = nVar
         nn = bsearchStrings(namesVar(ii), varNames, p)

         if(nn == 0) exit

         ! Transport variable is present. Determine the index in the
         ! solution file, the corresponding offset and read it. Store
         ! it at the correct location.

         nn = sorted2Or(nn)
         p  = indW(ii)
         P3D_Offset = sizeHeader + (nn-1)*sizeVolumeSol

         call readPlot3DVar(solID, turbScale(ii), p)
       enddo varLoop

       ! Check if all variables were present. If not, set all turbulence
       ! variables to the free-stream values.

       testPresent: if(ii <= 4) then

         ! Not all variables were present.
         ! Loop over the local number of blocks.

         domainLoop: do nn=1,nDom

           ! Store the end of the owned cells and the pointer offset
           ! a bit easier. Also abbreviate the pointer for w.

           il = flowDoms(nn,1,1)%il
           jl = flowDoms(nn,1,1)%jl
           kl = flowDoms(nn,1,1)%kl

           p  = IOVar(nn,solID)%pointerOffset
           w => IOVar(nn,solID)%w

           ! Loop over the owned cells and reset the turbulence variables
           ! to the free stream values. Take the possible pointer
           ! offset into account.

           do ii=nt1, nt2
             do k=2,kl
               kp = k+p
               do j=2,jl
                 jp = j+p
                 do i=2,il
                   ip = i+p
                   w(ip,jp,kp,ii) = wInf(ii)
                 enddo
               enddo
             enddo
           enddo

         enddo domainLoop

         ! Print a warning that the turbulence has been initialized to
         ! the free stream. Only processor 0 does this.

         if(myID == 0) then

           print "(a)", "#"
           print "(a)", "#                 Warning"
           print "(a)", "# Not all turbulence variables are present &
                        &for the v2f model."
           print "(a)", "# All turbulence variables have been &
                        &initialized to the free stream values."
           print "(a)", "#"

         endif

       endif testPresent

       end subroutine readTurbV2fPlot3D
