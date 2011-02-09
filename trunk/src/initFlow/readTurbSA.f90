!
!      ******************************************************************
!      *                                                                *
!      * File:          readTurbSA.F90                                  *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 05-05-2004                                      *
!      * Last modified: 09-10-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readTurbSA(nTypeMismatch)
!
!      ******************************************************************
!      *                                                                *
!      * readTurbSA reads or constructs the nu tilde transport          *
!      * variable for the Spalart-Allmaras type turbulence models.      *
!      * If no information could be retrieved some engineering guess of *
!      * the turbulent variables is made.                               *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use cgnsNames
       use constants
       use communication
       use flowVarRefState
       use IOModule
       use restartMod
       implicit none
!
!      Subroutine argument.
!
       integer(kind=intType), intent(inout) :: nTypeMismatch
!
!      Local variables.
!
       integer :: realTypeCGNS

       integer(kind=intType) :: i, j, k, nn, mm, po, ip, jp, kp
       integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd

       real(kind=realType) :: nuScale, ratio, nu

       logical :: eddyVisPresent
!
!      Function definitions.
!
       integer               :: setCGNSRealType
       integer(kind=intType) :: bsearchStrings

       real(kind=realType) :: saNuKnownEddyRatio
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set the cell range to be copied from the buffer.

       iBeg = lbound(buffer,1); iEnd = ubound(buffer,1)
       jBeg = lbound(buffer,2); jEnd = ubound(buffer,2)
       kBeg = lbound(buffer,3); kEnd = ubound(buffer,3)

       ! Set the cgns real type and abbreviate the solution variable and 
       ! the pointer offset to improve readability.
       ! Also compute the kinematic viscosity scale.

       realTypeCGNS = setCGNSRealType()

       po = IOVar(nbkLocal,solID)%pointerOffset
       w => IOVar(nbkLocal,solID)%w

       nuScale = muScale/rhoScale

       ! Check if the nu tilde variable is present.

       mm = nVar
       nn = bsearchStrings(cgnsTurbSANu, varNames, mm)

       nuTildePresent: if(nn > 0) then

         ! Nu tilde is present. First determine whether or not a type
         ! mismatch occurs. If so, update nTypeMismatch.

         if(realTypeCGNS /= varTypes(nn)) &
           nTypeMismatch = nTypeMismatch + 1

         ! Read nu tilde from the restart file and store
         ! it in buffer.

         call readRestartVariable(varNames(nn))

         ! Copy the variables from buffer into w and take the possible
         ! pointer offset into account.

         do k=kBeg,kEnd
           kp = k+po
           do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
               ip = i+po
               w(ip,jp,kp,itu1) = nuScale*buffer(i,j,k)
             enddo
           enddo
         enddo

         ! Variable is read, so a return can be made.

         return

       endif nuTildePresent

       ! NuTilde is not present. Try to construct the eddy viscosity.

       call readTurbEddyVis(nTypeMismatch, eddyVisPresent)

       ! Check if the eddy viscosity has been constructed.

       eddyPresent: if( eddyVisPresent ) then

         ! Eddy viscosity is present. As the laminar viscosity is not
         ! yet known, set it to the free-stream value.

         rlv = muInf

         ! Compute nuTilde from the known ratio of eddy and laminar
         ! viscosity. Take the possible pointer offset into account.

         do k=kBeg,kEnd
           kp = k+po
           do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
               ip = i+po

               ! Compute the eddy viscosity ratio and the laminar
               ! kinematic viscosity and call the function to
               ! compute the nu tilde variable.

               ratio            = rev(i,j,k)/rlv(i,j,k)
               nu               = rlv(i,j,k)/w(ip,jp,kp,irho)
               w(ip,jp,kp,itu1) = saNuKnownEddyRatio(ratio, nu)

             enddo
           enddo
         enddo

         ! Print a warning that nu tilde has been constructed and
         ! not read. Only processor 0 does this for block 1.

         if((myID == 0) .and. (nbkLocal == 1)) then

           print "(a)", "#"
           print "(a)", "#                 Warning"
           print "(a)", "# Nu tilde for Spalart-Allmaras model not &
                        &present in the restart file."
           print "(a)", "# Variable has been reconstructed from &
                        &the eddy viscosity ratio."
           print "(a)", "#"

         endif

         ! Variable is constructed, so a return can be made.

         return

       endif eddyPresent

       ! No turbulence info is present in the restart file.
       ! Initialize nu tilde to the free stream value.
       ! Take the possible pointer offset into account.

       do k=kBeg,kEnd
           kp = k+po
           do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
             ip = i+po
             w(ip,jp,kp,itu1) = wInf(itu1)
           enddo
         enddo
       enddo

       ! Print a warning that nu tilde has been set to the
       ! free stream values. Only processor 0 does this for block 1.

       if((myID == 0) .and. (nbkLocal == 1)) then

         print "(a)", "#"
         print "(a)", "#                 Warning"
         print "(a)", "# No turbulent info present in the restart file."
         print "(a)", "# Nu tilde for Spalart-Allmaras model has &
                      &been set to the free stream value."
         print "(a)", "#"

       endif

       end subroutine readTurbSA
