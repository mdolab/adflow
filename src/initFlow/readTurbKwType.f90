!
!      ******************************************************************
!      *                                                                *
!      * File:          readTurbKwType.f90                              *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 05-05-2004                                      *
!      * Last modified: 09-10-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readTurbKwType(nTypeMismatch)
!
!      ******************************************************************
!      *                                                                *
!      * readTurbKwType reads or constructs the k and omega values      *
!      * for two-equations turbulence models of the k-omega type.       *
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
       use inputPhysics
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

       real(kind=realType) :: nuScale, kScale, omegaScale, val

       logical :: turbKPresent, omegaPresent, eddyVisPresent
!
!      Function definitions.
!
       integer               :: setCGNSRealType
       integer(kind=intType) :: bsearchStrings
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

       realTypeCGNS = setCGNSRealType()

       po = IOVar(nbkLocal,solID)%pointerOffset
       w => IOVar(nbkLocal,solID)%w

       ! Compute the scales for nu, k and omega.

       nuScale    = muScale/rhoScale
       kScale     = velScale**2
       omegaScale = kScale/nuScale

       ! First check if k is present.

       turbKPresent = .false.

       mm = nVar
       nn = bsearchStrings(cgnsTurbK, varNames, mm)

       if(nn > 0) then

         ! K is present. First determine whether or not a type
         ! mismatch occurs. If so, update nTypeMismatch.

         if(realTypeCGNS /= varTypes(nn)) &
           nTypeMismatch = nTypeMismatch + 1

         ! Read k from the restart file and store it in buffer.

         call readRestartVariable(varNames(nn))

         ! Copy the variables from buffer into w. Multiply by the scale
         ! factor to obtain the correct non-dimensional value and take
         ! the possible pointer offset into account.

         do k=kBeg,kEnd
           kp = k+po
           do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
               ip = i+po
               w(ip,jp,kp,itu1) = kScale*buffer(i,j,k)
             enddo
           enddo
         enddo

         ! Set turbKPresent to .true.

         turbKPresent = .true.

       endif

       ! Check if omega is present.

       omegaPresent = .false.

       nn = bsearchStrings(cgnsTurbOmega, varNames, mm)

       if(nn > 0) then

         ! Omega is present. First determine whether or not a type
         ! mismatch occurs. If so, update nTypeMismatch.

         if(realTypeCGNS /= varTypes(nn)) &
           nTypeMismatch = nTypeMismatch + 1

         ! Read omega from the restart file and store it in buffer.

         call readRestartVariable(varNames(nn))

         ! Copy the variables from buffer into w. Multiply by the scale
         ! factor to obtain the correct non-dimensional value and take
         ! the possible pointer offset into account.

         do k=kBeg,kEnd
           kp = k+po
           do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
               ip = i+po
               w(ip,jp,kp,itu2) = omegaScale*buffer(i,j,k)
             enddo
           enddo
         enddo

         ! Set omegaPresent to .true.

         omegaPresent = .true.

       endif

       ! If omega is not present, check if tau is present and
       ! initialize omega accordingly.

       if(.not. omegaPresent) then

         nn = bsearchStrings(cgnsTurbTau, varNames, mm)

         if(nn > 0) then

           ! Tau is present. First determine whether or not a type
           ! mismatch occurs. If so, update nTypeMismatch.

           if(realTypeCGNS /= varTypes(nn)) &
             nTypeMismatch = nTypeMismatch + 1

           ! Read tau from the restart file and store it in buffer.

           call readRestartVariable(varNames(nn))

           ! Transform tau to omega and copy the variables from buffer
           ! into w. Multiply by the scale factor to obtain the correct
           ! non-dimensional value and take the possible pointer offset
           ! into account.

           do k=kBeg,kEnd
             kp = k+po
             do j=jBeg,jEnd
               jp = j+po
               do i=iBeg,iEnd
                 ip = i+po

                 val = buffer(i,j,k)
                 w(ip,jp,kp,itu2) = omegaScale/max(eps,val)
               enddo
             enddo
           enddo

           ! Set omegaPresent to .true.

           omegaPresent = .true.

         endif

       endif

       ! Check if both variables were present.
       ! If so go to the check to transform omega to tau.

       if(turbKPresent .and. omegaPresent) goto 1001

       ! K and omega are not both present. It is tried to construct
       ! their values with the information that is present.

       ! Try to read the eddy viscosity.

       call readTurbEddyVis(nTypeMismatch, eddyVisPresent)

       ! The eddy viscosity is either known or still initialized
       ! to the free stream value. In any case determine the
       ! situation we are dealing with and try to initialize k and
       ! omega accordingly.

       if( turbKPresent ) then

         ! K is present. Compute omega using the eddy viscosity.
         ! Assume that the standard k-omega formula is also valid
         ! for the SST-model.
         ! Take the possible pointer offset into account.

         do k=kBeg,kEnd
           kp = k+po
           do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
               ip = i+po
               w(ip,jp,kp,itu2) = w(ip,jp,kp,irho)*w(ip,jp,kp,itu1) &
                                / rev(i,j,k)
             enddo
           enddo
         enddo

         ! Print a warning that omega was not present and has been
         ! constructed. Only processor 0 does this for block 1.

         if((myID == 0) .and. (nbkLocal == 1)) then

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

       endif

       if( omegaPresent ) then

         ! Omega is present. Compute k using the eddy viscosity.
         ! Assume that the standard k-omega formula is also valid
         ! for the SST-model.
         ! Take the possible pointer offset into account.

         do k=kBeg,kEnd
           kp = k+po
           do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
               ip = i+po
               w(ip,jp,kp,itu1) = rev(i,j,k)*w(ip,jp,kp,itu2) &
                                / w(ip,jp,kp,irho)
             enddo
           enddo
         enddo

         ! Print a warning that k was not present and has been
         ! constructed. Only processor 0 does this for block 1.

         if((myID == 0) .and. (nbkLocal == 1)) then

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

       endif

       ! Both k and omega are not present. Use a guess for omega
       ! and compute k using the known value of the eddy viscosity.
       ! As the laminar viscosity is not yet known, set it to the
       ! free-stream value.

       rlv = muInf
       call initKOmega(po)

       ! Print a warning that both k and omega are not present in
       ! the restart file. Only processor 0 does this for block 1.

       if((myID == 0) .and. (nbkLocal == 1)) then

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
       ! Take the possible pointer offset into account.

 1001  select case (turbModel)

         case (ktau)

           do k=kBeg,kEnd
             kp = k+po
             do j=jBeg,jEnd
               jp = j+po
               do i=iBeg,iEnd
                 ip = i+po
                 w(ip,jp,kp,itu2) = one/w(ip,jp,kp,itu2)
               enddo
             enddo
           enddo

       end select

       end subroutine readTurbKwType
