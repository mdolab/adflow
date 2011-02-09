!
!      ******************************************************************
!      *                                                                *
!      * File:          readEnergy.f90                                  *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 05-09-2004                                      *
!      * Last modified: 09-10-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readEnergy(nTypeMismatch)
!
!      ******************************************************************
!      *                                                                *
!      * readEnergy reads the energy variable from the given place in   *
!      * the cgns file. If the energy is not stored then it is tried to *
!      * construct it from the pressure, density and velocities. If it  *
!      * is not possible to create the energy an error message is       *
!      * printed and the program will stop. It is assumed that the      *
!      * pointers in blockPointers already point to the correct block.  *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use cgnsNames
       use flowVarRefState
       use IOModule
       use restartMod
       implicit none
!
!      Subroutine argument.
!
       integer(kind=intType), intent(inout) :: nTypeMismatch
!
!      Local variables
!
       integer :: realTypeCGNS

       integer(kind=intType) :: i, j, k, nn, mm, po, ip, jp, kp
       integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd

       real(kind=realType) :: vx, vy, vz, dummyK, pres, rhoInv
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

       ! Find out if the total energy is present in the solution file.

       mm = nVar
       nn = bsearchStrings(cgnsEnergy, varNames, mm)

       testRhoEPresent: if(nn > 0) then

         ! Total energy is present. First determine whether or not a type
         ! mismatch occurs. If so, update nTypeMismatch.

         if(realTypeCGNS /= varTypes(nn)) &
           nTypeMismatch = nTypeMismatch + 1

         ! Read the energy from the restart file and store it in buffer.

         call readRestartVariable(varNames(nn))

         ! Copy the variables from buffer into w. Multiply by the scaling
         ! factor to obtain to correct non-dimensional value and take the
         ! possible pointer offset into account.

         do k=kBeg,kEnd
           kp = k+po
           do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
               ip = i+po
               w(ip,jp,kp,irhoE) = buffer(i,j,k)*pScale
             enddo
           enddo
         enddo

         ! Energy has been read, so a return can be made.

         return

       endif testRhoEPresent

       ! Total energy is not present. Check for the pressure.

       nn = bsearchStrings(cgnsPressure, varNames, mm)

       testPressure: if(nn > 0) then

         ! Pressure is present. First determine whether or not a type
         ! mismatch occurs. If so, update nTypeMismatch.

         if(realTypeCGNS /= varTypes(nn)) &
           nTypeMismatch = nTypeMismatch + 1

         ! Read the pressure from the restart file and store it in buffer.

         call readRestartVariable(varNames(nn))

         ! Compute the total energy. This depends whether or not
         ! a turbulent kinetic energy is present. Take the possible
         ! pointer offset into account.
         ! As this routine is only called to construct the states in
         ! the past for a time accurate computation, the momentum is
         ! stored and not the velocity.

         if( kPresent ) then

           do k=kBeg,kEnd
             kp = k+po
             do j=jBeg,jEnd
               jp = j+po
               do i=iBeg,iEnd
                 ip = i+po
                 rhoInv = one/w(ip,jp,kp,irho)
                 vx     = w(ip,jp,kp,imx)*rhoInv
                 vy     = w(ip,jp,kp,imy)*rhoInv
                 vz     = w(ip,jp,kp,imz)*rhoInv
                 pres   = buffer(i,j,k)*pScale
                 call etotArray(w(ip,jp,kp,irho), vx, vy, vz, pres,  &
                                w(ip,jp,kp,itu1), w(ip,jp,kp,irhoE), &
                                kPresent, 1_intType)
               enddo
             enddo
           enddo

         else

           dummyK = zero

           do k=kBeg,kEnd
             kp = k+po
             do j=jBeg,jEnd
               jp = j+po
               do i=iBeg,iEnd
                 ip = i+po
                 rhoInv = one/w(ip,jp,kp,irho)
                 vx     = w(ip,jp,kp,imx)*rhoInv
                 vy     = w(ip,jp,kp,imy)*rhoInv
                 vz     = w(ip,jp,kp,imz)*rhoInv
                 pres   = buffer(i,j,k)*pScale
                 call etotArray(w(ip,jp,kp,irho), vx, vy, vz, pres, &
                                dummyK, w(ip,jp,kp,irhoE),          &
                                kPresent, 1_intType)
               enddo
             enddo
           enddo

         endif

         ! Energy has been created. So a return can be made.

         return

       endif testPressure

       ! Energy could not be created. Terminate.

       call terminate("readEnergy", &
                      "Energy could not be created")

       end subroutine readEnergy
