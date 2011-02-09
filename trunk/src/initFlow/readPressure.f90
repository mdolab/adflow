!
!      ******************************************************************
!      *                                                                *
!      * File:          readPressure.f90                                *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-21-2003                                      *
!      * Last modified: 11-20-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readPressure(nTypeMismatch)
!
!      ******************************************************************
!      *                                                                *
!      * readPressure reads the pressure variable from the given place  *
!      * in the cgns file. If the pressure itself is not present it is  *
!      * tried to construct if from other variables. In that case it is *
!      * assumed that the density, velocity and turbulent variables are *
!      * already stored in the pointer variable w.                      *
!      * If it is not possible to create the pressure an error message  *
!      * is printed and the program will stop.                          *
!      * It is assumed that the pointers in blockPointers already       *
!      * point to the correct block.                                    *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use cgnsNames
       use inputPhysics
       use IOModule
       use flowVarRefState
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

       ! Find out if the pressure is present in the solution file.

       mm = nVar
       nn = bsearchStrings(cgnsPressure, varNames, mm)
       if(nn > 0) then

         ! Pressure is present. First determine whether or not a type
         ! mismatch occurs. If so, update nTypeMismatch.

         if(realTypeCGNS /= varTypes(nn)) &
           nTypeMismatch = nTypeMismatch + 1

         ! Read the pressure from the restart file and store
         ! it in buffer.

         call readRestartVariable(varNames(nn))

         ! Copy the variables from buffer into the position of rhoE
         ! in w. Multiply by the pressure scale factor to obtain the
         ! correct nondimensional value and take the possible pointer
         ! offset into account.

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

         ! Pressure is read, so a return can be made.

         return

       endif

       ! Pressure is not present. Check for the total energy.

       nn = bsearchStrings(cgnsEnergy, varNames, mm)
       if(nn > 0) then

         ! Total energy is present. First determine whether or not a type
         ! mismatch occurs. If so, update nTypeMismatch.

         if(realTypeCGNS /= varTypes(nn)) &
           nTypeMismatch = nTypeMismatch + 1

         ! Read the total energy from the restart file and store
         ! it in buffer.

         call readRestartVariable(varNames(nn))

         ! Copy the variables from buffer into w. Multiply by the
         ! pressure scale factor to obtain the correct nondimensional
         ! value and take the possible pointer offset into account.

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

         ! Compute the pressure from energy, density and velocities.
         ! This will still be stored in the irhoE position of w.

         call computePressure(iBeg,iEnd,jBeg,jEnd,kBeg,kEnd,po)

         ! Pressure is constructed, so a return can be made.

         return

       endif

       ! Not able to determine the pressure.
       ! Print an error message and exit.

       call terminate("readPressure", &
                      "Not able to retrieve the pressure from &
                      &the variables in the restart file.")

       end subroutine readPressure
