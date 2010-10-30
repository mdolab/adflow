!
!      ******************************************************************
!      *                                                                *
!      * File:          readZmomentum.f90                               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 05-09-2004                                      *
!      * Last modified: 09-10-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readZmomentum(nTypeMismatch)
!
!      ******************************************************************
!      *                                                                *
!      * readZmomentum reads the z-momentum variable from the given     *
!      * place in the cgns file. If the z-momentum itself is not stored *
!      * then it is tried to construct it from the z-velocity and       *
!      * density; it is assumed that the latter is already stored in    *
!      * the pointer variable w.                                        *
!      * If it is not possible to create the z-velocity an error        *
!      * message is printed and the program will stop.                  *
!      * It is assumed that the pointers in blockPointers already       *
!      * point to the correct block.                                    *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use cgnsNames
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

       real(kind=realType) :: momScale
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

       ! Compute the momentum scaling factor, set the cgns real type and
       ! abbreviate the solution variable and the pointer offset to
       ! improve readability.

       momScale     = rhoScale*velScale
       realTypeCGNS = setCGNSRealType()

       po = IOVar(nbkLocal,solID)%pointerOffset
       w => IOVar(nbkLocal,solID)%w

       ! Find out if the Z-momentum is present in the solution file.

       mm = nVar
       nn = bsearchStrings(cgnsMomZ, varNames, mm)

       testMzPresent: if(nn > 0) then

         ! Z-momentum is present. First determine whether or not a type
         ! mismatch occurs. If so, update nTypeMismatch.

         if(realTypeCGNS /= varTypes(nn)) &
           nTypeMismatch = nTypeMismatch + 1

         ! Read the z-momentum from the restart file and store it in buffer.

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
               w(ip,jp,kp,imz) = buffer(i,j,k)*momScale
             enddo
           enddo
         enddo

         ! Z-momentum is read, so a return can be made.

         return

       endif testMzPresent

       ! Z-momentum is not present. Check for z-velocity.

       nn = bsearchStrings(cgnsVelZ, varNames, mm)

       testVzPresent: if(nn > 0) then

         ! Z-velocity is present. First determine whether or not a type
         ! mismatch occurs. If so, update nTypeMismatch.

         if(realTypeCGNS /= varTypes(nn)) &
           nTypeMismatch = nTypeMismatch + 1

         ! Read the z-velocity from the restart file and store it in buffer.

         call readRestartVariable(varNames(nn))

         ! Copy the variables from buffer into w. Multiply by the
         ! density and velocity scaling factor to obtain to correct
         ! non-dimensional value. Take the possible pointer offset
         ! into account.

         do k=kBeg,kEnd
           kp = k+po
           do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
               ip = i+po
               w(ip,jp,kp,imz) = buffer(i,j,k)*w(ip,jp,kp,irho)*velScale
             enddo
           enddo
         enddo

         ! Z-momentum is constructed, so a return can be made.

         return

       endif testVzPresent

       ! Z-momentum could not be created. Terminate.

       call terminate("readZmomentum", &
                      "Z-Momentum could not be created")

       end subroutine readZmomentum
