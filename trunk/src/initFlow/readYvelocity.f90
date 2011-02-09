!
!      ******************************************************************
!      *                                                                *
!      * File:          readYvelocity.f90                               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-21-2003                                      *
!      * Last modified: 10-07-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readYvelocity(nTypeMismatch)
!
!      ******************************************************************
!      *                                                                *
!      * readYvelocity reads the y-velocity variable from the given     *
!      * place in the cgns file. If the y-velocity itself is not stored *
!      * then it is tried to construct it from the y-momentum and       *
!      * density; it is assumed that the latter is already stored in    *
!      * the pointer variable w.                                        *
!      * If it is not possible to create the y-velocity an error        *
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

       real(kind=realType) :: scale
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

       ! Find out if the y-velocity is present in the solution file.

       mm = nVar
       nn = bsearchStrings(cgnsVelY, varNames, mm)
       if(nn > 0) then

         ! Y-velocity is present. First determine whether or not a type
         ! mismatch occurs. If so, update nTypeMismatch.

         if(realTypeCGNS /= varTypes(nn)) &
           nTypeMismatch = nTypeMismatch + 1

         ! Read the y-velocity from the restart file and store it
         ! in buffer.

         call readRestartVariable(varNames(nn))

         ! Copy the variables from buffer into w. Multiply by the scale
         ! factor to obtain the correct nondimensional value and take
         ! the possible pointer offset into account.

         do k=kBeg,kEnd
           kp = k+po
           do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
               ip = i+po
               w(ip,jp,kp,ivy) = buffer(i,j,k)*velScale
             enddo
           enddo
         enddo

         ! Y-velocity is read, so a return can be made.

         return

       endif

       ! Y-velocity not present. Check for y-momentum.

       nn = bsearchStrings(cgnsMomY, varNames, mm)
       if(nn > 0) then

         ! Y-momentum is present. First determine whether or not a type
         ! mismatch occurs. If so, update nTypeMismatch.

         if(realTypeCGNS /= varTypes(nn)) &
           nTypeMismatch = nTypeMismatch + 1

         ! Read the y-momentum from the restart file and store
         ! it in buffer.

         call readRestartVariable(varNames(nn))

         ! Construct the y-velocity; it is assumed that the density is
         ! already stored in w. Multiply by the momentum scale factor
         ! to obtain the correct nondimensional value and take the
         ! possible pointer offset into account.

         scale = rhoScale*velScale

         do k=kBeg,kEnd
           kp = k+po
           do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
               ip = i+po
               w(ip,jp,kp,ivy) = buffer(i,j,k)*scale/w(ip,jp,kp,irho)
             enddo
           enddo
         enddo

         ! Y-velocity is constructed, so a return can be made.

         return

       endif

       ! Not able to determine the y-velocity.
       ! Print an error message and exit.

       call terminate("readYvelocity", &
                      "Not able to retrieve y-velocity from the &
                      &variables in the restart file.")

       end subroutine readYvelocity
