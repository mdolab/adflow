!
!      ******************************************************************
!      *                                                                *
!      * File:          readDensity.f90                                 *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-18-2003                                      *
!      * Last modified: 09-10-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readDensity(nTypeMismatch)
!
!      ******************************************************************
!      *                                                                *
!      * readDensity reads the density from the given place in the      *
!      * cgns file. If the density itself is not stored (unlikely),     *
!      * then it is tried to construct the density from other           *
!      * variables. If this is not possible an error message is printed *
!      * and the program will stop.                                     *
!      * It is assumed that the pointers in blockPointers already       *
!      * point to the correct block.                                    *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use cgnsNames
       use constants
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

       ! Find out if the density is present in the solution file.

       mm = nVar
       nn = bsearchStrings(cgnsDensity, varNames, mm)
       if(nn > 0) then

         ! Density is present. First determine whether or not a type
         ! mismatch occurs. If so, update nTypeMismatch.

         if(realTypeCGNS /= varTypes(nn)) &
           nTypeMismatch = nTypeMismatch + 1

         ! Read the density from the restart file and store it in buffer.

         call readRestartVariable(varNames(nn))

         ! Copy the variables from buffer into w. Multiply by the
         ! scaling factor to obtain to correct nondimensional value and
         ! take the possible pointer offset into account.

         do k=kBeg,kEnd
           kp = k+po
           do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
               ip = i+po
               w(ip,jp,kp,irho) = buffer(i,j,k)*rhoScale
             enddo
           enddo
         enddo

         ! Density is read, so a return can be made.

         return

       endif

       ! Not able to determine the density.
       ! Print an error message and exit.

       call terminate("readDensity", &
                      "Not able to retrieve density from the &
                      &variables in the restart file.")

       end subroutine readDensity
