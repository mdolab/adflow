       subroutine readXvelocity(nTypeMismatch)
!
!       readXvelocity reads the x-velocity variable from the given     
!       place in the cgns file. If the x-velocity itself is not stored 
!       then it is tried to construct it from the x-momentum and       
!       density; it is assumed that the latter is already stored in    
!       the pointer variable w.                                        
!       If it is not possible to create the x-velocity an error        
!       message is printed and the program will stop.                  
!       It is assumed that the pointers in blockPointers already       
!       point to the correct block.                                    
!
       use constants
       use cgnsNames
       use blockPointers, only : w, nbklocal
       use IOModule, only : IOVar
       use restartMod, only : nVar, solID, buffer, varNames, rhoScale, &
            velScale, varTypes
       use utils, only : setCGNSRealType, terminate
       use sorting, only : bsearchStrings
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

       ! Set the cell range to be copied from the buffer.

       iBeg = lbound(buffer,1); iEnd = ubound(buffer,1)
       jBeg = lbound(buffer,2); jEnd = ubound(buffer,2)
       kBeg = lbound(buffer,3); kEnd = ubound(buffer,3)

       ! Set the cgns real type and abbreviate the solution variable and
       ! the pointer offset to improve readability.

       realTypeCGNS = setCGNSRealType()

       po = IOVar(nbkLocal,solID)%pointerOffset
       w => IOVar(nbkLocal,solID)%w

       ! Find out if the x-velocity is present in the solution file.

       mm = nVar
       nn = bsearchStrings(cgnsVelX, varNames, mm)
       if(nn > 0) then

         ! X-velocity is present. First determine whether or not a type
         ! mismatch occurs. If so, update nTypeMismatch.

         if(realTypeCGNS /= varTypes(nn)) &
           nTypeMismatch = nTypeMismatch + 1

         ! Read the x-velocity from the restart file and store it
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
               w(ip,jp,kp,ivx) = buffer(i,j,k)*velScale
             enddo
           enddo
         enddo

         ! X-velocity is read, so a return can be made.

         return

       endif

       ! X-velocity not present. Check for x-momentum.

       nn = bsearchStrings(cgnsMomX, varNames, mm)
       if(nn > 0) then

         ! X-momentum is present. First determine whether or not a type
         ! mismatch occurs. If so, update nTypeMismatch.

         if(realTypeCGNS /= varTypes(nn)) &
           nTypeMismatch = nTypeMismatch + 1

         ! Read the x-momentum from the restart file and store
         ! it in buffer.

         call readRestartVariable(varNames(nn))

         ! Construct the x-velocity; it is assumed that the density is
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
               w(ip,jp,kp,ivx) = buffer(i,j,k)*scale/w(ip,jp,kp,irho)
             enddo
           enddo
         enddo

         ! X-velocity is constructed, so a return can be made.

         return

       endif

       ! Not able to determine the x-velocity.
       ! Print an error message and exit.

       call terminate("readXvelocity", &
                      "Not able to retrieve x-velocity from the &
                      &variables in the restart file.")

       end subroutine readXvelocity
