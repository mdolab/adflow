!
!      ******************************************************************
!      *                                                                *
!      * File:          readTurbV2f.f90                                 *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 05-05-2004                                      *
!      * Last modified: 09-10-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readTurbV2f(nTypeMismatch)
!
!      ******************************************************************
!      *                                                                *
!      * readTurbV2f reads or constructs the four transport variables   *
!      * for the v2f model. If no information could be retrieved some   *
!      * engineering guess of the turbulent variables is made.          *
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
       integer :: realTypeCGNS, itu
       integer, dimension(4) :: indW

       integer(kind=intType) :: i, j, k, ii, nn, mm, po, ip, jp, kp
       integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd

       real(kind=realType) :: nuScale, kScale, epsScale, fScale

       real(kind=realType), dimension(4) :: turbScale

       character(len=maxCGNSNameLen), dimension(4) :: namesVar
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

       ! Set the names and indices for the four variables.

       indW(1) = itu1; namesVar(1) = cgnsTurbK
       indW(2) = itu2; namesVar(2) = cgnsTurbEpsilon
       indW(3) = itu3; namesVar(3) = cgnsTurbV2
       indW(4) = itu4; namesVar(4) = cgnsTurbF

       ! Compute the scales for nu, k, epsilon and f; v2 has the same
       ! scaling as k.

       nuScale  = muScale/rhoScale
       kScale   = velScale**2
       fScale   = kScale/nuScale
       epsScale = kScale*fScale

       turbScale(1) = kScale
       turbScale(2) = epsScale
       turbScale(3) = kScale
       turbScale(4) = fScale

       ! Loop over the four variables of the v2f model.

       varLoop: do ii=1,4

         ! Find the index of the variable in the solution file and check
         ! if it is present. If not exit the loop.

         mm = nVar
         nn = bsearchStrings(namesVar(ii), varNames, mm)

         if(nn == 0) exit

         ! Variable is present. First determine whether or not a type
         ! mismatch occurs. If so, update nTypeMismatch.

         if(realTypeCGNS /= varTypes(nn)) &
           nTypeMismatch = nTypeMismatch + 1

         ! Read the variable from the restart file and store
         ! it in buffer.

         call readRestartVariable(varNames(nn))

         ! Copy the variables from buffer into w.
         ! Take the possible pointer offset into account.

         itu = indW(ii)

         do k=kBeg,kEnd
           kp = k+po
           do j=jBeg,jEnd
             jp = j+po
             do i=iBeg,iEnd
               ip = i+po
               w(ip,jp,kp,itu) = turbScale(ii)*buffer(i,j,k)
             enddo
           enddo
         enddo

       enddo varLoop

       ! Check if all variables were present. If not, set all turbulence
       ! variables to the free-stream values.

       testPresent: if(ii <= 4) then

         ! Not all variables are present. Set all 4 to the free-stream
         ! values. Take the possible pointer offset into account.

         do ii=nt1,nt2
           do k=kBeg,kEnd
             kp = k+po
             do j=jBeg,jEnd
               jp = j+po
               do i=iBeg,iEnd
                 ip = i+po
                 w(ip,jp,kp,ii) = wInf(ii)
               enddo
             enddo
           enddo
         enddo

         ! Print a warning that the turbulence has been initialized to
         ! the free-stream. Only processor 0 does this for block 1.

         if((myID == 0) .and. (nbkLocal == 1)) then

           print "(a)", "#"
           print "(a)", "#                 Warning"
           print "(a)", "# Not all turbulence variables are present &
                        &for the v2f model."
           print "(a)", "# They have been initialized to the free &
                        &stream values."
           print "(a)", "#"

         endif

       endif testPresent

       end subroutine readTurbV2f
