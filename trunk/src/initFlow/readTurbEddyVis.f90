!
!      ******************************************************************
!      *                                                                *
!      * File:          readTurbEddyVis.f90                             *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 05-05-2004                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readTurbEddyVis(nTypeMismatch, eddyVisPresent)
!
!      ******************************************************************
!      *                                                                *
!      * readTurbEddyVis tries to read the eddy viscosity from the      *
!      * restart file.                                                  *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use cgnsNames
       use constants
       use flowVarRefState
       use inputPhysics
       use restartMod
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(inout) :: nTypeMismatch
       logical, intent(out)                 :: eddyVisPresent
!
!      Local variables.
!
       integer :: realTypeCGNS

       integer(kind=intType) :: i, j, k, nn, mm
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
       ! Set the cgns real type

       realTypeCGNS = setCGNSRealType()

       ! Check if the eddy viscosity is present. If so, read it.

       mm = nVar
       nn = bsearchStrings(cgnsEddy, varNames, mm)

       if(nn > 0) then

         ! Eddy viscosity is present. Determine if a type mismatch
         ! occured, read the buffer from the file and set
         ! eddyVisPresent to .true.

         if(realTypeCGNS /= varTypes(nn)) &
           nTypeMismatch = nTypeMismatch + 1

         call readRestartVariable(varNames(nn))

         eddyVisPresent = .true.

         ! Scale the eddy viscosity such that it contains the
         ! correct nonDimensional value.

         do k=2,kl
           do j=2,jl
             do i=2,il
               rev(i,j,k) = muScale*buffer(i,j,k)
             enddo
           enddo
         enddo

         ! Eddy viscosity has been read, so make a return.

         return

       endif

       ! Eddy viscosity is not present. Check if the eddy viscosity
       ! ratio is present. If so read it and construct the eddy
       ! viscosity from it.

       nn = bsearchStrings(cgnsEddyRatio, varNames, mm)

       if(nn > 0) then

         ! Eddy viscosity ratio is present. Determine if a type
         ! mismatch occured, read the buffer from the file and set
         ! eddyVisPresent to .true.

         if(realTypeCGNS /= varTypes(nn)) &
           nTypeMismatch = nTypeMismatch + 1

         call readRestartVariable(varNames(nn))

         eddyVisPresent = .true.

         ! Multiply the eddy viscosity by the laminar viscosity such
         ! that it contains the correct nonDimensional value.
         ! As the laminar viscosity is not yet know, use the free
         ! stream viscosity.

         do k=2,kl
           do j=2,jl
             do i=2,il
               rev(i,j,k) = muInf*buffer(i,j,k)
             enddo
           enddo
         enddo

         ! Eddy viscosity has been read, so make a return.

         return

       endif

       ! Eddy viscosity cannot be constructed. Set it to the
       ! free stream eddy viscosity.

       do k=2,kl
         do j=2,jl
           do i=2,il
             rev(i,j,k) = muInf*eddyVisInfRatio
           enddo
         enddo
       enddo

       ! Eddy viscosity is not present, so set it to .false.

       eddyVisPresent = .false.

       end subroutine readTurbEddyVis
