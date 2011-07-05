!
!      ******************************************************************
!      *                                                                *
!      * File:          readRestartVariable.F90                         *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-20-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readRestartVariable(cgnsVarName)
!
!      ******************************************************************
!      *                                                                *
!      * readRestartVariable reads the given variable name from the     *
!      * cgns restart file.                                             *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use su_cgns
       use restartMod
       implicit none
!
!      Subroutine arguments.
!
       character(len=*), intent(in) :: cgnsVarName
!
!      Local variables.
!
       integer :: ierr, realTypeCGNS

       integer(kind=intType) :: i, j, k
!
!      Function definition.
!
       integer :: setCGNSRealType
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
#ifdef USE_NO_CGNS

       call terminate("readRestartVariable", &
                      "Routine should not be called if no cgns support &
                      &is selected.")

#else
       ! Set the cgns real type.

       realTypeCGNS = setCGNSRealType()

       ! Check where the solution variables are stored.

       locationTest: if(location == CellCenter) then

         ! Cell centered values. Read the values directly in the buffer.

         call cg_field_read_f(cgnsInd, cgnsBase, cgnsZone, cgnsSol, &
                              cgnsVarName, realTypeCGNS, rangeMin,  &
                              rangeMax, buffer, ierr)
         if(ierr /= all_ok)                        &
           call terminate("readRestartVariable", &
                          "Something wrong when calling cg_field_read_f")
       else locationTest

         ! Vertex centered values. First read the solution in the
         ! array bufferVertex.

         call cg_field_read_f(cgnsInd, cgnsBase, cgnsZone, cgnsSol, &
                              cgnsVarName, realTypeCGNS, rangeMin, &
                              rangeMax, bufferVertex, ierr)
         if(ierr /= all_ok)                        &
           call terminate("readRestartVariable", &
                          "Something wrong when calling cg_field_read_f")

         ! Create the cell centered values by averaging the vertex values.

         do k=2,kl
           do j=2,jl
             do i=2,il
               buffer(i,j,k) = eighth*(bufferVertex(i-1,j-1,k-1) &
                             +         bufferVertex(i,  j-1,k-1) &
                             +         bufferVertex(i-1,j,  k-1) &
                             +         bufferVertex(i,  j,  k-1) &
                             +         bufferVertex(i-1,j-1,k)   &
                             +         bufferVertex(i,  j-1,k)   &
                             +         bufferVertex(i-1,j,  k)   &
                             +         bufferVertex(i,  j,  k))
             enddo
           enddo
         enddo

       endif locationTest

#endif

       end subroutine readRestartVariable
