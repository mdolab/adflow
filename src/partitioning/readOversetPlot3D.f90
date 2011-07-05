!
!      ******************************************************************
!      *                                                                *
!      * File:          readOversetPlot3D.f90                           *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 06-10-2005                                      *
!      * Last modified: 08-13-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readOversetPlot3D
!
!      ******************************************************************
!      *                                                                *
!      * readOversetPlot3D does nothing yet.                            *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       oversetPresent = .false.

       ! Set the number of holes and overset connectivities to 0 for 
       ! each domain.

       do nn = 1,cgnsNDom
         cgnsDoms(nn)%nOverset      = 0
         cgnsDoms(nn)%nCellsOverset = 0
         cgnsDoms(nn)%nHoles        = 0

         nullify(cgnsDoms(nn)%connOver)
         nullify(cgnsDoms(nn)%hole)

         allocate(cgnsDoms(nn)%connOver(0), &
                  cgnsDoms(nn)%hole(0),     stat=ierr)
         if (ierr /= 0) &
           call terminate("readOversetPlot3D", &
                          "Memory allocation failure for overset data")
       end do

       end subroutine readOversetPlot3D
