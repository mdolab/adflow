!
!      ******************************************************************
!      *                                                                *
!      * File:          nullifyCGNSDomPointers.f90                      *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 08-13-2005                                      *
!      * Last modified: 11-08-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine nullifyCGNSDomPointers(nn)
!
!      ******************************************************************
!      *                                                                *
!      * nullifyCGNSDomPointers nullifies all the pointers of the       *
!      * given CGNS block.                                              *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nn
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       nullify(cgnsDoms(nn)%procStored)
       nullify(cgnsDoms(nn)%conn1to1)
       nullify(cgnsDoms(nn)%connNonMatchAbutting)
       nullify(cgnsDoms(nn)%connOver)
       nullify(cgnsDoms(nn)%hole)
       nullify(cgnsDoms(nn)%bocoInfo)

       end subroutine nullifyCGNSDomPointers
