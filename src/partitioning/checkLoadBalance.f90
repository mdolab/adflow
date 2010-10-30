!
!      ******************************************************************
!      *                                                                *
!      * File:          checkLoadBalance.f90                            *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-11-2003                                      *
!      * Last modified: 11-22-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine checkLoadBalance(cellsBalanced, facesBalanced)
!
!      ******************************************************************
!      *                                                                *
!      * checkLoadBalance determines whether or not the load balance    *
!      * for the cells and faces is met.                                *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use constants
       use partitionMod
       implicit none
!
!      Subroutine arguments.
!
       logical, intent(out) :: cellsBalanced, facesBalanced
!
!      Local variables
!
       integer(kind=intType) :: i, j
       integer(kind=intType) :: nCellMax, nCellTol
       integer(kind=intType) :: nFaceMax, nFaceTol

       integer(kind=intType), dimension(nProc) :: nCell, nFace

       integer(kind=8) :: nCellsEven, nFacesEven   ! 8 byte integers to
                                                   ! avoid overflow.
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize nCell and nFace to 0. These variables will contain
       ! the number of cells and faces per partition (== processor)
       ! respectively.

       nCell = 0
       nFace = 0

       ! Determine the number of cells and faces per partition.
       ! Note that part(i) is the processor id, which starts at 0.

       do i=1,nblocks
         j = part(i) + 1
         nCell(j) = nCell(j) + blocks(i)%nCell
         nFace(j) = nFace(j) + blocks(i)%nFace
       enddo

       ! Determine the desirable number of cells and faces per processor.

       nCellsEven = nCell(1)
       nFacesEven = nFace(1)

       do i=2,nProc
         nCellsEven = nCellsEven + nCell(i)
         nFacesEven = nFacesEven + nFace(i)
       enddo

       nCellsEven = nCellsEven/nProc
       nFacesEven = nFacesEven/nProc

       ! Determine the maximum value of nCell and nFace
       ! and substract the optimal value.

       nCellMax = abs(maxval(nCell) - nCellsEven)
       nFaceMax = abs(maxval(nFace) - nFacesEven)

       ! Determine the tolerance for the cells and faces.

       nCellTol = (ubvec(1) - one)*nCellsEven
       nFaceTol = (ubvec(2) - one)*nFacesEven

       ! Check whether the load balance values for the cells and faces
       ! are met.

       cellsBalanced = .true.
       facesBalanced = .true.

       if(nCellMax > nCellTol) cellsBalanced = .false.
       if(nFaceMax > nFaceTol) facesBalanced = .false.

       ! Determine the load imbalances for the cells and faces
       ! and store it in ubvec.

       ubvec(1) = real(nCellMax,realType) &
                / real(nCellsEven,realType)
       ubvec(2) = real(nFaceMax,realType) &
                / real(nFacesEven,realType)

       end subroutine checkLoadBalance
