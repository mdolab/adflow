!
!      ******************************************************************
!      *                                                                *
!      * File:          changeIblanks.f90                               *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 02-23-2005                                      *
!      * Last modified: 10-16-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine changeIblanks(setHalos, bndryOption)
!
!      ******************************************************************
!      *                                                                *
!      * ChangeIblanks is a utility routine for the iblank array. The   *
!      * input arguments control the following functions:               *
!      *                                                                *
!      * setHalos:     If true, the halos of all blocks on the given    *
!      *               level are assigned a value based on their        *
!      *               boundary condition type. The value assigned is 2 *
!      *               for all bocos types except OversetOuterBound     *
!                      which gets -1.                                   *
!      * bndryOption:  Input 0, 1, or 9. 0 will set the iblank of the   *
!      *               overset boundary cells to 0 as needed by the     *
!      *               solver. 1 will set them to a multiple of 10      *
!      *               where the multiplier is the index into the array *
!      *               ibndry, idonor, etc. 10 sets only the values for *
!      *               the boundary cells with donors.                  *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       implicit none
!
!      Subroutine arguments
!
       integer(kind=intType), intent(in) :: bndryOption
       logical,               intent(in) :: setHalos
!
!      Local variables.
!
       integer(kind=intType) :: m, ii
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set the halos if needed.
 
       halos: if (setHalos) then

         ! Intialize all the halos to -1. This is the value reserved for
         ! bcType = oversetOuterBound. Then loop over the subfaces and
         ! allow all other bc's to overwrite the corner halos. This way
         ! the priority is such that 1-to-1 faces have highest priority
         ! to corner halos, all other bocos have lesser, and overset
         ! outer boundaries have the least priority.

         iblank(0:1, : , : ) = -1
         iblank( : ,0:1, : ) = -1
         iblank( : , : ,0:1) = -1

         iblank(ie:ib,  :  ,  :  ) = -1
         iblank(  :  ,je:jb,  :  ) = -1
         iblank(  :  ,  :  ,ke:kb) = -1

         ! Loop over the number of subfaces for this block.

         subfaces: do m = 1,nSubface
 
           ! Skip the subface if it's an overset outer boundary.

           if (BCType(m) == OversetOuterBound) cycle
 
           ! Set the iblank.
 
           iblank(icBeg(m):icEnd(m), &
                  jcBeg(m):jcEnd(m), &
                  kcBeg(m):kcEnd(m)) = 2

         end do subfaces

         ! Extend the iblank values to the 2nd level halos.

         ! I faces.

         iblank( 0,1:je,1:ke) = iblank( 1,1:je,1:ke)
         iblank(ib,1:je,1:ke) = iblank(ie,1:je,1:ke)

         ! J faces.

         iblank(0:ib, 0,1:ke) = iblank(0:ib, 1,1:ke)
         iblank(0:ib,jb,1:ke) = iblank(0:ib,je,1:ke)

         ! K faces.

         iblank(0:ib,0:jb, 0) = iblank(0:ib,0:jb, 1)
         iblank(0:ib,0:jb,kb) = iblank(0:ib,0:jb,ke)

       end if halos

       ! Set the boundary cells based on the input option.

       select case (bndryOption)

         case (0_intType)

           ! This is the solver option, set to 0 so residuals are
           ! zeroed on the boundary.

           do m = 1,nCellsOverset + nOrphans
             iblank(ibndry(1,m),ibndry(2,m),ibndry(3,m)) = 0
           end do

         case (1_intType)

           ! Multiple of 10 for easy access to array position via
           ! iblank value. Orphans are set to 9.

           do m = 1,nCellsOverset
             iblank(ibndry(1,m),ibndry(2,m),ibndry(3,m)) = m*10
           end do

           do m = 1,nOrphans
             ii = nCellsOverset + m
             iblank(ibndry(1,ii),ibndry(2,ii),ibndry(3,ii)) = 9
           end do

         case (10_intType)

           ! Multiple of 10 for easy access to array position via
           ! iblank value. Ignore the orphans for this option.

           do m = 1,nCellsOverset
             iblank(ibndry(1,m),ibndry(2,m),ibndry(3,m)) = m*10
           end do

       end select

       end subroutine changeIblanks
