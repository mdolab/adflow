!
!      ******************************************************************
!      *                                                                *
!      * File:          estimateCoarseDonor.F90                         *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 04-10-2005                                      *
!      * Last modified: 10-16-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine estimateCoarseDonor(b,l,t, i,j,k, xc,yc,zc, &
                                      donorProc,donorBlock,di,dj,dk)
!
!      ******************************************************************
!      *                                                                *
!      * estimateCoarseDonor initializes the donor processor, block,    *
!      * and indices for a coarse boundary cell based on the donor info *
!      * of the next finer level. The coarse cell on block b, level l,  *
!      * spectral solution t, and with indices i,j,k is searched for    *
!      * finer grid boundary cells and the donor info for the one       *
!      * closest to xc,yc,zc is copied for output. Note that the input  *
!      * coordinates should be 8 times the actual. If the finer cell    *
!      * does not contain any boundary cells, the donor estimate cannot *
!      * be determined and -1 is returned in donorProc.                 *
!      *                                                                *
!      ******************************************************************
!
       use block
       implicit none
!
!      Subroutine arguments
!
       integer(kind=intType), intent(in)  :: b, l, t, i, j, k
       integer(kind=intType), intent(out) :: donorProc, donorBlock
       integer(kind=intType), intent(out) :: di, dj, dk

       real(kind=realType),   intent(in)  :: xc, yc, zc
!
!      Local variables.
!
       integer(kind=intType) :: ss, fl, ii, jj, kk
 
       real(kind=realType) :: xs, ys, zs, dist, minDist
 
       logical :: foundOne
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Next finer level.
 
       fl = l - 1

       ! Set donor process to -1 and the block and indices to 0 in case 
       ! donor info cannot be found.

       donorProc  = -1
       donorBlock = 0
       di = 0; dj = 0; dk = 0
 
       ! Loop over the fine grid cells inside the coarse cell given by
       ! i,j,k and copy the donor information of the fine grid cell whose
       ! centroid is closest to that of the target coarse cell xc,yc,zc.
 
       foundOne = .false.
       minDist  = zero
 
       do kk=flowDoms(b,l,1)%mgKFine(k,1),flowDoms(b,l,1)%mgKFine(k,2)
         do jj=flowDoms(b,l,1)%mgJFine(j,1),flowDoms(b,l,1)%mgJFine(j,2)
           do ii=flowDoms(b,l,1)%mgIFine(i,1),flowDoms(b,l,1)%mgIFine(i,2)
 
             ! Test if this is a fine grid boundary cell.
 
             fineBndry: if (flowDoms(b,fl,t)%iblank(ii,jj,kk) >= 10) then
 
               ! Calculate 8 times the centroid and the distance
               ! to the coarse cell centroid.

               call cellCentroid(b, fl, t, ii, jj, kk, xs, ys, zs) 
 
               dist = (xc - xs)**2 + (yc - ys)**2 + (zc - zs)**2
 
               if (foundOne .and. dist >= minDist) cycle
 
               ! This is a better cell. Find the index into the
               ! finer grid boundary info and copy to outputs.
 
               foundOne = .true.
               minDist  = dist
 
               ss = flowDoms(b,fl,t)%iblank(ii,jj,kk)/10
 
               donorProc  = flowDoms(b,fl,t)%neighProcOver(ss)
               donorBlock = flowDoms(b,fl,t)%neighBlockOver(ss)
 
               di = flowDoms(b,fl,t)%idonor(1,ss)
               dj = flowDoms(b,fl,t)%idonor(2,ss)
               dk = flowDoms(b,fl,t)%idonor(3,ss)
 
             end if fineBndry

           end do
         end do
       end do

       end subroutine estimateCoarseDonor
