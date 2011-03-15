!
!      ******************************************************************
!      *                                                                *
!      * File:          cellCentroid.f90                                *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 08-15-2005                                      *
!      * Last modified: 08-28-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine cellCentroid(blk, level, sps, i, j, k, xc, yc, zc)
!
!      ******************************************************************
!      *                                                                *
!      * cellCentroid computes 8 times the centroid coordinates for the *
!      * cell given by (i,j,k) on the given block, level, and spectral  *
!      * mode. This is a simple calculation in most cases except when   *
!      * there is a possibility that the cell is a second level halo in *
!      * which case extrapolation is required and handled here.         *
!      *                                                                *
!      ******************************************************************
!
       use block
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: blk, level, sps, i, j, k

       real(kind=realType), intent(out) :: xc, yc, zc
!
!      Local variables.
!
       integer(kind=intType) :: ir, jr, kr
 
       real(kind=realType), dimension(2,2,2,3) :: xn
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Restrict the boundary cell indices to a cell where the
       ! nodal coordinates making up the cell exist (i.e. a halo with
       ! level of indirectness at most 1). Then store the nodes
       ! making up the cell in xn.

       ir = max(1_intType, min(i,flowDoms(blk,level,1)%ie))
       jr = max(1_intType, min(j,flowDoms(blk,level,1)%je))
       kr = max(1_intType, min(k,flowDoms(blk,level,1)%ke))

       xn = flowDoms(blk,level,sps)%x(ir-1:ir,jr-1:jr,kr-1:kr,:)

       ! Check if the boundary cell is a 2nd layer halo, and 
       ! extrapolate the coordinates if it is.

       if (i == flowDoms(blk,level,1)%ib) then
         xn(1,:,:,:) = two*xn(2,:,:,:) - xn(1,:,:,:)
       else if (i == 0) then
         xn(2,:,:,:) = two*xn(1,:,:,:) - xn(2,:,:,:)
       end if

       if (j == flowDoms(blk,level,1)%jb) then
         xn(:,1,:,:) = two*xn(:,2,:,:) - xn(:,1,:,:)
       else if (j == 0) then
         xn(:,2,:,:) = two*xn(:,1,:,:) - xn(:,2,:,:)
       end if

       if (k == flowDoms(blk,level,1)%kb) then
         xn(:,:,1,:) = two*xn(:,:,2,:) - xn(:,:,1,:)
       else if (k == 0) then
         xn(:,:,2,:) = two*xn(:,:,1,:) - xn(:,:,2,:)
       end if

       ! Compute eight times the centroid and then store the actual
       ! coordinates in the halo list.

       xc = sum(xn(:,:,:,1))
       yc = sum(xn(:,:,:,2))
       zc = sum(xn(:,:,:,3))

       end subroutine cellCentroid
