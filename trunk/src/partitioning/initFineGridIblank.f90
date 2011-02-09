!
!      ******************************************************************
!      *                                                                *
!      * File:          initFineGridIblank.f90                          *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 03-27-2005                                      *
!      * Last modified: 07-22-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine initFineGridIblank
!
!      ******************************************************************
!      *                                                                *
!      * InitFineGridIblank allocates the fine grid iblank array and    *
!      * initializes the values for the holes, boundary, and halos. The *
!      * holes read into the cgns domains are distributed amongst its   *
!      * sublocks in the form of iblanks. That is, we do not store a    *
!      * list of indices for the holes of the flow domains as done in   *
!      * the CGNS. The number of holes in each domain are also counted. *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use block
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: i, j, k, l, m, n, cgnsId
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the local blocks.

       domains: do n = 1,nDom

         ! Allocate memory for the iblank array of this block.

         i = flowDoms(n,1,1)%ib
         j = flowDoms(n,1,1)%jb
         k = flowDoms(n,1,1)%kb

         allocate(flowDoms(n,1,1)%iblank(0:i,0:j,0:k), &
                  stat=ierr)
         if(ierr /= 0)                             &
           call terminate("initFineGridIblank", &
                          "Memory allocation failure for iblank")
 
         ! Initialize iblank to 1 everywhere, and the number of holes
         ! for this domain to 0.
 
         flowDoms(n,1,1)%iblank = 1
         flowDoms(n,1,1)%nHoles = 0
 
         ! Loop over the holes of the cgns block that this domain belongs
         ! to and the indices of that hole.
 
         cgnsID = flowDoms(n,1,1)%cgnsBlockId
 
         holes: do m = 1,cgnsDoms(cgnsId)%nholes
           indices: do l = 1,cgnsDoms(cgnsId)%hole(m)%npnts
 
             ! Store the indices a little easier.
 
             i = cgnsDoms(cgnsId)%hole(m)%indices(1,l)
             j = cgnsDoms(cgnsId)%hole(m)%indices(2,l)
             k = cgnsDoms(cgnsId)%hole(m)%indices(3,l)
 
             ! Check if these cell indices belong to this domain.
 
             belongToMe: if (i <  flowDoms(n,1,1)%iEndor .and. &
                             i >= flowDoms(n,1,1)%iBegor .and. &
                             j <  flowDoms(n,1,1)%jEndor .and. &
                             j >= flowDoms(n,1,1)%jBegor .and. &
                             k <  flowDoms(n,1,1)%kEndor .and. &
                             k >= flowDoms(n,1,1)%kBegor) then
 
               ! Convert the indices from the original domain to
               ! this domain due to a possible split. Also, note
               ! the additional +1 offset because owned cells
               ! start at index 2.
 
               i = i - flowDoms(n,1,1)%iBegor + 2
               j = j - flowDoms(n,1,1)%jBegor + 2
               k = k - flowDoms(n,1,1)%kBegor + 2
 
               ! Set the iblank to 0 and update the number of holes 
               ! in this domain.
 
               flowDoms(n,1,1)%iblank(i,j,k) = 0
               flowDoms(n,1,1)%nHoles = flowDoms(n,1,1)%nHoles + 1

             end if belongToMe
           end do indices
         end do holes

         ! Now loop over the boundary for this domain and subtract for
         ! any boundary cells that were also input as holes. This is just
         ! a robustness thing because the CGNS is somewhat ambiguous in
         ! saying what is a hole. In this code, boundary cells are not
         ! holes, although their residuals are blanked.

         boundary: do m = 1,flowDoms(n,1,1)%nCellsOverset

           i = flowDoms(n,1,1)%ibndry(1,m)
           j = flowDoms(n,1,1)%ibndry(2,m)
           k = flowDoms(n,1,1)%ibndry(3,m)

           if (flowDoms(n,1,1)%iblank(i,j,k) == 0) then
             flowDoms(n,1,1)%nHoles = flowDoms(n,1,1)%nHoles - 1
           end if

         end do boundary
       end do domains
 
       end subroutine initFineGridIblank
