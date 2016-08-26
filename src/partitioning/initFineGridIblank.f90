       subroutine initFineGridIblank
!
!       InitFineGridIblank allocates the fine grid iblank array and    
!       initializes the values for the holes, boundary, and halos. The 
!       holes read into the cgns domains are distributed amongst its   
!       sublocks in the form of iblanks. That is, we do not store a    
!       list of indices for the holes of the flow domains as done in   
!       the CGNS. The number of holes in each domain are also counted. 
!
       use constants
       use block, only : nDom, flowDoms
       use utils, only : terminate
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: i, j, k, l, m, n, cgnsId

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
 
      end do domains
 
       end subroutine initFineGridIblank
