!
!      ******************************************************************
!      *                                                                *
!      * File:          initWallDistance.f90                            *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-27-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine initWallDistance(level, sps, allocMem)
!
!      ******************************************************************
!      *                                                                *
!      * initWallDistance allocates the memory for the wall distance,   *
!      * if needed, and initializes the wall distance to a large value. *
!      *                                                                *
!      ******************************************************************
!
       use block
       use constants
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, sps
       logical, intent(in)               :: allocMem
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn, il, jl, kl
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the domains.

       domain: do nn=1, nDom

         ! Allocate the memory for d2Wall, if desired.

         if( allocMem ) then

           il = flowDoms(nn,level,sps)%il
           jl = flowDoms(nn,level,sps)%jl
           kl = flowDoms(nn,level,sps)%kl

           allocate(flowDoms(nn,level,sps)%d2Wall(2:il,2:jl,2:kl), &
                    stat=ierr)
           if(ierr /= 0)                          &
             call terminate("initWallDistance", &
                            "Memory allocation failure for d2Wall")
         endif

         ! Initialize the wall distances to a large value.

         flowDoms(nn,level,sps)%d2Wall = large

       enddo domain

       end subroutine initWallDistance
