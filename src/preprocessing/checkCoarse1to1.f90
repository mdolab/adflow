!
!      ******************************************************************
!      *                                                                *
!      * File:          checkCoarse1to1.f90                             *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-09-2003                                      *
!      * Last modified: 10-26-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine checkCoarse1to1(level)
!
!      ******************************************************************
!      *                                                                *
!      * checkCoarse1to1 removes the nonmatching block boundaries       *
!      * from the list of 1 to 1 matching ones. They are in there,      *
!      * because they are 1 to 1 matching on the finer grids.           *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use block
       use inputTimeSpectral
       use coarseningInfo
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level
!
!      Local variables.
!
       character(len=maxStringLen) :: errorMessage

       integer(kind=intType) :: i, nn, mm, kk, ll
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of domains.

       domains: do nn=1,nDom

         ! Loop over the number of 1 to 1 subfaces. As this number can
         ! change during the loop, the control statement is done via
         ! an exit.

         i  = 1
         n1to1: do

           ! Exit the loop if the counter is larger than the number
           ! of 1 to 1 subfaces.

           if(i > flowDoms(nn,level,1)%n1to1) exit

           ! Add the offset to i to store the correct place in the
           ! subface info.

           mm = i + flowDoms(nn,level,1)%nBocos

           ! Test if this is still a 1 to 1 subface on the coarse grid.

           is1to1: if( coarseInfo(nn)%coarseIs1to1(i) ) then

             ! Subface is still a 1 to 1 subface on the coarse grid.
             ! Only the counter i must be updated.

             i = i + 1

           else is1to1

             ! Due to the coarsening the subface is not a 1 to 1 block
             ! boundary anymore. Swap the current entry with the entry
             ! kk.

             kk = flowDoms(nn,level,1)%nBocos + flowDoms(nn,level,1)%n1to1

             ! The range info.

             ll                             = flowDoms(nn,level,1)%inBeg(mm)
             flowDoms(nn,level,1)%inBeg(mm) = flowDoms(nn,level,1)%inBeg(kk)
             flowDoms(nn,level,1)%inBeg(kk) = ll

             ll                             = flowDoms(nn,level,1)%jnBeg(mm)
             flowDoms(nn,level,1)%jnBeg(mm) = flowDoms(nn,level,1)%jnBeg(kk)
             flowDoms(nn,level,1)%jnBeg(kk) = ll

             ll                             = flowDoms(nn,level,1)%knBeg(mm)
             flowDoms(nn,level,1)%knBeg(mm) = flowDoms(nn,level,1)%knBeg(kk)
             flowDoms(nn,level,1)%knBeg(kk) = ll

             ll                             = flowDoms(nn,level,1)%inEnd(mm)
             flowDoms(nn,level,1)%inEnd(mm) = flowDoms(nn,level,1)%inEnd(kk)
             flowDoms(nn,level,1)%inEnd(kk) = ll

             ll                             = flowDoms(nn,level,1)%jnEnd(mm)
             flowDoms(nn,level,1)%jnEnd(mm) = flowDoms(nn,level,1)%jnEnd(kk)
             flowDoms(nn,level,1)%jnEnd(kk) = ll

             ll                             = flowDoms(nn,level,1)%knEnd(mm)
             flowDoms(nn,level,1)%knEnd(mm) = flowDoms(nn,level,1)%knEnd(kk)
             flowDoms(nn,level,1)%knEnd(kk) = ll

             ! The donor info.

             ll                              = flowDoms(nn,level,1)%dinBeg(mm)
             flowDoms(nn,level,1)%dinBeg(mm) = flowDoms(nn,level,1)%dinBeg(kk)
             flowDoms(nn,level,1)%dinBeg(kk) = ll

             ll                              = flowDoms(nn,level,1)%djnBeg(mm)
             flowDoms(nn,level,1)%djnBeg(mm) = flowDoms(nn,level,1)%djnBeg(kk)
             flowDoms(nn,level,1)%djnBeg(kk) = ll

             ll                              = flowDoms(nn,level,1)%dknBeg(mm)
             flowDoms(nn,level,1)%dknBeg(mm) = flowDoms(nn,level,1)%dknBeg(kk)
             flowDoms(nn,level,1)%dknBeg(kk) = ll

             ll                              = flowDoms(nn,level,1)%dinEnd(mm)
             flowDoms(nn,level,1)%dinEnd(mm) = flowDoms(nn,level,1)%dinEnd(kk)
             flowDoms(nn,level,1)%dinEnd(kk) = ll

             ll                              = flowDoms(nn,level,1)%djnEnd(mm)
             flowDoms(nn,level,1)%djnEnd(mm) = flowDoms(nn,level,1)%djnEnd(kk)
             flowDoms(nn,level,1)%djnEnd(kk) = ll

             ll                              = flowDoms(nn,level,1)%dknEnd(mm)
             flowDoms(nn,level,1)%dknEnd(mm) = flowDoms(nn,level,1)%dknEnd(kk)
             flowDoms(nn,level,1)%dknEnd(kk) = ll

             ! The transformation matrix.

             ll                          = flowDoms(nn,level,1)%l1(mm)
             flowDoms(nn,level,1)%l1(mm) = flowDoms(nn,level,1)%l1(kk)
             flowDoms(nn,level,1)%l1(kk) = ll

             ll                          = flowDoms(nn,level,1)%l2(mm)
             flowDoms(nn,level,1)%l2(mm) = flowDoms(nn,level,1)%l2(kk)
             flowDoms(nn,level,1)%l2(kk) = ll

             ll                          = flowDoms(nn,level,1)%l3(mm)
             flowDoms(nn,level,1)%l3(mm) = flowDoms(nn,level,1)%l3(kk)
             flowDoms(nn,level,1)%l3(kk) = ll

             ! The rest of the subface info.

             ll = flowDoms(nn,level,1)%BCType(mm)
             flowDoms(nn,level,1)%BCType(mm) = &
                                 flowDoms(nn,level,1)%BCType(kk)
             flowDoms(nn,level,1)%BCType(kk) = ll

             ll = flowDoms(nn,level,1)%BCFaceID(mm)
             flowDoms(nn,level,1)%BCFaceID(mm) = &
                                 flowDoms(nn,level,1)%BCFaceID(kk)
             flowDoms(nn,level,1)%BCFaceID(kk) = ll

             ll = flowDoms(nn,level,1)%neighBlock(mm)
             flowDoms(nn,level,1)%neighBlock(mm) = &
                                  flowDoms(nn,level,1)%neighBlock(kk)
             flowDoms(nn,level,1)%neighBlock(kk) = ll

             ll = flowDoms(nn,level,1)%neighProc(mm)
             flowDoms(nn,level,1)%neighProc(mm) = &
                                 flowDoms(nn,level,1)%neighProc(kk)
             flowDoms(nn,level,1)%neighProc(kk) = ll

             ll = flowDoms(nn,level,1)%groupNum(mm)
             flowDoms(nn,level,1)%groupNum(mm) = &
                                 flowDoms(nn,level,1)%groupNum(kk)
             flowDoms(nn,level,1)%groupNum(kk) = ll

             ll = flowDoms(nn,level,1)%cgnsSubface(mm)
             flowDoms(nn,level,1)%cgnsSubface(mm) = &
                                 flowDoms(nn,level,1)%cgnsSubface(kk)
             flowDoms(nn,level,1)%cgnsSubface(kk) = ll

             ! Decrease the number of 1 to 1 block boundaries. Note the
             ! counter i should not be updated.

             flowDoms(nn,level,1)%n1to1 = flowDoms(nn,level,1)%n1to1 - 1

 101         format("Non-matching block-to-block face on zone ",a, &
                    "...Support not implemented yet.")
             ll = flowDoms(nn,1,1)%cgnsBlockID
             write(errorMessage,101) cgnsDoms(ll)%zoneName
             call terminate("checkCoarse1to1", errorMessage)

           endif is1to1

         enddo n1to1

         ! Copy the number of internal subfaces for the rest of
         ! the spectral solutions.

         do i=2,nTimeIntervalsSpectral
           flowDoms(nn,level,i)%n1to1 = flowDoms(nn,level,1)%n1to1
         enddo

       enddo domains

       end subroutine checkCoarse1to1
