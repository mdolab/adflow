!
!      ******************************************************************
!      *                                                                *
!      * File:          update1to1Coarse.f90                            *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-08-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine update1to1Coarse(level, subface)
!
!      ******************************************************************
!      *                                                                *
!      * update1to1Coarse determines whether or not the given 1 to 1    *
!      * block boundary subface on the fine level is still a 1 to 1     *
!      * subface on the coarser grid level.                             *
!      *                                                                *
!      ******************************************************************
!
       use block
       use coarse1to1Subface
       use coarseningInfo
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType),        intent(in) :: level
       type(coarse1to1SubfaceType), intent(in) :: subface
!
!      Local variables.
!
       integer(kind=intType) :: levm1, nn, mm, ii, jj, kk, ll
       integer(kind=intType) :: i, j, k
       integer(kind=intType) :: iBeg, jBeg, kBeg, iEnd, jEnd, kEnd

       logical :: subfaceFound, idir, jdir, kdir
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Store the finer grid level in levm1

       levm1 = level -1

       ! Store the local block id a bit easier and store the minimum and
       ! maximum index in the three coordinate directions for the fine
       ! grid.

       nn = subface%neighBlock

       mm = subface%ndi
       iBeg = min(subface%idfine(1), subface%idfine(mm))
       iEnd = max(subface%idfine(1), subface%idfine(mm))

       mm = subface%ndj
       jBeg = min(subface%jdfine(1), subface%jdfine(mm))
       jEnd = max(subface%jdfine(1), subface%jdfine(mm))

       mm = subface%ndk
       kBeg = min(subface%kdfine(1), subface%kdfine(mm))
       kEnd = max(subface%kdfine(1), subface%kdfine(mm))

       ! Find the subface id of the block, which corresponds to the
       ! given subface of the subroutine header.

       subfaceFound = .false.
       findSubface: do ii=1,flowDoms(nn,levm1,1)%n1to1

         mm = ii + flowDoms(nn,levm1,1)%nBocos

         ! Check if this subface is the correct one.

         kk = min(flowDoms(nn,levm1,1)%inBeg(mm), &
                  flowDoms(nn,levm1,1)%inEnd(mm))
         ll = max(flowDoms(nn,levm1,1)%inBeg(mm), &
                  flowDoms(nn,levm1,1)%inEnd(mm))

         if(kk == iBeg .and. ll == iEnd) then
           kk = min(flowDoms(nn,levm1,1)%jnBeg(mm), &
                    flowDoms(nn,levm1,1)%jnEnd(mm))
           ll = max(flowDoms(nn,levm1,1)%jnBeg(mm), &
                    flowDoms(nn,levm1,1)%jnEnd(mm))

           if(kk == jBeg .and. ll == jEnd) then
             kk = min(flowDoms(nn,levm1,1)%knBeg(mm), &
                      flowDoms(nn,levm1,1)%knEnd(mm))
             ll = max(flowDoms(nn,levm1,1)%knBeg(mm), &
                      flowDoms(nn,levm1,1)%knEnd(mm))

             if(kk == kBeg .and. ll == kEnd) subfaceFound = .true.
           endif
         endif

         ! Exit the loop if this is indeed the subface.

         if( subfaceFound ) exit
       enddo findSubface

       ! The subface must be found on the fine grid. Check this.

       if(.not. subfaceFound)               &
         call terminate("update1to1Coarse", &
                        "Invalid fine grid 1 to 1 subface connectivity")

       ! Check the i-direction of the coarse grid subface to see if it
       ! is still a 1 to 1 subface. First check the number of grid lines.
       ! If these are identical, check each coarse grid line.

       idir = .true.
       ii = abs(flowDoms(nn,level,1)%inEnd(mm) &
          -     flowDoms(nn,level,1)%inBeg(mm)) + 1

       if(ii == subface%ndi) then
         do i=1,subface%ndi
           ii = subface%idfine(i)
           if(.not. flowDoms(nn,levm1,1)%ico(ii) ) idir = .false.
         enddo
       else
         idir = .false.
       endif

       ! Check the j-direction of the coarse grid subface to see if it
       ! is still a 1 to 1 subface. First check the number of grid lines.
       ! If these are identical, check each coarse grid line.

       jdir = .true.
       jj = abs(flowDoms(nn,level,1)%jnEnd(mm) &
          -     flowDoms(nn,level,1)%jnBeg(mm)) + 1

       if(jj == subface%ndj) then
         do j=1,subface%ndj
           jj = subface%jdfine(j)
           if(.not. flowDoms(nn,levm1,1)%jco(jj) ) jdir = .false.
         enddo
       else
         jdir = .false.
       endif

       ! Check the k-direction of the coarse grid subface to see if it
       ! is still a 1 to 1 subface. First check the number of grid lines.
       ! If these are identical, check each coarse grid line.

       kdir = .true.
       kk = abs(flowDoms(nn,level,1)%knEnd(mm) &
          -     flowDoms(nn,level,1)%knBeg(mm)) + 1

       if(kk == subface%ndk) then
         do k=1,subface%ndk
           kk = subface%kdfine(k)
           if(.not. flowDoms(nn,levm1,1)%kco(kk) ) kdir = .false.
         enddo
       else
         kdir = .false.
       endif

       ! Set coarseIs1to1 to .true. if the subface on the coarse grid
       ! is still a 1 to 1 subface. Otherwise set it to .false.

       ii = mm - flowDoms(nn,levm1,1)%nBocos
       if(idir .and. jdir .and. kdir) then
         coarseInfo(nn)%coarseIs1to1(ii) = .true.
       else
         coarseInfo(nn)%coarseIs1to1(ii) = .false.
       endif

       ! Store the donor range. It is possible that the lower and
       ! upper boundary must be reversed.

       if(flowDoms(nn,levm1,1)%dinBeg(mm) < &
          flowDoms(nn,levm1,1)%dinEnd(mm)) then
         flowDoms(nn,level,1)%dinBeg(mm) = subface%iBeg
         flowDoms(nn,level,1)%dinEnd(mm) = subface%iEnd
       else
         flowDoms(nn,level,1)%dinBeg(mm) = subface%iEnd
         flowDoms(nn,level,1)%dinEnd(mm) = subface%iBeg
       endif

       if(flowDoms(nn,levm1,1)%djnBeg(mm) < &
          flowDoms(nn,levm1,1)%djnEnd(mm)) then
         flowDoms(nn,level,1)%djnBeg(mm) = subface%jBeg
         flowDoms(nn,level,1)%djnEnd(mm) = subface%jEnd
       else
         flowDoms(nn,level,1)%djnBeg(mm) = subface%jEnd
         flowDoms(nn,level,1)%djnEnd(mm) = subface%jBeg
       endif

       if(flowDoms(nn,levm1,1)%dknBeg(mm) < &
          flowDoms(nn,levm1,1)%dknEnd(mm)) then
         flowDoms(nn,level,1)%dknBeg(mm) = subface%kBeg
         flowDoms(nn,level,1)%dknEnd(mm) = subface%kEnd
       else
         flowDoms(nn,level,1)%dknBeg(mm) = subface%kEnd
         flowDoms(nn,level,1)%dknEnd(mm) = subface%kBeg
       endif

       end subroutine update1to1Coarse
