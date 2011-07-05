!
!      ******************************************************************
!      *                                                                *
!      * File:          cellRangeSubface.f90                            *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-07-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine cellRangeSubface
!
!      ******************************************************************
!      *                                                                *
!      * cellRangeSubface determines the cell range for every subface   *
!      * of every block all grid levels. This subrange can include one  *
!      * cell of overlap if the boundary coincides with the block       *
!      * boundary.                                                      *
!      *                                                                *
!      ******************************************************************
!
       use block
       use BCTypes
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nLevels, level
       integer(kind=intType) :: nn, mm, il, jl, kl, ie, je, ke
       integer(kind=intType) :: iBeg, jBeg, kBeg, iEnd, jEnd, kEnd
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number of grid levels.

       nLevels = ubound(flowDoms,2)

       ! Loop over the number of grid levels.

       levelLoop: do level=1,nLevels

         ! Loop over the blocks.

         domains: do nn=1,nDom

           ! Allocate the memory for the variables defining the cell
           ! range of the subfaces. Only allocated for the 1st spectral
           ! solution, because this info is identical for all of them.

           mm = flowDoms(nn,level,1)%nSubface
           allocate(flowDoms(nn,level,1)%icBeg(mm), &
                    flowDoms(nn,level,1)%jcBeg(mm), &
                    flowDoms(nn,level,1)%kcBeg(mm), &
                    flowDoms(nn,level,1)%icEnd(mm), &
                    flowDoms(nn,level,1)%jcEnd(mm), &
                    flowDoms(nn,level,1)%kcEnd(mm), stat=ierr)
           if(ierr /= 0)                        &
             call terminate("cellRangeSubface", &
                            "Memory allocation failure for &
                            &cell subranges")

           ! Store the nodal dimensions of the block a bit easier.

           il = flowDoms(nn,level,1)%il
           jl = flowDoms(nn,level,1)%jl
           kl = flowDoms(nn,level,1)%kl

           ie = flowDoms(nn,level,1)%ie
           je = flowDoms(nn,level,1)%je
           ke = flowDoms(nn,level,1)%ke

           ! Loop over the number of subfaces for this block.

           subfaces: do mm=1,flowDoms(nn,level,1)%nSubface

             ! Store the nodal range of the subface a bit easier.
             ! Make sure that iBeg, jBeg and kBeg contain the lowest and
             ! iEnd, jEnd and kEnd the highest node numbers.

             iBeg = min(flowDoms(nn,level,1)%inBeg(mm), &
                        flowDoms(nn,level,1)%inEnd(mm))
             iEnd = max(flowDoms(nn,level,1)%inBeg(mm), &
                        flowDoms(nn,level,1)%inEnd(mm))

             jBeg = min(flowDoms(nn,level,1)%jnBeg(mm), &
                        flowDoms(nn,level,1)%jnEnd(mm))
             jEnd = max(flowDoms(nn,level,1)%jnBeg(mm), &
                        flowDoms(nn,level,1)%jnEnd(mm))

             kBeg = min(flowDoms(nn,level,1)%knBeg(mm), &
                        flowDoms(nn,level,1)%knEnd(mm))
             kEnd = max(flowDoms(nn,level,1)%knBeg(mm), &
                        flowDoms(nn,level,1)%knEnd(mm))

             ! Determine the block face on which the subface is located
             ! and set the range accordingly.

             select case (flowDoms(nn,level,1)%BCFaceID(mm))

               case (iMin)
                 flowDoms(nn,level,1)%icBeg(mm) = 1
                 flowDoms(nn,level,1)%icEnd(mm) = 1

                 flowDoms(nn,level,1)%jcBeg(mm) = jBeg +1
                 if(jBeg == 1) flowDoms(nn,level,1)%jcBeg(mm) = 1

                 flowDoms(nn,level,1)%jcEnd(mm) = jEnd
                 if(jEnd == jl) flowDoms(nn,level,1)%jcEnd(mm) = je

                 flowDoms(nn,level,1)%kcBeg(mm) = kBeg +1
                 if(kBeg == 1) flowDoms(nn,level,1)%kcBeg(mm) = 1

                 flowDoms(nn,level,1)%kcEnd(mm) = kEnd
                 if(kEnd == kl) flowDoms(nn,level,1)%kcEnd(mm) = ke

               !=========================================================

               case (iMax)
                 flowDoms(nn,level,1)%icBeg(mm) = ie
                 flowDoms(nn,level,1)%icEnd(mm) = ie

                 flowDoms(nn,level,1)%jcBeg(mm) = jBeg +1
                 if(jBeg == 1) flowDoms(nn,level,1)%jcBeg(mm) = 1

                 flowDoms(nn,level,1)%jcEnd(mm) = jEnd
                 if(jEnd == jl) flowDoms(nn,level,1)%jcEnd(mm) = je

                 flowDoms(nn,level,1)%kcBeg(mm) = kBeg +1
                 if(kBeg == 1) flowDoms(nn,level,1)%kcBeg(mm) = 1

                 flowDoms(nn,level,1)%kcEnd(mm) = kEnd
                 if(kEnd == kl) flowDoms(nn,level,1)%kcEnd(mm) = ke

               !=========================================================

               case (jMin)
                 flowDoms(nn,level,1)%icBeg(mm) = iBeg +1
                 if(iBeg == 1) flowDoms(nn,level,1)%icBeg(mm) = 1

                 flowDoms(nn,level,1)%icEnd(mm) = iEnd
                 if(iEnd == il) flowDoms(nn,level,1)%icEnd(mm) = ie

                 flowDoms(nn,level,1)%jcBeg(mm) = 1
                 flowDoms(nn,level,1)%jcEnd(mm) = 1

                 flowDoms(nn,level,1)%kcBeg(mm) = kBeg +1
                 if(kBeg == 1) flowDoms(nn,level,1)%kcBeg(mm) = 1

                 flowDoms(nn,level,1)%kcEnd(mm) = kEnd
                 if(kEnd == kl) flowDoms(nn,level,1)%kcEnd(mm) = ke

               !=========================================================

               case (jMax)
                 flowDoms(nn,level,1)%icBeg(mm) = iBeg +1
                 if(iBeg == 1) flowDoms(nn,level,1)%icBeg(mm) = 1

                 flowDoms(nn,level,1)%icEnd(mm) = iEnd
                 if(iEnd == il) flowDoms(nn,level,1)%icEnd(mm) = ie

                 flowDoms(nn,level,1)%jcBeg(mm) = je
                 flowDoms(nn,level,1)%jcEnd(mm) = je

                 flowDoms(nn,level,1)%kcBeg(mm) = kBeg +1
                 if(kBeg == 1) flowDoms(nn,level,1)%kcBeg(mm) = 1

                 flowDoms(nn,level,1)%kcEnd(mm) = kEnd
                 if(kEnd == kl) flowDoms(nn,level,1)%kcEnd(mm) = ke

               !=========================================================

               case (kMin)
                 flowDoms(nn,level,1)%icBeg(mm) = iBeg +1
                 if(iBeg == 1) flowDoms(nn,level,1)%icBeg(mm) = 1

                 flowDoms(nn,level,1)%icEnd(mm) = iEnd
                 if(iEnd == il) flowDoms(nn,level,1)%icEnd(mm) = ie

                 flowDoms(nn,level,1)%jcBeg(mm) = jBeg +1
                 if(jBeg == 1) flowDoms(nn,level,1)%jcBeg(mm) = 1

                 flowDoms(nn,level,1)%jcEnd(mm) = jEnd
                 if(jEnd == jl) flowDoms(nn,level,1)%jcEnd(mm) = je

                 flowDoms(nn,level,1)%kcBeg(mm) = 1
                 flowDoms(nn,level,1)%kcEnd(mm) = 1

               !=========================================================

               case (kMax)
                 flowDoms(nn,level,1)%icBeg(mm) = iBeg +1
                 if(iBeg == 1) flowDoms(nn,level,1)%icBeg(mm) = 1

                 flowDoms(nn,level,1)%icEnd(mm) = iEnd
                 if(iEnd == il) flowDoms(nn,level,1)%icEnd(mm) = ie

                 flowDoms(nn,level,1)%jcBeg(mm) = jBeg +1
                 if(jBeg == 1) flowDoms(nn,level,1)%jcBeg(mm) = 1

                 flowDoms(nn,level,1)%jcEnd(mm) = jEnd
                 if(jEnd == jl) flowDoms(nn,level,1)%jcEnd(mm) = je

                 flowDoms(nn,level,1)%kcBeg(mm) = ke
                 flowDoms(nn,level,1)%kcEnd(mm) = ke

             end select

           enddo subfaces
         enddo domains
       enddo levelLoop

       end subroutine cellRangeSubface
