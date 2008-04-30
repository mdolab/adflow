!
!      ******************************************************************
!      *                                                                *
!      * File:          extractBCStatesAdj.f90                             *
!      * Author:        C.A.(Sandy) Mader                               *
!      * Starting date: 04-16-2008                                      *
!      * Last modified: 04-16-2008                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine extractBCStatesAdj(nn,wAdj,pAdj,wAdj0, wAdj1, wAdj2,wAdj3,&
            pAdj0,pAdj1, pAdj2,pAdj3,&
            rlvAdj, revAdj,rlvAdj1, rlvAdj2,revAdj1, revAdj2,iOffset,&
            jOffset, kOffset,iCell, jCell,kCell,&
            isbeg,jsbeg,ksbeg,isend,jsend,ksend,ibbeg,jbbeg,kbbeg,ibend,&
            jbend,kbend,icbeg,jcbeg,icend,jcend,secondHalo)

!      ******************************************************************
!      *                                                                *
!      * copyBCStatesAdj copies the states needed for the boundary      *
!      * condition treatment on a general face, such that the boundary  *
!      * routines are only implemented once instead of 6 times.         *
!      * This routine takes the place of setBCPointers in the orginal   *
!      * Code                                                           *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers, only : ie, ib, je, jb, ke, kb, nBocos, &
                                  BCFaceID, BCType, BCData
       use flowVarRefState
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nn!, offset
       integer(kind=intType), intent(in) ::iCell, jCell,kCell
       integer(kind=intType), intent(in) ::isbeg,jsbeg,ksbeg,isend,jsend,ksend
       integer(kind=intType), intent(in) ::ibbeg,jbbeg,kbbeg,ibend,jbend,kbend
       integer(kind=intType), intent(out) ::icbeg,jcbeg,icend,jcend
       integer(kind=intType) ::irbeg,jrbeg,krbeg,irend,jrend,krend

       integer(kind=intType) :: iOffset, jOffset, kOffset

       real(kind=realType), dimension(-2:2,-2:2,-2:2,nw), &
                   intent(in) :: wAdj
       real(kind=realType), dimension(-2:2,-2:2,-2:2),intent(in) :: pAdj
       

       real(kind=realType), dimension(-2:2,-2:2,nw) :: wAdj0, wAdj1
       real(kind=realType), dimension(-2:2,-2:2,nw) :: wAdj2, wAdj3
       real(kind=realType), dimension(-2:2,-2:2)    :: pAdj0, pAdj1
       real(kind=realType), dimension(-2:2,-2:2)    :: pAdj2, pAdj3

       real(kind=realType), dimension(-2:2,-2:2,-2:2)::rlvAdj, revAdj
       real(kind=realType), dimension(-2:2,-2:2)::rlvAdj1, rlvAdj2
       real(kind=realType), dimension(-2:2,-2:2)::revAdj1, revAdj2
       
       logical :: iOverlap, jOverlap, kOverlap
       logical,intent(inout) :: secondHalo
       

!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
!!$       ! Determine the range of the stencil for the given cell.
!!$
!!$       iSBeg = iCell - 2; iSEnd = iCell + 2
!!$       jSBeg = jCell - 2; jSEnd = jCell + 2
!!$       kSBeg = kCell - 2; kSEnd = kCell + 2
!!$
!!$       iOffset = zero
!!$       jOffset = zero
!!$       kOffset = zero
!!$

!       ! Loop over the number of physical boundary subfaces of the block.
!
!       do nn=1,nBocos
!!$
!!$       ! Determine the range of halo cells which this boundary subface
!!$       ! will change.
!!$
!!$       !print *,'in extract states',icell,jcell,kcell,nn,isbeg,jsbeg,ksbeg
!!$       
!!$       select case (BCFaceID(nn))
!!$       case (iMin)
!!$          iBBeg = 0;                iBEnd = 1
!!$          jBBeg = BCData(nn)%icBeg; jBEnd = BCData(nn)%icEnd
!!$          kBBeg = BCData(nn)%jcBeg; kBEnd = BCData(nn)%jcEnd
!!$          
!!$          !=============================================================
!!$          
!!$       case (iMax)
!!$          iBBeg = ie;               iBEnd = ib
!!$          jBBeg = BCData(nn)%icBeg; jBEnd = BCData(nn)%icEnd
!!$          kBBeg = BCData(nn)%jcBeg; kBEnd = BCData(nn)%jcEnd
!!$          
!!$          !=============================================================
!!$          
!!$       case (jMin)
!!$          iBBeg = BCData(nn)%icBeg; iBEnd = BCData(nn)%icEnd
!!$          jBBeg = 0;                jBEnd = 1
!!$          kBBeg = BCData(nn)%jcBeg; kBEnd = BCData(nn)%jcEnd
!!$          
!!$          !=============================================================
!!$          
!!$       case (jMax)
!!$          iBBeg = BCData(nn)%icBeg; iBEnd = BCData(nn)%icEnd
!!$          jBBeg = je;               jBEnd = jb
!!$          kBBeg = BCData(nn)%jcBeg; kBEnd = BCData(nn)%jcEnd
!!$          
!!$          !=============================================================
!!$          
!!$       case (kMin)
!!$          iBBeg = BCData(nn)%icBeg; iBEnd = BCData(nn)%icEnd
!!$          jBBeg = BCData(nn)%jcBeg; jBEnd = BCData(nn)%jcEnd
!!$          kBBeg = 0;                kBEnd = 1
!!$          
!!$          !=============================================================
!!$          
!!$       case (kMax)
!!$          iBBeg = BCData(nn)%icBeg; iBEnd = BCData(nn)%icEnd
!!$          jBBeg = BCData(nn)%jcBeg; jBEnd = BCData(nn)%jcEnd
!!$          kBBeg = ke;               kBEnd = kb
!!$          
!!$       end select
!!$
!!$       !print *,'MID extract states',icell,jcell,kcell,nn,ibbeg,jbbeg,kbbeg,ibend,jbend,kbend,icbeg,jcbeg,icend,jcend
!!$       ! Check for an overlap between the stencil range and the
!!$       ! halo cells influenced by this boundary subface.
!!$       
!!$       iOverlap = .false.
!!$       if(iSBeg <= iBEnd .and. iSEnd >= iBBeg) iOverlap = .true.
!!$
!!$       jOverlap = .false.
!!$       if(jSBeg <= jBEnd .and. jSEnd >= jBBeg) jOverlap = .true.
!!$       
!!$       kOverlap = .false.
!!$       if(kSBeg <= kBEnd .and. kSEnd >= kBBeg) kOverlap = .true.
!!$
!!$       !print *,'overlap?',ioverlap,joverlap,koverlap

!       checkOverlap: if(iOverlap .and. jOverlap .and. kOverlap) then

          ! There is an overlap between the two ranges.
          ! Determine the actual overlap region.
          
          iRBeg = max(iSBeg, iBBeg); iREnd = min(iSEnd, iBEnd)
          jRBeg = max(jSBeg, jBBeg); jREnd = min(jSEnd, jBEnd)
          kRBeg = max(kSBeg, kBBeg); kREnd = min(kSEnd, kBEnd)
          
          ! Determine the offset values to be used in the computation
          ! of the halo cells. Either the cells with index -2 or with
          ! -1 corresponds to the starting cell ID of the subrange
          ! to be considered.
          
          if(iCell==iBBeg)then
             iOffset = iRBeg
          else
             iOffset = iRBeg + 1
          end if
          if(iSBeg == iRBeg) iOffset = iOffset + 1
          
          if(jCell==jBBeg)then
             jOffset = jRBeg
          else
             jOffset = jRBeg + 1
          end if
          if(jSBeg == jRBeg) jOffset = jOffset + 1
          
          if(kCell==kBBeg)then
             kOffset = kRBeg
          else
             kOffset = kRBeg + 1
          end if
          if(kSBeg == kRBeg) kOffset = kOffset + 1  
           

          ! Set some pointers depending on the situation.
          
          select case (BCFaceID(nn))
          case (iMin)

             ! At most the halo's of nearest neighbors enter the
             ! stencil. Make sure to limit properly.
             
             icBeg = max(jRBeg,jCell-1)
             icEnd = min(jREnd,jCell+1)
             jcBeg = max(kRBeg,kCell-1)
             jcEnd = min(kREnd,kCell+1)
             
             ! Other straightforward stuff.
             
             iOffset = jOffset
             jOffset = kOffset
             
             secondHalo = .true.
             if(iRBeg == iREnd) secondHalo = .false.
        
             if( secondHalo ) then
                wAdj0(-2:2,-2:2,1:nw) = wAdj(-2,-2:2,-2:2,1:nw)
                wAdj1(-2:2,-2:2,1:nw) = wAdj(-1,-2:2,-2:2,1:nw)
                wAdj2(-2:2,-2:2,1:nw) = wAdj( 0,-2:2,-2:2,1:nw)
                wAdj3(-2:2,-2:2,1:nw) = wAdj( 1,-2:2,-2:2,1:nw)
                
                pAdj0(-2:2,-2:2) = pAdj(-2,-2:2,-2:2)
                pAdj1(-2:2,-2:2) = pAdj(-1,-2:2,-2:2)
                pAdj2(-2:2,-2:2) = pAdj( 0,-2:2,-2:2)
                pAdj3(-2:2,-2:2) = pAdj( 1,-2:2,-2:2)
                
                if( viscous ) then
                   rlvAdj1(-2:2,-2:2) = rlvAdj(-1,-2:2,-2:2)
                   rlvAdj2(-2:2,-2:2) = rlvAdj( 0,-2:2,-2:2)
                endif
                if( eddyModel ) then
                   revAdj1(-2:2,-2:2) = revAdj(-1,-2:2,-2:2)
                   revAdj2(-2:2,-2:2) = revAdj( 0,-2:2,-2:2)
                endif
                
             else
                wAdj1(-2:2,-2:2,1:nw) = wAdj(-2,-2:2,-2:2,1:nw)
                wAdj2(-2:2,-2:2,1:nw) = wAdj(-1,-2:2,-2:2,1:nw)
                wAdj3(-2:2,-2:2,1:nw) = wAdj( 0,-2:2,-2:2,1:nw)
                
                pAdj1(-2:2,-2:2) = pAdj(-2,-2:2,-2:2)
                pAdj2(-2:2,-2:2) = pAdj(-1,-2:2,-2:2)
                pAdj3(-2:2,-2:2) = pAdj( 0,-2:2,-2:2)

                if( viscous ) then
                   rlvAdj1(-2:2,-2:2) = rlvAdj(-2,-2:2,-2:2)
                   rlvAdj2(-2:2,-2:2) = rlvAdj(-1,-2:2,-2:2)
                endif
                
                if( eddyModel ) then
                   revAdj1(-2:2,-2:2) = revAdj(-2,-2:2,-2:2)
                   revAdj2(-2:2,-2:2) = revAdj(-1,-2:2,-2:2)
                endif

             endif
             
             !===========================================================
             
          case (iMax)
             
             ! At most the halo's of nearest neighbors enter the
             ! stencil. Make sure to limit properly.
             
             icBeg = max(jRBeg,jCell-1)
             icEnd = min(jREnd,jCell+1)
             jcBeg = max(kRBeg,kCell-1)
             jcEnd = min(kREnd,kCell+1)
             
             ! Other straightforward stuff.
             
             iOffset = jOffset
             jOffset = kOffset
             
             secondHalo = .true.
             if(iRBeg == iREnd) secondHalo = .false.
             
             if( secondHalo ) then
                wAdj0(-2:2,-2:2,1:nw) = wAdj( 2,-2:2,-2:2,1:nw)
                wAdj1(-2:2,-2:2,1:nw) = wAdj( 1,-2:2,-2:2,1:nw)
                wAdj2(-2:2,-2:2,1:nw) = wAdj( 0,-2:2,-2:2,1:nw)
                wAdj3(-2:2,-2:2,1:nw) = wAdj(-1,-2:2,-2:2,1:nw)
                
                pAdj0(-2:2,-2:2) = pAdj( 2,-2:2,-2:2)
                pAdj1(-2:2,-2:2) = pAdj( 1,-2:2,-2:2)
                pAdj2(-2:2,-2:2) = pAdj( 0,-2:2,-2:2)
                pAdj3(-2:2,-2:2) = pAdj(-1,-2:2,-2:2)

                if( viscous ) then
                   rlvAdj1(-2:2,-2:2) = rlvAdj(1,-2:2,-2:2)
                   rlvAdj2(-2:2,-2:2) = rlvAdj(0,-2:2,-2:2)
                endif
                if( eddyModel ) then
                   revAdj1(-2:2,-2:2) = revAdj(1,-2:2,-2:2)
                   revAdj2(-2:2,-2:2) = revAdj(0,-2:2,-2:2)
                endif
             else
                wAdj1(-2:2,-2:2,1:nw) = wAdj( 2,-2:2,-2:2,1:nw)
                wAdj2(-2:2,-2:2,1:nw) = wAdj( 1,-2:2,-2:2,1:nw)
                wAdj3(-2:2,-2:2,1:nw) = wAdj( 0,-2:2,-2:2,1:nw)
                
                pAdj1(-2:2,-2:2) = pAdj( 2,-2:2,-2:2)
                pAdj2(-2:2,-2:2) = pAdj( 1,-2:2,-2:2)
                pAdj3(-2:2,-2:2) = pAdj( 0,-2:2,-2:2)

                if( viscous ) then
                   rlvAdj1(-2:2,-2:2) = rlvAdj(2,-2:2,-2:2)
                   rlvAdj2(-2:2,-2:2) = rlvAdj(1,-2:2,-2:2)
                endif

                if( eddyModel ) then
                   revAdj1(-2:2,-2:2) = revAdj(2,-2:2,-2:2)
                   revAdj2(-2:2,-2:2) = revAdj(1,-2:2,-2:2)
                endif
             endif
             
             !===========================================================
             
          case (jMin)
             
             ! At most the halo's of nearest neighbors enter the
             ! stencil. Make sure to limit properly.
             
             icBeg = max(iRBeg,iCell-1)
             icEnd = min(iREnd,iCell+1)
             jcBeg = max(kRBeg,kCell-1)
             jcEnd = min(kREnd,kCell+1)
             
             ! Other straightforward stuff.
             
             jOffset = kOffset
             
             secondHalo = .true.
             if(jRBeg == jREnd) secondHalo = .false.
             
             if( secondHalo ) then
                wAdj0(-2:2,-2:2,1:nw) = wAdj(-2:2,-2,-2:2,1:nw)
                wAdj1(-2:2,-2:2,1:nw) = wAdj(-2:2,-1,-2:2,1:nw)
                wAdj2(-2:2,-2:2,1:nw) = wAdj(-2:2, 0,-2:2,1:nw)
                wAdj3(-2:2,-2:2,1:nw) = wAdj(-2:2, 1,-2:2,1:nw)
                
                pAdj0(-2:2,-2:2) = pAdj(-2:2,-2,-2:2)
                pAdj1(-2:2,-2:2) = pAdj(-2:2,-1,-2:2)
                pAdj2(-2:2,-2:2) = pAdj(-2:2, 0,-2:2)
                pAdj3(-2:2,-2:2) = pAdj(-2:2, 1,-2:2)

                if( viscous ) then
                   rlvAdj1(-2:2,-2:2) = rlvAdj(-2:2,-1,-2:2)
                   rlvAdj2(-2:2,-2:2) = rlvAdj(-2:2,0,-2:2)
                endif
                if( eddyModel ) then
                   revAdj1(-2:2,-2:2) = revAdj(-2:2,-1,-2:2)
                   revAdj2(-2:2,-2:2) = revAdj(-2:2,0,-2:2)
                endif
             else
                wAdj1(-2:2,-2:2,1:nw) = wAdj(-2:2,-2,-2:2,1:nw)
                wAdj2(-2:2,-2:2,1:nw) = wAdj(-2:2,-1,-2:2,1:nw)
                wAdj3(-2:2,-2:2,1:nw) = wAdj(-2:2, 0,-2:2,1:nw)
                
                pAdj1(-2:2,-2:2) = pAdj(-2:2,-2,-2:2)
                pAdj2(-2:2,-2:2) = pAdj(-2:2,-1,-2:2)
                pAdj3(-2:2,-2:2) = pAdj(-2:2, 0,-2:2)
                
                if( viscous ) then
                   rlvAdj1(-2:2,-2:2) = rlvAdj(-2:2,-2,-2:2)
                   rlvAdj2(-2:2,-2:2) = rlvAdj(-2:2,-1,-2:2)
                endif
                if( eddyModel ) then
                   revAdj1(-2:2,-2:2) = revAdj(-2:2,-2,-2:2)
                   revAdj2(-2:2,-2:2) = revAdj(-2:2,-1,-2:2)
                endif

             endif
             
             !===========================================================
             
          case (jMax)
             
             ! At most the halo's of nearest neighbors enter the
             ! stencil. Make sure to limit properly.
             
             icBeg = max(iRBeg,iCell-1)
             icEnd = min(iREnd,iCell+1)
             jcBeg = max(kRBeg,kCell-1)
             jcEnd = min(kREnd,kCell+1)
             
             ! Other straightforward stuff.
             
             jOffset = kOffset
             
             secondHalo = .true.
             if(jRBeg == jREnd) secondHalo = .false.
             
             if( secondHalo ) then
                wAdj0(-2:2,-2:2,1:nw) = wAdj(-2:2, 2,-2:2,1:nw)
                wAdj1(-2:2,-2:2,1:nw) = wAdj(-2:2, 1,-2:2,1:nw)
                wAdj2(-2:2,-2:2,1:nw) = wAdj(-2:2, 0,-2:2,1:nw)
                wAdj3(-2:2,-2:2,1:nw) = wAdj(-2:2,-1,-2:2,1:nw)
                
                pAdj0(-2:2,-2:2) = pAdj(-2:2, 2,-2:2)
                pAdj1(-2:2,-2:2) = pAdj(-2:2, 1,-2:2)
                pAdj2(-2:2,-2:2) = pAdj(-2:2, 0,-2:2)
                pAdj3(-2:2,-2:2) = pAdj(-2:2,-1,-2:2)

                if( viscous ) then
                   rlvAdj1(-2:2,-2:2) = rlvAdj(-2:2,1,-2:2)
                   rlvAdj2(-2:2,-2:2) = rlvAdj(-2:2,0,-2:2)
                endif
                if( eddyModel ) then
                   revAdj1(-2:2,-2:2) = revAdj(-2:2,1,-2:2)
                   revAdj2(-2:2,-2:2) = revAdj(-2:2,0,-2:2)
                endif
             else
                wAdj1(-2:2,-2:2,1:nw) = wAdj(-2:2, 2,-2:2,1:nw)
                wAdj2(-2:2,-2:2,1:nw) = wAdj(-2:2, 1,-2:2,1:nw)
                wAdj3(-2:2,-2:2,1:nw) = wAdj(-2:2, 0,-2:2,1:nw)
                
                pAdj1(-2:2,-2:2) = pAdj(-2:2, 2,-2:2)
                pAdj2(-2:2,-2:2) = pAdj(-2:2, 1,-2:2)
                pAdj3(-2:2,-2:2) = pAdj(-2:2, 0,-2:2)

                if( viscous ) then
                   rlvAdj1(-2:2,-2:2) = rlvAdj(-2:2,2,-2:2)
                   rlvAdj2(-2:2,-2:2) = rlvAdj(-2:2,1,-2:2)
                endif
                if( eddyModel ) then
                   revAdj1(-2:2,-2:2) = revAdj(-2:2,2,-2:2)
                   revAdj2(-2:2,-2:2) = revAdj(-2:2,1,-2:2)
                endif
             endif
             
             !===========================================================
             
          case (kMin)
             
             ! At most the halo's of nearest neighbors enter the
             ! stencil. Make sure to limit properly.
             
             icBeg = max(iRBeg,iCell-1)
             icEnd = min(iREnd,iCell+1)
             jcBeg = max(jRBeg,jCell-1)
             jcEnd = min(jREnd,jCell+1)
             
             ! Other straightforward stuff.
             
             secondHalo = .true.
             if(kRBeg == kREnd) secondHalo = .false.
             
             if( secondHalo ) then
                wAdj0(-2:2,-2:2,1:nw) = wAdj(-2:2,-2:2,-2,1:nw)
                wAdj1(-2:2,-2:2,1:nw) = wAdj(-2:2,-2:2,-1,1:nw)
                wAdj2(-2:2,-2:2,1:nw) = wAdj(-2:2,-2:2, 0,1:nw)
                wAdj3(-2:2,-2:2,1:nw) = wAdj(-2:2,-2:2, 1,1:nw)
                
                pAdj0(-2:2,-2:2) = pAdj(-2:2,-2:2,-2)
                pAdj1(-2:2,-2:2) = pAdj(-2:2,-2:2,-1)
                pAdj2(-2:2,-2:2) = pAdj(-2:2,-2:2, 0)
                pAdj3(-2:2,-2:2) = pAdj(-2:2,-2:2, 1)

                if( viscous ) then
                   rlvAdj1(-2:2,-2:2) = rlvAdj(-2:2,-2:2,-1)
                   rlvAdj2(-2:2,-2:2) = rlvAdj(-2:2,-2:2,0)
                endif
                if( eddyModel ) then
                   revAdj1(-2:2,-2:2) = revAdj(-2:2,-2:2,-1)
                   revAdj2(-2:2,-2:2) = revAdj(-2:2,-2:2,0)
                endif
             else
                wAdj1(-2:2,-2:2,1:nw) = wAdj(-2:2,-2:2,-2,1:nw)
                wAdj2(-2:2,-2:2,1:nw) = wAdj(-2:2,-2:2,-1,1:nw)
                wAdj3(-2:2,-2:2,1:nw) = wAdj(-2:2,-2:2, 0,1:nw)
                
                pAdj1(-2:2,-2:2) = pAdj(-2:2,-2:2,-2)
                pAdj2(-2:2,-2:2) = pAdj(-2:2,-2:2,-1)
                pAdj3(-2:2,-2:2) = pAdj(-2:2,-2:2, 0)

                
                if( viscous ) then
                   rlvAdj1(-2:2,-2:2) = rlvAdj(-2:2,-2:2,-2)
                   rlvAdj2(-2:2,-2:2) = rlvAdj(-2:2,-2:2,-1)
                endif
                if( eddyModel ) then
                   revAdj1(-2:2,-2:2) = revAdj(-2:2,-2:2,-2)
                   revAdj2(-2:2,-2:2) = revAdj(-2:2,-2:2,-1)
                endif
             endif
             
             !===========================================================
             
          case (kMax)
             
             ! At most the halo's of nearest neighbors enter the
             ! stencil. Make sure to limit properly.
             
             icBeg = max(iRBeg,iCell-1)
             icEnd = min(iREnd,iCell+1)
             jcBeg = max(jRBeg,jCell-1)
             jcEnd = min(jREnd,jCell+1)
             
             ! Other straightforward stuff.
             
             secondHalo = .true.
             if(kRBeg == kREnd) secondHalo = .false.
             
             if( secondHalo ) then
                wAdj0(-2:2,-2:2,1:nw) = wAdj(-2:2,-2:2, 2,1:nw)
                wAdj1(-2:2,-2:2,1:nw) = wAdj(-2:2,-2:2, 1,1:nw)
                wAdj2(-2:2,-2:2,1:nw) = wAdj(-2:2,-2:2, 0,1:nw)
                wAdj3(-2:2,-2:2,1:nw) = wAdj(-2:2,-2:2,-1,1:nw)
                
                pAdj0(-2:2,-2:2) = pAdj(-2:2,-2:2, 2)
                pAdj1(-2:2,-2:2) = pAdj(-2:2,-2:2, 1)
                pAdj2(-2:2,-2:2) = pAdj(-2:2,-2:2, 0)
                pAdj3(-2:2,-2:2) = pAdj(-2:2,-2:2,-1)

                
                if( viscous ) then
                   rlvAdj1(-2:2,-2:2) = rlvAdj(-2:2,-2:2,1)
                   rlvAdj2(-2:2,-2:2) = rlvAdj(-2:2,-2:2,0)
                endif
                if( eddyModel ) then
                   revAdj1(-2:2,-2:2) = revAdj(-2:2,-2:2,1)
                   revAdj2(-2:2,-2:2) = revAdj(-2:2,-2:2,0)
                endif
             else
                wAdj1(-2:2,-2:2,1:nw) = wAdj(-2:2,-2:2, 2,1:nw)
                wAdj2(-2:2,-2:2,1:nw) = wAdj(-2:2,-2:2, 1,1:nw)
                wAdj3(-2:2,-2:2,1:nw) = wAdj(-2:2,-2:2, 0,1:nw)
                
                pAdj1(-2:2,-2:2) = pAdj(-2:2,-2:2, 2)
                pAdj2(-2:2,-2:2) = pAdj(-2:2,-2:2, 1)
                pAdj3(-2:2,-2:2) = pAdj(-2:2,-2:2, 0)

                
                if( viscous ) then
                   rlvAdj1(-2:2,-2:2) = rlvAdj(-2:2,-2:2,2)
                   rlvAdj2(-2:2,-2:2) = rlvAdj(-2:2,-2:2,1)
                endif
                if( eddyModel ) then
                   revAdj1(-2:2,-2:2) = revAdj(-2:2,-2:2,2)
                   revAdj2(-2:2,-2:2) = revAdj(-2:2,-2:2,1)
                endif
             endif
             
          end select
          
          
   !    end if checkOverlap
 
  !  enddo
    
     end subroutine extractBCStatesAdj
  
