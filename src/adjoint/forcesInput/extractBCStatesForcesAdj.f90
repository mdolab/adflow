!
!      ******************************************************************
!      *                                                                *
!      * File:          extractBCStatesForcesAdj.f90                    *
!      * Author:        C.A.(Sandy) Mader                               *
!      * Starting date: 04-16-2008                                      *
!      * Last modified: 06-09-2008                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine extractBCStatesForcesAdj(nn,wAdj,pAdj,wAdj0, wAdj1,&
            wAdj2,wAdj3,pAdj0,pAdj1, pAdj2,pAdj3,&
            rlvAdj, revAdj,rlvAdj1, rlvAdj2,revAdj1, revAdj2,&
            icbeg,jcbeg,icend,jcend,secondHalo)

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
       integer(kind=intType) ::icbeg,jcbeg,icend,jcend


       real(kind=realType), dimension(0:ib,0:jb,0:kb,nw), intent(in) :: wAdj
       real(kind=realType), dimension(0:ib,0:jb,0:kb), intent(in) :: pAdj
       logical, intent(in) :: secondHalo

       real(kind=realType), dimension(icbeg:icend,jcbeg:jcend,nw):: wAdj0, wAdj1
       real(kind=realType), dimension(icbeg:icend,jcbeg:jcend,nw):: wAdj2, wAdj3
       real(kind=realType), dimension(icbeg:icend,jcbeg:jcend) :: pAdj0, pAdj1
       real(kind=realType), dimension(icbeg:icend,jcbeg:jcend) :: pAdj2, pAdj3

       real(kind=realType), dimension(0:ib,0:jb,0:kb)::rlvAdj, revAdj
       real(kind=realType), dimension(icbeg:icend,jcbeg:jcend)::rlvAdj1, rlvAdj2
       real(kind=realType), dimension(icbeg:icend,jcbeg:jcend)::revAdj1, revAdj2
       
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!

       ! Set some pointers depending on the situation.
       !print *,'selecting case,bcfacid',BCFaceID(nn)
       select case (BCFaceID(nn))
       case (iMin)
          
          ! At most the halo's of nearest neighbors enter the
          ! stencil. Make sure to limit properly.
          
          !icBeg = max(jRBeg,jCell-1)
          !icEnd = min(jREnd,jCell+1)
          !jcBeg = max(kRBeg,kCell-1)
          !jcEnd = min(kREnd,kCell+1)
             
          ! Other straightforward stuff.
          
          !iOffset = jOffset
          !jOffset = kOffset
             
          !secondHalo = .true.
          !if(iRBeg == iREnd) secondHalo = .false.
        
          if( secondHalo ) then
             wAdj0(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(0,icbeg:icend,jcbeg:jcend,1:nw)
             wAdj1(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(0,icbeg:icend,jcbeg:jcend,1:nw)
             wAdj2(icbeg:icend,jcbeg:jcend,1:nw) = wAdj( 2,icbeg:icend,jcbeg:jcend,1:nw)
             wAdj3(icbeg:icend,jcbeg:jcend,1:nw) = wAdj( 3,icbeg:icend,jcbeg:jcend,1:nw)
             
             pAdj0(icbeg:icend,jcbeg:jcend) = pAdj(0,icbeg:icend,jcbeg:jcend)
             pAdj1(icbeg:icend,jcbeg:jcend) = pAdj(1,icbeg:icend,jcbeg:jcend)
             pAdj2(icbeg:icend,jcbeg:jcend) = pAdj(2,icbeg:icend,jcbeg:jcend)
             pAdj3(icbeg:icend,jcbeg:jcend) = pAdj(3,icbeg:icend,jcbeg:jcend)
             
             if( viscous ) then
                rlvAdj1(icbeg:icend,jcbeg:jcend) = rlvAdj(1,icbeg:icend,jcbeg:jcend)
                rlvAdj2(icbeg:icend,jcbeg:jcend) = rlvAdj(2,icbeg:icend,jcbeg:jcend)
             endif
             if( eddyModel ) then
                revAdj1(icbeg:icend,jcbeg:jcend) = revAdj(1,icbeg:icend,jcbeg:jcend)
                revAdj2(icbeg:icend,jcbeg:jcend) = revAdj(2,icbeg:icend,jcbeg:jcend)
             endif
             
          else
                wAdj1(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(0,icbeg:icend,jcbeg:jcend,1:nw)
                wAdj2(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(1,icbeg:icend,jcbeg:jcend,1:nw)
                wAdj3(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(2,icbeg:icend,jcbeg:jcend,1:nw)
                
                pAdj1(icbeg:icend,jcbeg:jcend) = pAdj(0,icbeg:icend,jcbeg:jcend)
                pAdj2(icbeg:icend,jcbeg:jcend) = pAdj(1,icbeg:icend,jcbeg:jcend)
                pAdj3(icbeg:icend,jcbeg:jcend) = pAdj(2,icbeg:icend,jcbeg:jcend)

                if( viscous ) then
                   rlvAdj1(icbeg:icend,jcbeg:jcend) = rlvAdj(0,icbeg:icend,jcbeg:jcend)
                   rlvAdj2(icbeg:icend,jcbeg:jcend) = rlvAdj(1,icbeg:icend,jcbeg:jcend)
                endif
                
                if( eddyModel ) then
                   revAdj1(icbeg:icend,jcbeg:jcend) = revAdj(0,icbeg:icend,jcbeg:jcend)
                   revAdj2(icbeg:icend,jcbeg:jcend) = revAdj(1,icbeg:icend,jcbeg:jcend)
                endif

             endif
             
             !===========================================================
             
          case (iMax)
             
             ! At most the halo's of nearest neighbors enter the
             ! stencil. Make sure to limit properly.
             
!!$             icBeg = max(jRBeg,jCell-1)
!!$             icEnd = min(jREnd,jCell+1)
!!$             jcBeg = max(kRBeg,kCell-1)
!!$             jcEnd = min(kREnd,kCell+1)
!!$             
!!$             ! Other straightforward stuff.
!!$             
!!$             iOffset = jOffset
!!$             jOffset = kOffset
!!$             
!!$             secondHalo = .true.
!!$             if(iRBeg == iREnd) secondHalo = .false.
             
             if( secondHalo ) then
                wAdj0(icbeg:icend,jcbeg:jcend,1:nw) = wAdj( ib,icbeg:icend,jcbeg:jcend,1:nw)
                wAdj1(icbeg:icend,jcbeg:jcend,1:nw) = wAdj( ib-1,icbeg:icend,jcbeg:jcend,1:nw)
                wAdj2(icbeg:icend,jcbeg:jcend,1:nw) = wAdj( ib-2,icbeg:icend,jcbeg:jcend,1:nw)
                wAdj3(icbeg:icend,jcbeg:jcend,1:nw) = wAdj( ib-3,icbeg:icend,jcbeg:jcend,1:nw)
                
                pAdj0(icbeg:icend,jcbeg:jcend) = pAdj( ib,icbeg:icend,jcbeg:jcend)
                pAdj1(icbeg:icend,jcbeg:jcend) = pAdj( ib-1,icbeg:icend,jcbeg:jcend)
                pAdj2(icbeg:icend,jcbeg:jcend) = pAdj( ib-2,icbeg:icend,jcbeg:jcend)
                pAdj3(icbeg:icend,jcbeg:jcend) = pAdj( ib-3,icbeg:icend,jcbeg:jcend)

                if( viscous ) then
                   rlvAdj1(icbeg:icend,jcbeg:jcend) = rlvAdj(ib-1,icbeg:icend,jcbeg:jcend)
                   rlvAdj2(icbeg:icend,jcbeg:jcend) = rlvAdj(ib-2,icbeg:icend,jcbeg:jcend)
                endif
                if( eddyModel ) then
                   revAdj1(icbeg:icend,jcbeg:jcend) = revAdj(ib-1,icbeg:icend,jcbeg:jcend)
                   revAdj2(icbeg:icend,jcbeg:jcend) = revAdj(ib-2,icbeg:icend,jcbeg:jcend)
                endif
             else
                wAdj1(icbeg:icend,jcbeg:jcend,1:nw) = wAdj( ib,icbeg:icend,jcbeg:jcend,1:nw)
                wAdj2(icbeg:icend,jcbeg:jcend,1:nw) = wAdj( ib-1,icbeg:icend,jcbeg:jcend,1:nw)
                wAdj3(icbeg:icend,jcbeg:jcend,1:nw) = wAdj( ib-2,icbeg:icend,jcbeg:jcend,1:nw)
                
                pAdj1(icbeg:icend,jcbeg:jcend) = pAdj( ib,icbeg:icend,jcbeg:jcend)
                pAdj2(icbeg:icend,jcbeg:jcend) = pAdj( ib-1,icbeg:icend,jcbeg:jcend)
                pAdj3(icbeg:icend,jcbeg:jcend) = pAdj( ib-2,icbeg:icend,jcbeg:jcend)

                if( viscous ) then
                   rlvAdj1(icbeg:icend,jcbeg:jcend) = rlvAdj(ib,icbeg:icend,jcbeg:jcend)
                   rlvAdj2(icbeg:icend,jcbeg:jcend) = rlvAdj(ib-1,icbeg:icend,jcbeg:jcend)
                endif

                if( eddyModel ) then
                   revAdj1(icbeg:icend,jcbeg:jcend) = revAdj(ib,icbeg:icend,jcbeg:jcend)
                   revAdj2(icbeg:icend,jcbeg:jcend) = revAdj(ib-1,icbeg:icend,jcbeg:jcend)
                endif
             endif
             
             !===========================================================
             
          case (jMin)
             !print *,' in jmin'
             ! At most the halo's of nearest neighbors enter the
             ! stencil. Make sure to limit properly.
             
!!$             icBeg = max(iRBeg,iCell-1)
!!$             icEnd = min(iREnd,iCell+1)
!!$             jcBeg = max(kRBeg,kCell-1)
!!$             jcEnd = min(kREnd,kCell+1)
!!$             
!!$             ! Other straightforward stuff.
!!$             
!!$             jOffset = kOffset
!!$             
!!$             secondHalo = .true.
!!$             if(jRBeg == jREnd) secondHalo = .false.
             
             if( secondHalo ) then
                !print *,' in second halo'
                wAdj0(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(icbeg:icend,0,jcbeg:jcend,1:nw)
                !print *,'wadj0'
                wAdj1(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(icbeg:icend,1,jcbeg:jcend,1:nw)
                !print *,'wadj1'
                wAdj2(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(icbeg:icend,2,jcbeg:jcend,1:nw)
                !print *,'wadj2'
                wAdj3(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(icbeg:icend,3,jcbeg:jcend,1:nw)
                !print *,'wadj3'
                
                pAdj0(icbeg:icend,jcbeg:jcend) = pAdj(icbeg:icend,0,jcbeg:jcend)
                !print *,'padj0'
                pAdj1(icbeg:icend,jcbeg:jcend) = pAdj(icbeg:icend,1,jcbeg:jcend)
                !print *,'padj1'
                pAdj2(icbeg:icend,jcbeg:jcend) = pAdj(icbeg:icend,2,jcbeg:jcend)
                !print *,'padj2'
                pAdj3(icbeg:icend,jcbeg:jcend) = pAdj(icbeg:icend,3,jcbeg:jcend)
                !print *,'padj3'

                if( viscous ) then
                   rlvAdj1(icbeg:icend,jcbeg:jcend) = rlvAdj(icbeg:icend,1,jcbeg:jcend)
                   rlvAdj2(icbeg:icend,jcbeg:jcend) = rlvAdj(icbeg:icend,2,jcbeg:jcend)
                endif
                if( eddyModel ) then
                   revAdj1(icbeg:icend,jcbeg:jcend) = revAdj(icbeg:icend,1,jcbeg:jcend)
                   revAdj2(icbeg:icend,jcbeg:jcend) = revAdj(icbeg:icend,2,jcbeg:jcend)
                endif
             else
                wAdj1(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(icbeg:icend,0,jcbeg:jcend,1:nw)
                wAdj2(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(icbeg:icend,1,jcbeg:jcend,1:nw)
                wAdj3(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(icbeg:icend,2,jcbeg:jcend,1:nw)
                
                pAdj1(icbeg:icend,jcbeg:jcend) = pAdj(icbeg:icend,0,jcbeg:jcend)
                pAdj2(icbeg:icend,jcbeg:jcend) = pAdj(icbeg:icend,1,jcbeg:jcend)
                pAdj3(icbeg:icend,jcbeg:jcend) = pAdj(icbeg:icend,2,jcbeg:jcend)
                
                if( viscous ) then
                   rlvAdj1(icbeg:icend,jcbeg:jcend) = rlvAdj(icbeg:icend,0,jcbeg:jcend)
                   rlvAdj2(icbeg:icend,jcbeg:jcend) = rlvAdj(icbeg:icend,1,jcbeg:jcend)
                endif
                if( eddyModel ) then
                   revAdj1(icbeg:icend,jcbeg:jcend) = revAdj(icbeg:icend,0,jcbeg:jcend)
                   revAdj2(icbeg:icend,jcbeg:jcend) = revAdj(icbeg:icend,1,jcbeg:jcend)
                endif

             endif
             
             !===========================================================
             
          case (jMax)
             
             ! At most the halo's of nearest neighbors enter the
             ! stencil. Make sure to limit properly.
             
!!$             icBeg = max(iRBeg,iCell-1)
!!$             icEnd = min(iREnd,iCell+1)
!!$             jcBeg = max(kRBeg,kCell-1)
!!$             jcEnd = min(kREnd,kCell+1)
!!$             
!!$             ! Other straightforward stuff.
!!$             
!!$             jOffset = kOffset
!!$             
!!$             secondHalo = .true.
!!$             if(jRBeg == jREnd) secondHalo = .false.
             
             if( secondHalo ) then
                wAdj0(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(icbeg:icend, jb,jcbeg:jcend,1:nw)
                wAdj1(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(icbeg:icend, jb-1,jcbeg:jcend,1:nw)
                wAdj2(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(icbeg:icend, jb-2,jcbeg:jcend,1:nw)
                wAdj3(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(icbeg:icend, jb-3,jcbeg:jcend,1:nw)
                
                pAdj0(icbeg:icend,jcbeg:jcend) = pAdj(icbeg:icend, jb,jcbeg:jcend)
                pAdj1(icbeg:icend,jcbeg:jcend) = pAdj(icbeg:icend, jb-1,jcbeg:jcend)
                pAdj2(icbeg:icend,jcbeg:jcend) = pAdj(icbeg:icend, jb-2,jcbeg:jcend)
                pAdj3(icbeg:icend,jcbeg:jcend) = pAdj(icbeg:icend, jb-3,jcbeg:jcend)

                if( viscous ) then
                   rlvAdj1(icbeg:icend,jcbeg:jcend) = rlvAdj(icbeg:icend,jb-1,jcbeg:jcend)
                   rlvAdj2(icbeg:icend,jcbeg:jcend) = rlvAdj(icbeg:icend,jb-2,jcbeg:jcend)
                endif
                if( eddyModel ) then
                   revAdj1(icbeg:icend,jcbeg:jcend) = revAdj(icbeg:icend,jb-1,jcbeg:jcend)
                   revAdj2(icbeg:icend,jcbeg:jcend) = revAdj(icbeg:icend,jb-2,jcbeg:jcend)
                endif
             else
                wAdj1(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(icbeg:icend, jb,jcbeg:jcend,1:nw)
                wAdj2(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(icbeg:icend, jb-1,jcbeg:jcend,1:nw)
                wAdj3(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(icbeg:icend, jb-2,jcbeg:jcend,1:nw)
                
                pAdj1(icbeg:icend,jcbeg:jcend) = pAdj(icbeg:icend, jb,jcbeg:jcend)
                pAdj2(icbeg:icend,jcbeg:jcend) = pAdj(icbeg:icend, jb-1,jcbeg:jcend)
                pAdj3(icbeg:icend,jcbeg:jcend) = pAdj(icbeg:icend, jb-2,jcbeg:jcend)

                if( viscous ) then
                   rlvAdj1(icbeg:icend,jcbeg:jcend) = rlvAdj(icbeg:icend,jb,jcbeg:jcend)
                   rlvAdj2(icbeg:icend,jcbeg:jcend) = rlvAdj(icbeg:icend,jb-1,jcbeg:jcend)
                endif
                if( eddyModel ) then
                   revAdj1(icbeg:icend,jcbeg:jcend) = revAdj(icbeg:icend,jb,jcbeg:jcend)
                   revAdj2(icbeg:icend,jcbeg:jcend) = revAdj(icbeg:icend,jb-1,jcbeg:jcend)
                endif
             endif
             
             !===========================================================
             
          case (kMin)
             
             ! At most the halo's of nearest neighbors enter the
             ! stencil. Make sure to limit properly.
             
!!$             icBeg = max(iRBeg,iCell-1)
!!$             icEnd = min(iREnd,iCell+1)
!!$             jcBeg = max(jRBeg,jCell-1)
!!$             jcEnd = min(jREnd,jCell+1)
!!$             
!!$             ! Other straightforward stuff.
!!$             
!!$             secondHalo = .true.
!!$             if(kRBeg == kREnd) secondHalo = .false.
             
             if( secondHalo ) then
                wAdj0(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(icbeg:icend,jcbeg:jcend,0,1:nw)
                wAdj1(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(icbeg:icend,jcbeg:jcend,1,1:nw)
                wAdj2(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(icbeg:icend,jcbeg:jcend,2,1:nw)
                wAdj3(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(icbeg:icend,jcbeg:jcend,3,1:nw)
                
                pAdj0(icbeg:icend,jcbeg:jcend) = pAdj(icbeg:icend,jcbeg:jcend,0)
                pAdj1(icbeg:icend,jcbeg:jcend) = pAdj(icbeg:icend,jcbeg:jcend,1)
                pAdj2(icbeg:icend,jcbeg:jcend) = pAdj(icbeg:icend,jcbeg:jcend,2)
                pAdj3(icbeg:icend,jcbeg:jcend) = pAdj(icbeg:icend,jcbeg:jcend,3)

                if( viscous ) then
                   rlvAdj1(icbeg:icend,jcbeg:jcend) = rlvAdj(icbeg:icend,jcbeg:jcend,1)
                   rlvAdj2(icbeg:icend,jcbeg:jcend) = rlvAdj(icbeg:icend,jcbeg:jcend,2)
                endif
                if( eddyModel ) then
                   revAdj1(icbeg:icend,jcbeg:jcend) = revAdj(icbeg:icend,jcbeg:jcend,1)
                   revAdj2(icbeg:icend,jcbeg:jcend) = revAdj(icbeg:icend,jcbeg:jcend,2)
                endif
             else
                wAdj1(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(icbeg:icend,jcbeg:jcend,0,1:nw)
                wAdj2(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(icbeg:icend,jcbeg:jcend,1,1:nw)
                wAdj3(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(icbeg:icend,jcbeg:jcend,2,1:nw)
                
                pAdj1(icbeg:icend,jcbeg:jcend) = pAdj(icbeg:icend,jcbeg:jcend,0)
                pAdj2(icbeg:icend,jcbeg:jcend) = pAdj(icbeg:icend,jcbeg:jcend,1)
                pAdj3(icbeg:icend,jcbeg:jcend) = pAdj(icbeg:icend,jcbeg:jcend,2)

                
                if( viscous ) then
                   rlvAdj1(icbeg:icend,jcbeg:jcend) = rlvAdj(icbeg:icend,jcbeg:jcend,0)
                   rlvAdj2(icbeg:icend,jcbeg:jcend) = rlvAdj(icbeg:icend,jcbeg:jcend,1)
                endif
                if( eddyModel ) then
                   revAdj1(icbeg:icend,jcbeg:jcend) = revAdj(icbeg:icend,jcbeg:jcend,0)
                   revAdj2(icbeg:icend,jcbeg:jcend) = revAdj(icbeg:icend,jcbeg:jcend,1)
                endif
             endif
             
             !===========================================================
             
          case (kMax)
             
             ! At most the halo's of nearest neighbors enter the
             ! stencil. Make sure to limit properly.
             
!!$             icBeg = max(iRBeg,iCell-1)
!!$             icEnd = min(iREnd,iCell+1)
!!$             jcBeg = max(jRBeg,jCell-1)
!!$             jcEnd = min(jREnd,jCell+1)
!!$             
!!$             ! Other straightforward stuff.
!!$             
!!$             secondHalo = .true.
!!$             if(kRBeg == kREnd) secondHalo = .false.
             
             if( secondHalo ) then
                wAdj0(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(icbeg:icend,jcbeg:jcend, kb,1:nw)
                wAdj1(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(icbeg:icend,jcbeg:jcend, kb-1,1:nw)
                wAdj2(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(icbeg:icend,jcbeg:jcend, kb-2,1:nw)
                wAdj3(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(icbeg:icend,jcbeg:jcend,kb-3,1:nw)
                
                pAdj0(icbeg:icend,jcbeg:jcend) = pAdj(icbeg:icend,jcbeg:jcend, kb)
                pAdj1(icbeg:icend,jcbeg:jcend) = pAdj(icbeg:icend,jcbeg:jcend, kb-1)
                pAdj2(icbeg:icend,jcbeg:jcend) = pAdj(icbeg:icend,jcbeg:jcend, kb-2)
                pAdj3(icbeg:icend,jcbeg:jcend) = pAdj(icbeg:icend,jcbeg:jcend, kb-3)

                
                if( viscous ) then
                   rlvAdj1(icbeg:icend,jcbeg:jcend) = rlvAdj(icbeg:icend,jcbeg:jcend,kb-1)
                   rlvAdj2(icbeg:icend,jcbeg:jcend) = rlvAdj(icbeg:icend,jcbeg:jcend,kb-2)
                endif
                if( eddyModel ) then
                   revAdj1(icbeg:icend,jcbeg:jcend) = revAdj(icbeg:icend,jcbeg:jcend,kb-1)
                   revAdj2(icbeg:icend,jcbeg:jcend) = revAdj(icbeg:icend,jcbeg:jcend,kb-2)
                endif
             else
                wAdj1(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(icbeg:icend,jcbeg:jcend, kb,1:nw)
                wAdj2(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(icbeg:icend,jcbeg:jcend, kb-1,1:nw)
                wAdj3(icbeg:icend,jcbeg:jcend,1:nw) = wAdj(icbeg:icend,jcbeg:jcend, kb-2,1:nw)
                
                pAdj1(icbeg:icend,jcbeg:jcend) = pAdj(icbeg:icend,jcbeg:jcend, kb)
                pAdj2(icbeg:icend,jcbeg:jcend) = pAdj(icbeg:icend,jcbeg:jcend, kb-1)
                pAdj3(icbeg:icend,jcbeg:jcend) = pAdj(icbeg:icend,jcbeg:jcend, kb-2)

                
                if( viscous ) then
                   rlvAdj1(icbeg:icend,jcbeg:jcend) = rlvAdj(icbeg:icend,jcbeg:jcend,kb)
                   rlvAdj2(icbeg:icend,jcbeg:jcend) = rlvAdj(icbeg:icend,jcbeg:jcend,kb-1)
                endif
                if( eddyModel ) then
                   revAdj1(icbeg:icend,jcbeg:jcend) = revAdj(icbeg:icend,jcbeg:jcend,kb)
                   revAdj2(icbeg:icend,jcbeg:jcend) = revAdj(icbeg:icend,jcbeg:jcend,kb-1)
                endif
             endif
             
          end select
          
         
        end subroutine extractBCStatesForcesAdj
  
