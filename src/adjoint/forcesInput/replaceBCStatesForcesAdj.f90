!
!      ******************************************************************
!      *                                                                *
!      * File:          replaceBCStatesForcesAdj.f90                    *
!      * Author:        C.A.(Sandy) Mader                               *
!      * Starting date: 04-17-2008                                      *
!      * Last modified: 06-09-2008                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine replaceBCStatesForcesAdj(mm,  wAdj0,wAdj1, wAdj2, wAdj3,&
            pAdj0,pAdj1, pAdj2, pAdj3,rlvAdj1, rlvAdj2,revAdj1, revAdj2,&
            wAdj,pAdj,rlvAdj,revAdj,&
            icbeg,jcbeg,icend,jcend,secondHalo)

!      modules
       use BCTypes
       use blockPointers, only : ie, ib, je, jb, ke, kb, nBocos, &
                                  BCFaceID, BCType, BCData
       use flowVarRefState
       implicit none

       integer(kind=intType), intent(in) :: mm
       integer(kind=intType) ::icbeg,jcbeg,icend,jcend
 
      
       real(kind=realType), dimension(0:ib,0:jb,0:kb,nw), intent(inout) :: wAdj
       real(kind=realType), dimension(0:ib,0:jb,0:kb), intent(inout) :: pAdj
       logical, intent(in) :: secondHalo

       real(kind=realType), dimension(icbeg:icend,jcbeg:jcend,nw):: wAdj0, wAdj1
       real(kind=realType), dimension(icbeg:icend,jcbeg:jcend,nw):: wAdj2, wAdj3
       real(kind=realType), dimension(icbeg:icend,jcbeg:jcend) :: pAdj0, pAdj1
       real(kind=realType), dimension(icbeg:icend,jcbeg:jcend) :: pAdj2, pAdj3

       real(kind=realType), dimension(0:ib,0:jb,0:kb)::rlvAdj, revAdj
       real(kind=realType), dimension(icbeg:icend,jcbeg:jcend)::rlvAdj1, rlvAdj2
       real(kind=realType), dimension(icbeg:icend,jcbeg:jcend)::revAdj1, revAdj2
!!$       real(kind=realType), dimension(-2:2,-2:2,nw),intent(in) :: wAdj0, wAdj1
!!$       real(kind=realType), dimension(-2:2,-2:2,nw),intent(in) :: wAdj2, wAdj3
!!$       real(kind=realType), dimension(-2:2,-2:2),intent(in)    :: pAdj0, pAdj1
!!$       real(kind=realType), dimension(-2:2,-2:2),intent(in)    :: pAdj2, pAdj3
!!$
!!$       real(kind=realType), dimension(-2:2,-2:2,-2:2,nw),intent(inout) :: wAdj
!!$       real(kind=realType), dimension(-2:2,-2:2,-2:2),intent(inout) :: pAdj
       
       integer(kind=intType) ::i,j,k,l
!!$       real(kind=realType), dimension(-2:2,-2:2,-2:2)::rlvAdj, revAdj
!!$       real(kind=realType), dimension(-2:2,-2:2)::rlvAdj1, rlvAdj2
!!$       real(kind=realType), dimension(-2:2,-2:2)::revAdj1, revAdj2
       
       !logical :: secondHalo

       ! Copy the information back to the original arrays wAdj and pAdj.

       select case (BCFaceID(mm))
       case (iMin)
          
          if( secondHalo ) then

             wAdj(0,icbeg:icend,jcbeg:jcend,1:nw) = wAdj0(icbeg:icend,jcbeg:jcend,1:nw) 
             wAdj(1,icbeg:icend,jcbeg:jcend,1:nw) = wAdj1(icbeg:icend,jcbeg:jcend,1:nw)
             
             pAdj(0,icbeg:icend,jcbeg:jcend) = pAdj0(icbeg:icend,jcbeg:jcend) 
             pAdj(1,icbeg:icend,jcbeg:jcend) = pAdj1(icbeg:icend,jcbeg:jcend) 

             if( viscous ) then
                rlvAdj(1,icbeg:icend,jcbeg:jcend) = rlvAdj1(icbeg:icend,jcbeg:jcend) 
                rlvAdj(2,icbeg:icend,jcbeg:jcend) = rlvAdj2(icbeg:icend,jcbeg:jcend) 
             endif
             if( eddyModel ) then
                revAdj(1,icbeg:icend,jcbeg:jcend) = revAdj1(icbeg:icend,jcbeg:jcend) 
                revAdj(2,icbeg:icend,jcbeg:jcend) = revAdj2(icbeg:icend,jcbeg:jcend)
             endif
          else
             wAdj(0,icbeg:icend,jcbeg:jcend,1:nw) = wAdj1(icbeg:icend,jcbeg:jcend,1:nw)
             pAdj(0,icbeg:icend,jcbeg:jcend)      = pAdj1(icbeg:icend,jcbeg:jcend)
          endif
          
          !===========================================================
          
       case (iMax)
          
          if( secondHalo ) then

             wAdj( ib,icbeg:icend,jcbeg:jcend,1:nw) = wAdj0(icbeg:icend,jcbeg:jcend,1:nw)
             wAdj( ib-1,icbeg:icend,jcbeg:jcend,1:nw) = wAdj1(icbeg:icend,jcbeg:jcend,1:nw)
             
             pAdj( ib,icbeg:icend,jcbeg:jcend) = pAdj0(icbeg:icend,jcbeg:jcend)
             pAdj( ib-1,icbeg:icend,jcbeg:jcend) = pAdj1(icbeg:icend,jcbeg:jcend)
          else
             wAdj( ib,icbeg:icend,jcbeg:jcend,1:nw) = wAdj1(icbeg:icend,jcbeg:jcend,1:nw)
             pAdj( ib,icbeg:icend,jcbeg:jcend)      = pAdj1(icbeg:icend,jcbeg:jcend)
          endif
          
          !===========================================================
          
       case (jMin)

          if( secondHalo ) then
             wAdj(icbeg:icend,0,jcbeg:jcend,1:nw) = wAdj0(icbeg:icend,jcbeg:jcend,1:nw)
             wAdj(icbeg:icend,1,jcbeg:jcend,1:nw) = wAdj1(icbeg:icend,jcbeg:jcend,1:nw)
             
             pAdj(icbeg:icend,0,jcbeg:jcend) = pAdj0(icbeg:icend,jcbeg:jcend)
             pAdj(icbeg:icend,1,jcbeg:jcend) = pAdj1(icbeg:icend,jcbeg:jcend)
          else
             wAdj(icbeg:icend,0,jcbeg:jcend,1:nw) = wAdj1(icbeg:icend,jcbeg:jcend,1:nw)
             pAdj(icbeg:icend,0,jcbeg:jcend)      = pAdj1(icbeg:icend,jcbeg:jcend)
          endif
          
          !===========================================================
          
       case (jMax)
          
          if( secondHalo ) then
             wAdj(icbeg:icend, jb,jcbeg:jcend,1:nw) = wAdj0(icbeg:icend,jcbeg:jcend,1:nw)
             wAdj(icbeg:icend, jb-1,jcbeg:jcend,1:nw) = wAdj1(icbeg:icend,jcbeg:jcend,1:nw)
             
             pAdj(icbeg:icend, jb,jcbeg:jcend) = pAdj0(icbeg:icend,jcbeg:jcend)
             pAdj(icbeg:icend, jb-1,jcbeg:jcend) = pAdj1(icbeg:icend,jcbeg:jcend)
          else
             wAdj(icbeg:icend, jb,jcbeg:jcend,1:nw) = wAdj1(icbeg:icend,jcbeg:jcend,1:nw)
             pAdj(icbeg:icend, jb,jcbeg:jcend)      = pAdj1(icbeg:icend,jcbeg:jcend)
          endif
          
          !===========================================================
          
       case (kMin)
          
          if( secondHalo ) then
             wAdj(icbeg:icend,jcbeg:jcend,0,1:nw) = wAdj0(icbeg:icend,jcbeg:jcend,1:nw)
             wAdj(icbeg:icend,jcbeg:jcend,1,1:nw) = wAdj1(icbeg:icend,jcbeg:jcend,1:nw)

             pAdj(icbeg:icend,jcbeg:jcend,0) = pAdj0(icbeg:icend,jcbeg:jcend)
             pAdj(icbeg:icend,jcbeg:jcend,1) = pAdj1(icbeg:icend,jcbeg:jcend)
          else
             wAdj(icbeg:icend,jcbeg:jcend,0,1:nw) = wAdj1(icbeg:icend,jcbeg:jcend,1:nw)
             pAdj(icbeg:icend,jcbeg:jcend,0)      = pAdj1(icbeg:icend,jcbeg:jcend)
          endif
          
          !===========================================================
          
       case (kMax)
          
          if( secondHalo ) then
             wAdj(icbeg:icend,jcbeg:jcend, kb,1:nw) = wAdj0(icbeg:icend,jcbeg:jcend,1:nw)
             wAdj(icbeg:icend,jcbeg:jcend, kb-1,1:nw) = wAdj1(icbeg:icend,jcbeg:jcend,1:nw)
             
             pAdj(icbeg:icend,jcbeg:jcend, kb) = pAdj0(icbeg:icend,jcbeg:jcend)
             pAdj(icbeg:icend,jcbeg:jcend, kb-1) = pAdj1(icbeg:icend,jcbeg:jcend)
          else
             wAdj(icbeg:icend,jcbeg:jcend, kb,1:nw) = wAdj1(icbeg:icend,jcbeg:jcend,1:nw)
             pAdj(icbeg:icend,jcbeg:jcend, kb)      = pAdj1(icbeg:icend,jcbeg:jcend)
          endif
          
       end select
       
     end subroutine replaceBCStatesForcesAdj
