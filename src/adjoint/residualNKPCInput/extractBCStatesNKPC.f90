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
subroutine extractBCStatesNKPC(nn,wAdj,pAdj,wAdj0, wAdj1, wAdj2,wAdj3,&
     pAdj0,pAdj1, pAdj2,pAdj3,&
     rlvAdj, revAdj,rlvAdj1, rlvAdj2,revAdj1, revAdj2,iOffset,&
     jOffset, kOffset,iCell, jCell,kCell,&
     isbeg,jsbeg,ksbeg,isend,jsend,ksend,ibbeg,jbbeg,kbbeg,ibend,&
     jbend,kbend,icbeg,jcbeg,icend,jcend,secondHalo,nnn,level,sps,sps2)

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
  use inputTimeSpectral !nIntervalTimespectral

  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: nn,nnn,level,sps,sps2!, offset
  integer(kind=intType), intent(in) ::iCell, jCell,kCell
  integer(kind=intType), intent(in) ::isbeg,jsbeg,ksbeg,isend,jsend,ksend
  integer(kind=intType), intent(in) ::ibbeg,jbbeg,kbbeg,ibend,jbend,kbend
  integer(kind=intType) ::icbeg,jcbeg,icend,jcend
  integer(kind=intType) ::irbeg,jrbeg,krbeg,irend,jrend,krend

  integer(kind=intType) :: iOffset, jOffset, kOffset

  real(kind=realType), dimension(-2:2,-2:2,-2:2,nw,nTimeIntervalsSpectral), &
       intent(in) :: wAdj
  real(kind=realType), dimension(-2:2,-2:2,-2:2,nTimeIntervalsSpectral),intent(in) :: pAdj


  real(kind=realType), dimension(-2:2,-2:2,nw) :: wAdj0, wAdj1
  real(kind=realType), dimension(-2:2,-2:2,nw) :: wAdj2, wAdj3
  real(kind=realType), dimension(-2:2,-2:2)    :: pAdj0, pAdj1
  real(kind=realType), dimension(-2:2,-2:2)    :: pAdj2, pAdj3

  real(kind=realType), dimension(-2:2,-2:2,-2:2)::rlvAdj, revAdj
  real(kind=realType), dimension(-2:2,-2:2)::rlvAdj1, rlvAdj2
  real(kind=realType), dimension(-2:2,-2:2)::revAdj1, revAdj2

  logical :: iOverlap, jOverlap, kOverlap
  !logical,intent(inout) :: secondHalo
  logical :: secondHalo


  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! There is an overlap between the two ranges.
  ! Determine the actual overlap region.

  iRBeg = max(iSBeg, iBBeg); iREnd = min(iSEnd, iBEnd)
  jRBeg = max(jSBeg, jBBeg); jREnd = min(jSEnd, jBEnd)
  kRBeg = max(kSBeg, kBBeg); kREnd = min(kSEnd, kBEnd)
  !          print *,'stencil',kSBeg, kBBeg,krbeg,kSEnd, kBEnd,krend
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

     icBeg = max(jRBeg,jCell-2)
     icEnd = min(jREnd,jCell+2)
     jcBeg = max(kRBeg,kCell-2)
     jcEnd = min(kREnd,kCell+2)


     ! Other straightforward stuff.

     iOffset = jOffset
     jOffset = kOffset

     secondHalo = .true.
     if(iRBeg == iREnd) secondHalo = .false.

     if( secondHalo ) then
        wAdj0(-2:2,-2:2,1:nw) = wAdj(-2,-2:2,-2:2,1:nw,sps2)
        wAdj1(-2:2,-2:2,1:nw) = wAdj(-1,-2:2,-2:2,1:nw,sps2)
        wAdj2(-2:2,-2:2,1:nw) = wAdj( 0,-2:2,-2:2,1:nw,sps2)
        wAdj3(-2:2,-2:2,1:nw) = wAdj( 1,-2:2,-2:2,1:nw,sps2)

        pAdj0(-2:2,-2:2) = pAdj(-2,-2:2,-2:2,sps2)
        pAdj1(-2:2,-2:2) = pAdj(-1,-2:2,-2:2,sps2)
        pAdj2(-2:2,-2:2) = pAdj( 0,-2:2,-2:2,sps2)
        pAdj3(-2:2,-2:2) = pAdj( 1,-2:2,-2:2,sps2)

        if( viscous ) then
           rlvAdj1(-2:2,-2:2) = rlvAdj(-1,-2:2,-2:2)
           rlvAdj2(-2:2,-2:2) = rlvAdj( 0,-2:2,-2:2)
        endif
        if( eddyModel ) then
           revAdj1(-2:2,-2:2) = revAdj(-1,-2:2,-2:2)
           revAdj2(-2:2,-2:2) = revAdj( 0,-2:2,-2:2)
        endif

     else
        wAdj1(-2:2,-2:2,1:nw) = wAdj(-2,-2:2,-2:2,1:nw,sps2)
        wAdj2(-2:2,-2:2,1:nw) = wAdj(-1,-2:2,-2:2,1:nw,sps2)
        wAdj3(-2:2,-2:2,1:nw) = wAdj( 0,-2:2,-2:2,1:nw,sps2)

        pAdj1(-2:2,-2:2) = pAdj(-2,-2:2,-2:2,sps2)
        pAdj2(-2:2,-2:2) = pAdj(-1,-2:2,-2:2,sps2)
        pAdj3(-2:2,-2:2) = pAdj( 0,-2:2,-2:2,sps2)

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

     icBeg = max(jRBeg,jCell-2)
     icEnd = min(jREnd,jCell+2)
     jcBeg = max(kRBeg,kCell-2)
     jcEnd = min(kREnd,kCell+2)

     ! Other straightforward stuff.

     iOffset = jOffset
     jOffset = kOffset

     secondHalo = .true.
     if(iRBeg == iREnd) secondHalo = .false.

     if( secondHalo ) then
        wAdj0(-2:2,-2:2,1:nw) = wAdj( 2,-2:2,-2:2,1:nw,sps2)
        wAdj1(-2:2,-2:2,1:nw) = wAdj( 1,-2:2,-2:2,1:nw,sps2)
        wAdj2(-2:2,-2:2,1:nw) = wAdj( 0,-2:2,-2:2,1:nw,sps2)
        wAdj3(-2:2,-2:2,1:nw) = wAdj(-1,-2:2,-2:2,1:nw,sps2)

        pAdj0(-2:2,-2:2) = pAdj( 2,-2:2,-2:2,sps2)
        pAdj1(-2:2,-2:2) = pAdj( 1,-2:2,-2:2,sps2)
        pAdj2(-2:2,-2:2) = pAdj( 0,-2:2,-2:2,sps2)
        pAdj3(-2:2,-2:2) = pAdj(-1,-2:2,-2:2,sps2)

        if( viscous ) then
           rlvAdj1(-2:2,-2:2) = rlvAdj(1,-2:2,-2:2)
           rlvAdj2(-2:2,-2:2) = rlvAdj(0,-2:2,-2:2)
        endif
        if( eddyModel ) then
           revAdj1(-2:2,-2:2) = revAdj(1,-2:2,-2:2)
           revAdj2(-2:2,-2:2) = revAdj(0,-2:2,-2:2)
        endif
     else
        wAdj1(-2:2,-2:2,1:nw) = wAdj( 2,-2:2,-2:2,1:nw,sps2)
        wAdj2(-2:2,-2:2,1:nw) = wAdj( 1,-2:2,-2:2,1:nw,sps2)
        wAdj3(-2:2,-2:2,1:nw) = wAdj( 0,-2:2,-2:2,1:nw,sps2)

        pAdj1(-2:2,-2:2) = pAdj( 2,-2:2,-2:2,sps2)
        pAdj2(-2:2,-2:2) = pAdj( 1,-2:2,-2:2,sps2)
        pAdj3(-2:2,-2:2) = pAdj( 0,-2:2,-2:2,sps2)

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

     icBeg = max(iRBeg,iCell-2)
     icEnd = min(iREnd,iCell+2)
     jcBeg = max(kRBeg,kCell-2)
     jcEnd = min(kREnd,kCell+2)
     ! Other straightforward stuff.

     jOffset = kOffset

     secondHalo = .true.
     if(jRBeg == jREnd) secondHalo = .false.

     if( secondHalo ) then
        wAdj0(-2:2,-2:2,1:nw) = wAdj(-2:2,-2,-2:2,1:nw,sps2)
        wAdj1(-2:2,-2:2,1:nw) = wAdj(-2:2,-1,-2:2,1:nw,sps2)
        wAdj2(-2:2,-2:2,1:nw) = wAdj(-2:2, 0,-2:2,1:nw,sps2)
        wAdj3(-2:2,-2:2,1:nw) = wAdj(-2:2, 1,-2:2,1:nw,sps2)

        pAdj0(-2:2,-2:2) = pAdj(-2:2,-2,-2:2,sps2)
        pAdj1(-2:2,-2:2) = pAdj(-2:2,-1,-2:2,sps2)
        pAdj2(-2:2,-2:2) = pAdj(-2:2, 0,-2:2,sps2)
        pAdj3(-2:2,-2:2) = pAdj(-2:2, 1,-2:2,sps2)

        if( viscous ) then
           rlvAdj1(-2:2,-2:2) = rlvAdj(-2:2,-1,-2:2)
           rlvAdj2(-2:2,-2:2) = rlvAdj(-2:2,0,-2:2)
        endif
        if( eddyModel ) then
           revAdj1(-2:2,-2:2) = revAdj(-2:2,-1,-2:2)
           revAdj2(-2:2,-2:2) = revAdj(-2:2,0,-2:2)
        endif
     else
        wAdj1(-2:2,-2:2,1:nw) = wAdj(-2:2,-2,-2:2,1:nw,sps2)
        wAdj2(-2:2,-2:2,1:nw) = wAdj(-2:2,-1,-2:2,1:nw,sps2)
        wAdj3(-2:2,-2:2,1:nw) = wAdj(-2:2, 0,-2:2,1:nw,sps2)

        pAdj1(-2:2,-2:2) = pAdj(-2:2,-2,-2:2,sps2)
        pAdj2(-2:2,-2:2) = pAdj(-2:2,-1,-2:2,sps2)
        pAdj3(-2:2,-2:2) = pAdj(-2:2, 0,-2:2,sps2)

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

     icBeg = max(iRBeg,iCell-2)
     icEnd = min(iREnd,iCell+2)
     jcBeg = max(kRBeg,kCell-2)
     jcEnd = min(kREnd,kCell+2)

     ! Other straightforward stuff.

     jOffset = kOffset

     secondHalo = .true.
     if(jRBeg == jREnd) secondHalo = .false.

     if( secondHalo ) then
        wAdj0(-2:2,-2:2,1:nw) = wAdj(-2:2, 2,-2:2,1:nw,sps2)
        wAdj1(-2:2,-2:2,1:nw) = wAdj(-2:2, 1,-2:2,1:nw,sps2)
        wAdj2(-2:2,-2:2,1:nw) = wAdj(-2:2, 0,-2:2,1:nw,sps2)
        wAdj3(-2:2,-2:2,1:nw) = wAdj(-2:2,-1,-2:2,1:nw,sps2)

        pAdj0(-2:2,-2:2) = pAdj(-2:2, 2,-2:2,sps2)
        pAdj1(-2:2,-2:2) = pAdj(-2:2, 1,-2:2,sps2)
        pAdj2(-2:2,-2:2) = pAdj(-2:2, 0,-2:2,sps2)
        pAdj3(-2:2,-2:2) = pAdj(-2:2,-1,-2:2,sps2)

        if( viscous ) then
           rlvAdj1(-2:2,-2:2) = rlvAdj(-2:2,1,-2:2)
           rlvAdj2(-2:2,-2:2) = rlvAdj(-2:2,0,-2:2)
        endif
        if( eddyModel ) then
           revAdj1(-2:2,-2:2) = revAdj(-2:2,1,-2:2)
           revAdj2(-2:2,-2:2) = revAdj(-2:2,0,-2:2)
        endif
     else
        wAdj1(-2:2,-2:2,1:nw) = wAdj(-2:2, 2,-2:2,1:nw,sps2)
        wAdj2(-2:2,-2:2,1:nw) = wAdj(-2:2, 1,-2:2,1:nw,sps2)
        wAdj3(-2:2,-2:2,1:nw) = wAdj(-2:2, 0,-2:2,1:nw,sps2)

        pAdj1(-2:2,-2:2) = pAdj(-2:2, 2,-2:2,sps2)
        pAdj2(-2:2,-2:2) = pAdj(-2:2, 1,-2:2,sps2)
        pAdj3(-2:2,-2:2) = pAdj(-2:2, 0,-2:2,sps2)

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


     icBeg = max(iRBeg,iCell-2)
     icEnd = min(iREnd,iCell+2)
     jcBeg = max(jRBeg,jCell-2)
     jcEnd = min(jREnd,jCell+2)

     ! Other straightforward stuff.

     secondHalo = .true.
     if(kRBeg == kREnd) secondHalo = .false.

     if( secondHalo ) then
        wAdj0(-2:2,-2:2,1:nw) = wAdj(-2:2,-2:2,-2,1:nw,sps2)
        wAdj1(-2:2,-2:2,1:nw) = wAdj(-2:2,-2:2,-1,1:nw,sps2)
        wAdj2(-2:2,-2:2,1:nw) = wAdj(-2:2,-2:2, 0,1:nw,sps2)
        wAdj3(-2:2,-2:2,1:nw) = wAdj(-2:2,-2:2, 1,1:nw,sps2)

        pAdj0(-2:2,-2:2) = pAdj(-2:2,-2:2,-2,sps2)
        pAdj1(-2:2,-2:2) = pAdj(-2:2,-2:2,-1,sps2)
        pAdj2(-2:2,-2:2) = pAdj(-2:2,-2:2, 0,sps2)
        pAdj3(-2:2,-2:2) = pAdj(-2:2,-2:2, 1,sps2)

        if( viscous ) then
           rlvAdj1(-2:2,-2:2) = rlvAdj(-2:2,-2:2,-1)
           rlvAdj2(-2:2,-2:2) = rlvAdj(-2:2,-2:2,0)
        endif
        if( eddyModel ) then
           revAdj1(-2:2,-2:2) = revAdj(-2:2,-2:2,-1)
           revAdj2(-2:2,-2:2) = revAdj(-2:2,-2:2,0)
        endif
     else
        wAdj1(-2:2,-2:2,1:nw) = wAdj(-2:2,-2:2,-2,1:nw,sps2)
        wAdj2(-2:2,-2:2,1:nw) = wAdj(-2:2,-2:2,-1,1:nw,sps2)
        wAdj3(-2:2,-2:2,1:nw) = wAdj(-2:2,-2:2, 0,1:nw,sps2)

        pAdj1(-2:2,-2:2) = pAdj(-2:2,-2:2,-2,sps2)
        pAdj2(-2:2,-2:2) = pAdj(-2:2,-2:2,-1,sps2)
        pAdj3(-2:2,-2:2) = pAdj(-2:2,-2:2, 0,sps2)


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

     icBeg = max(iRBeg,iCell-2)
     icEnd = min(iREnd,iCell+2)
     jcBeg = max(jRBeg,jCell-2)
     jcEnd = min(jREnd,jCell+2)

     ! Other straightforward stuff.

     secondHalo = .true.
     if(kRBeg == kREnd) secondHalo = .false.

     if( secondHalo ) then
        wAdj0(-2:2,-2:2,1:nw) = wAdj(-2:2,-2:2, 2,1:nw,sps2)
        wAdj1(-2:2,-2:2,1:nw) = wAdj(-2:2,-2:2, 1,1:nw,sps2)
        wAdj2(-2:2,-2:2,1:nw) = wAdj(-2:2,-2:2, 0,1:nw,sps2)
        wAdj3(-2:2,-2:2,1:nw) = wAdj(-2:2,-2:2,-1,1:nw,sps2)

        pAdj0(-2:2,-2:2) = pAdj(-2:2,-2:2, 2,sps2)
        pAdj1(-2:2,-2:2) = pAdj(-2:2,-2:2, 1,sps2)
        pAdj2(-2:2,-2:2) = pAdj(-2:2,-2:2, 0,sps2)
        pAdj3(-2:2,-2:2) = pAdj(-2:2,-2:2,-1,sps2)


        if( viscous ) then
           rlvAdj1(-2:2,-2:2) = rlvAdj(-2:2,-2:2,1)
           rlvAdj2(-2:2,-2:2) = rlvAdj(-2:2,-2:2,0)
        endif
        if( eddyModel ) then
           revAdj1(-2:2,-2:2) = revAdj(-2:2,-2:2,1)
           revAdj2(-2:2,-2:2) = revAdj(-2:2,-2:2,0)
        endif
     else
        wAdj1(-2:2,-2:2,1:nw) = wAdj(-2:2,-2:2, 2,1:nw,sps2)
        wAdj2(-2:2,-2:2,1:nw) = wAdj(-2:2,-2:2, 1,1:nw,sps2)
        wAdj3(-2:2,-2:2,1:nw) = wAdj(-2:2,-2:2, 0,1:nw,sps2)

        pAdj1(-2:2,-2:2) = pAdj(-2:2,-2:2, 2,sps2)
        pAdj2(-2:2,-2:2) = pAdj(-2:2,-2:2, 1,sps2)
        pAdj3(-2:2,-2:2) = pAdj(-2:2,-2:2, 0,sps2)


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

end subroutine extractBCStatesNKPC

