!
!      ******************************************************************
!      *                                                                *
!      * File:          replaceBCStatesAdj.f90                          *
!      * Author:        C.A.(Sandy) Mader                               *
!      * Starting date: 04-17-2008                                      *
!      * Last modified: 04-17-2008                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine replaceBCStatesNKPC(nn,  wAdj0,wAdj1, wAdj2, wAdj3,&
     pAdj0,pAdj1, pAdj2, pAdj3,rlvAdj1, rlvAdj2,revAdj1, revAdj2,&
     iCell, jCell,kCell,&
     wAdj,pAdj,rlvAdj,revAdj,secondHalo,nnn,level,sps,sps2)

  !      modules
  use BCTypes
  use blockPointers, only : ie, ib, je, jb, ke, kb, nBocos, &
       BCFaceID, BCType, BCData
  use flowVarRefState
  use inputTimeSpectral !nIntervalTimespectral
  implicit none

  integer(kind=intType), intent(in) :: nn,nnn,level,sps,sps2
  real(kind=realType), dimension(-2:2,-2:2,nw),intent(in) :: wAdj0, wAdj1
  real(kind=realType), dimension(-2:2,-2:2,nw),intent(in) :: wAdj2, wAdj3
  real(kind=realType), dimension(-2:2,-2:2),intent(in)    :: pAdj0, pAdj1
  real(kind=realType), dimension(-2:2,-2:2),intent(in)    :: pAdj2, pAdj3

  real(kind=realType), dimension(-2:2,-2:2,-2:2,nw,nTimeIntervalsSpectral),intent(inout) :: wAdj
  real(kind=realType), dimension(-2:2,-2:2,-2:2,nTimeIntervalsSpectral),intent(inout) :: pAdj

  integer(kind=intType) ::iCell, jCell,kCell,i,j,k,l
  real(kind=realType), dimension(-2:2,-2:2,-2:2)::rlvAdj, revAdj
  real(kind=realType), dimension(-2:2,-2:2)::rlvAdj1, rlvAdj2
  real(kind=realType), dimension(-2:2,-2:2)::revAdj1, revAdj2

  logical :: secondHalo

  ! Copy the information back to the original arrays wAdj and pAdj.

  select case (BCFaceID(nn))
  case (iMin)

     if( secondHalo ) then

        wAdj(-2,-2:2,-2:2,1:nw,sps2) = wAdj0(-2:2,-2:2,1:nw) 
        wAdj(-1,-2:2,-2:2,1:nw,sps2) = wAdj1(-2:2,-2:2,1:nw)

        pAdj(-2,-2:2,-2:2,sps2) = pAdj0(-2:2,-2:2) 
        pAdj(-1,-2:2,-2:2,sps2) = pAdj1(-2:2,-2:2) 

        if( viscous ) then
           rlvAdj(-1,-2:2,-2:2) = rlvAdj1(-2:2,-2:2) 
           rlvAdj( 0,-2:2,-2:2) = rlvAdj2(-2:2,-2:2) 
        endif
        if( eddyModel ) then
           revAdj(-1,-2:2,-2:2) = revAdj1(-2:2,-2:2) 
           revAdj( 0,-2:2,-2:2) = revAdj2(-2:2,-2:2)
        endif
     else
        wAdj(-2,-2:2,-2:2,1:nw,sps2) = wAdj1(-2:2,-2:2,1:nw)
        pAdj(-2,-2:2,-2:2,sps2)      = pAdj1(-2:2,-2:2)
     endif

     !===========================================================

  case (iMax)

     if( secondHalo ) then

        wAdj( 2,-2:2,-2:2,1:nw,sps2) = wAdj0(-2:2,-2:2,1:nw)
        wAdj( 1,-2:2,-2:2,1:nw,sps2) = wAdj1(-2:2,-2:2,1:nw)

        pAdj( 2,-2:2,-2:2,sps2) = pAdj0(-2:2,-2:2)
        pAdj( 1,-2:2,-2:2,sps2) = pAdj1(-2:2,-2:2)
     else
        wAdj( 2,-2:2,-2:2,1:nw,sps2) = wAdj1(-2:2,-2:2,1:nw)
        pAdj( 2,-2:2,-2:2,sps2)      = pAdj1(-2:2,-2:2)
     endif

     !===========================================================

  case (jMin)

     if( secondHalo ) then
        wAdj(-2:2,-2,-2:2,1:nw,sps2) = wAdj0(-2:2,-2:2,1:nw)
        wAdj(-2:2,-1,-2:2,1:nw,sps2) = wAdj1(-2:2,-2:2,1:nw)

        pAdj(-2:2,-2,-2:2,sps2) = pAdj0(-2:2,-2:2)
        pAdj(-2:2,-1,-2:2,sps2) = pAdj1(-2:2,-2:2)
     else
        wAdj(-2:2,-2,-2:2,1:nw,sps2) = wAdj1(-2:2,-2:2,1:nw)
        pAdj(-2:2,-2,-2:2,sps2)      = pAdj1(-2:2,-2:2)
     endif

     !===========================================================

  case (jMax)

     if( secondHalo ) then
        wAdj(-2:2, 2,-2:2,1:nw,sps2) = wAdj0(-2:2,-2:2,1:nw)
        wAdj(-2:2, 1,-2:2,1:nw,sps2) = wAdj1(-2:2,-2:2,1:nw)

        pAdj(-2:2, 2,-2:2,sps2) = pAdj0(-2:2,-2:2)
        pAdj(-2:2, 1,-2:2,sps2) = pAdj1(-2:2,-2:2)
     else
        wAdj(-2:2, 2,-2:2,1:nw,sps2) = wAdj1(-2:2,-2:2,1:nw)
        pAdj(-2:2, 2,-2:2,sps2)      = pAdj1(-2:2,-2:2)
     endif

     !===========================================================

  case (kMin)

     if( secondHalo ) then
        wAdj(-2:2,-2:2,-2,1:nw,sps2) = wAdj0(-2:2,-2:2,1:nw)
        wAdj(-2:2,-2:2,-1,1:nw,sps2) = wAdj1(-2:2,-2:2,1:nw)

        pAdj(-2:2,-2:2,-2,sps2) = pAdj0(-2:2,-2:2)
        pAdj(-2:2,-2:2,-1,sps2) = pAdj1(-2:2,-2:2)
     else
        wAdj(-2:2,-2:2,-2,1:nw,sps2) = wAdj1(-2:2,-2:2,1:nw)
        pAdj(-2:2,-2:2,-2,sps2)      = pAdj1(-2:2,-2:2)
     endif

     !===========================================================

  case (kMax)

     if( secondHalo ) then
        wAdj(-2:2,-2:2, 2,1:nw,sps2) = wAdj0(-2:2,-2:2,1:nw)
        wAdj(-2:2,-2:2, 1,1:nw,sps2) = wAdj1(-2:2,-2:2,1:nw)

        pAdj(-2:2,-2:2, 2,sps2) = pAdj0(-2:2,-2:2)
        pAdj(-2:2,-2:2, 1,sps2) = pAdj1(-2:2,-2:2)
     else
        wAdj(-2:2,-2:2, 2,1:nw,sps2) = wAdj1(-2:2,-2:2,1:nw)
        pAdj(-2:2,-2:2, 2,sps2)      = pAdj1(-2:2,-2:2)
     endif

  end select

end subroutine replaceBCStatesNKPC
