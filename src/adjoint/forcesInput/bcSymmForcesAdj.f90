!
!      ******************************************************************
!      *                                                                *
!      * File:          bcSymmAdj.f90                                   *
!      * Author:        Edwin van der Weide,C.A.(Sandy) Mader           *
!      * Starting date: 04-17-2008                                      *
!      * Last modified: 06-12-2008                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine bcSymmForcesAdj(secondHalo, wAdj,pAdj,normAdj,nn,&
     iiBeg,iiEnd,jjBeg,jjEnd,i2Beg,i2End,j2Beg,j2End)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * bcSymmAdj applies the symmetry boundary conditions to a single *
  !      * cell stencil.
  !      * It is assumed that the pointers in blockPointers are already   *
  !      * set to the correct block on the correct grid level.            *
  !      *                                                                *
  !      * In case also the second halo must be set the loop over the     *
  !      * boundary subfaces is executed twice. This is the only correct  *
  !      * way in case the block contains only 1 cell between two         *
  !      * symmetry planes, i.e. a 2D problem.                            *
  !      *                                                                *
  !      ******************************************************************
  !
  use blockPointers, only : ie, ib,il, je, jb,jl, ke, kb,kl, nBocos, &
       gamma,BCFaceID, BCType, BCData
  use BCTypes
  use constants
  use flowVarRefState  !nw
  use iteration        !nt1mg,nt2mg
  implicit none
  !
  !      Subroutine arguments.
  !
  logical, intent(in) :: secondHalo
  integer(kind=intType), intent(in) :: iiBeg,iiEnd,jjBeg,jjEnd
  integer(kind=intType), intent(in) :: i2Beg,i2End,j2Beg,j2End
  
  !       integer(kind=intType) :: iCell, jCell, kCell
  integer(kind=intType):: icBeg,icEnd,jcBeg,jcEnd
  
  real(kind=realType), dimension(0:ib,0:jb,0:kb,nw), intent(in) :: wAdj
  real(kind=realType), dimension(0:ib,0:jb,0:kb), intent(in) :: pAdj
  
  real(kind=realType), dimension(iiBeg:iiEnd,jjBeg:jjEnd,3) :: normAdj

  real(kind=realType), dimension(iibeg:iiend,jjbeg:jjend,nw):: wAdj0, wAdj1
  real(kind=realType), dimension(iibeg:iiend,jjbeg:jjend,nw):: wAdj2, wAdj3
  real(kind=realType), dimension(iibeg:iiend,jjbeg:jjend) :: pAdj0, pAdj1
  real(kind=realType), dimension(iibeg:iiend,jjbeg:jjend) :: pAdj2, pAdj3
  
  real(kind=realType), dimension(0:ib,0:jb,0:kb)::rlvAdj, revAdj
  real(kind=realType), dimension(iibeg:iiend,jjbeg:jjend)::rlvAdj1, rlvAdj2
  real(kind=realType), dimension(iibeg:iiend,jjbeg:jjend)::revAdj1, revAdj2


!!$  real(kind=realType), dimension(-2:2,-2:2,-2:2,nw), &
!!$       intent(in) :: wAdj
!!$  real(kind=realType), dimension(-2:2,-2:2,-2:2),intent(in) :: pAdj
!!$  real(kind=realType), dimension(nBocos,-2:2,-2:2,3), intent(in) :: normAdj

  !
  !      Local variables.
  !
  integer(kind=intType) :: nn, i, j, l

  real(kind=realType) :: vn, nnx, nny, nnz

  real(kind=realType), dimension(iibeg:iiend,jjbeg:jjend):: gamma1, gamma2

!!$  real(kind=realType), dimension(-2:2,-2:2,nw) :: wAdj0, wAdj1
!!$  real(kind=realType), dimension(-2:2,-2:2,nw) :: wAdj2, wAdj3
!!$  real(kind=realType), dimension(-2:2,-2:2)    :: pAdj0, pAdj1
!!$  real(kind=realType), dimension(-2:2,-2:2)    :: pAdj2, pAdj3
!!$
!!$  real(kind=realType), dimension(-2:2,-2:2,-2:2)::rlvAdj, revAdj
!!$  real(kind=realType), dimension(-2:2,-2:2)::rlvAdj1, rlvAdj2
!!$  real(kind=realType), dimension(-2:2,-2:2)::revAdj1, revAdj2
!!$
!!$  integer(kind=intType) ::iCell, jCell,kCell
!!$  integer(kind=intType) ::isbeg,jsbeg,ksbeg,isend,jsend,ksend
!!$  integer(kind=intType) ::ibbeg,jbbeg,kbbeg,ibend,jbend,kbend
!!$  integer(kind=intType) ::icbeg,jcbeg,kcbeg,icend,jcend,kcend
!!$  integer(kind=intType) :: iOffset, jOffset, kOffset
!!$  logical :: computeBC
!!$!
!!$!      Interfaces
!!$!
!!$       interface
!!$         subroutine setBcPointers(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
!!$                                  rev1, rev2, offset)
!!$           use blockPointers
!!$           implicit none
!!$
!!$           integer(kind=intType), intent(in) :: nn, offset
!!$           real(kind=realType), dimension(:,:,:), pointer :: ww1, ww2
!!$           real(kind=realType), dimension(:,:),   pointer :: pp1, pp2
!!$           real(kind=realType), dimension(:,:),   pointer :: rlv1, rlv2
!!$           real(kind=realType), dimension(:,:),   pointer :: rev1, rev2
!!$         end subroutine setBcPointers
!!$       end interface
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
!!$  ! Set the value of kk; kk == 0 means only single halo, kk == 1
!!$  ! double halo.
!!$
!!$  kk = 0
!!$  if( secondHalo ) kk = 1
!!$
!!$  ! Loop over the number of times the halo computation must be done.
!!$
!!$  nHalo: do mm=0,kk

     ! Loop over the boundary condition subfaces of this block.

!!$     bocos: do nn=1,nBocos
!!$
!!$        call checkOverlapAdj(nn,icell,jcell,kcell,isbeg,jsbeg,&
!!$            ksbeg,isend,jsend,ksend,ibbeg,jbbeg,kbbeg,ibend,jbend,kbend,&
!!$            computeBC)
!!$        !(nn,icell,jcell,kcell,isbeg,jsbeg,&
!!$        !     ksbeg,isend,jsend,ksend,ibbeg,jbbeg,kbbeg,ibend,jbend,kbend,&
!!$        !     computeBC)
!!$
!!$
!!$        if (computeBC) then

           ! Check for symmetry boundary condition.

  symmetry: if(BCType(nn) == symm) then

     icbeg = iibeg
     icend = iiend
     jcbeg = jjbeg
     jcend = jjend

!!$             ! Nullify the pointers, because some compilers require that.
!!$
!!$             nullify(ww1, ww2, pp1, pp2, rlv1, rlv2, rev1, rev2)

!!$              ! Set the pointers to the correct subface.
!!$              call setBcPointers(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
!!$                                rev1, rev2, mm)

              !Copy the states and other parameters to subfaces
     call extractBCStatesForcesAdj(nn,wAdj,pAdj,wAdj0, wAdj1, wAdj2,wAdj3,&
          pAdj0,pAdj1, pAdj2,pAdj3,&
          rlvAdj, revAdj,rlvAdj1, rlvAdj2,revAdj1, revAdj2,&
          iibeg,jjbeg,iiend,jjend,secondHalo)

     ! Set the additional pointers for gamma1 and gamma2.
     
     select case (BCFaceID(nn))
     case (iMin)
        gamma1 = gamma(1,iibeg:iiend,jjbeg:jjend)
        gamma2 = gamma(2,iibeg:iiend,jjbeg:jjend)
     case (iMax)
        gamma1 = gamma(ie,iibeg:iiend,jjbeg:jjend)
        gamma2 = gamma(il,iibeg:iiend,jjbeg:jjend)
     case (jMin)
        gamma1 = gamma(iibeg:iiend,1, jjbeg:jjend)
        gamma2 = gamma(iibeg:iiend,2, jjbeg:jjend)
     case (jMax)
        gamma1 = gamma(iibeg:iiend,je,jjbeg:jjend)
        gamma2 = gamma(iibeg:iiend,jl,jjbeg:jjend)
     case (kMin)
        gamma1 = gamma(iibeg:iiend,jjbeg:jjend,1 )
        gamma2 = gamma(iibeg:iiend,jjbeg:jjend,2 )
     case (kMax)
        gamma1 = gamma(iibeg:iiend,jjbeg:jjend,ke)
        gamma2 = gamma(iibeg:iiend,jjbeg:jjend,kl)
     end select

              ! Loop over the generic subface to set the state in the
              ! halo cells.

!!$             do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
!!$               do i=BCData(nn)%icBeg, BCData(nn)%icEnd
     do j=jcBeg,jcEnd
        do i=icBeg,icEnd

           ! Store the three components of the unit normal a
           ! bit easier.
           
           nnx = normAdj(i,j,1)!BCData(nn)%norm(i,j,1)
           nny = normAdj(i,j,2)!BCData(nn)%norm(i,j,2)
           nnz = normAdj(i,j,3)!BCData(nn)%norm(i,j,3)
           
           ! Determine twice the normal velocity component,
           ! which must be substracted from the donor velocity
           ! to obtain the halo velocity.
           
           vn = two*(wAdj2(i,j,ivx)*nnx + wAdj2(i,j,ivy)*nny &
                +      wAdj2(i,j,ivz)*nnz)
           
           ! Determine the flow variables in the halo cell.
           
           wAdj1(i,j,irho) = wAdj2(i,j,irho)
           
           wAdj1(i,j,ivx) = wAdj2(i,j,ivx) - vn*nnx
           wAdj1(i,j,ivy) = wAdj2(i,j,ivy) - vn*nny
           wAdj1(i,j,ivz) = wAdj2(i,j,ivz) - vn*nnz
           
           wAdj1(i,j,irhoE) = wAdj2(i,j,irhoE)
           
           ! Simply copy the turbulent variables.
           
           do l=nt1MG,nt2MG
              wAdj1(i,j,l) = wAdj2(i,j,l)
           enddo
           
           ! Set the pressure and gamma and possibly the
           ! laminar and eddy viscosity in the halo.
           
           gamma1(i,j) = gamma2(i,j)
           pAdj1(i,j)    = pAdj2(i,j)
           if( viscous )   rlvAdj1(i,j) = rlvAdj2(i,j)
           if( eddyModel ) revAdj1(i,j) = revAdj2(i,j)
           
           if(secondHalo)then
              vn = two*(wAdj3(i,j,ivx)*nnx + wAdj3(i,j,ivy)*nny &
                   +      wAdj3(i,j,ivz)*nnz)
              
              ! Determine the flow variables in the halo cell.
              
              wAdj0(i,j,irho) = wAdj3(i,j,irho)
              
              wAdj0(i,j,ivx) = wAdj3(i,j,ivx) - vn*nnx
              wAdj0(i,j,ivy) = wAdj3(i,j,ivy) - vn*nny
              wAdj0(i,j,ivz) = wAdj3(i,j,ivz) - vn*nnz
              
              wAdj0(i,j,irhoE) = wAdj3(i,j,irhoE)
              
              ! Simply copy the turbulent variables.
              
              do l=nt1MG,nt2MG
                 wAdj0(i,j,l) = wAdj3(i,j,l)
              enddo
              
              ! Set the pressure and gamma and possibly the
              ! laminar and eddy viscosity in the halo.
              
              gamma1(i,j) = gamma2(i,j)
              pAdj0(i,j)    = pAdj3(i,j)
!!$              if( viscous )   rlvAdj1(i,j) = rlvAdj2(i,j)
!!$              if( eddyModel ) revAdj1(i,j) = revAdj2(i,j)
           endif
        enddo
     enddo
     
     call replaceBCStatesForcesAdj(nn,  wAdj0,wAdj1, wAdj2, wAdj3,&
                pAdj0,pAdj1, pAdj2, pAdj3,rlvAdj1, rlvAdj2,revAdj1, revAdj2,&
                wAdj,pAdj,rlvAdj,revAdj,&
                iibeg,jjbeg,iiend,jjend,secondHalo)
         
  endif symmetry

end subroutine bcSymmForcesAdj
