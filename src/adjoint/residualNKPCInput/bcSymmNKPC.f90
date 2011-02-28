!
!      ******************************************************************
!      *                                                                *
!      * File:          bcSymmAdj.f90                                   *
!      * Author:        Edwin van der Weide,C.A.(Sandy) Mader           *
!      * Starting date: 04-17-2008                                      *
!      * Last modified: 04-17-2008                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine bcSymmNKPC(wAdj,pAdj,normAdj,iCell,jCell,kCell,secondHalo,nnn,level,sps,sps2)
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
  use inputTimeSpectral !nIntervalTimespectral
  implicit none
  !
  !      Subroutine arguments.
  !
  logical:: secondHalo
  integer(kind=intType)::nnn,level,sps,sps2
  real(kind=realType), dimension(-2:2,-2:2,-2:2,nw,nTimeIntervalsSpectral), &
       intent(in) :: wAdj
  real(kind=realType), dimension(-2:2,-2:2,-2:2,nTimeIntervalsSpectral),intent(in) :: pAdj
  real(kind=realType), dimension(nBocos,-2:2,-2:2,3,nTimeIntervalsSpectral), intent(in) :: normAdj

  !
  !      Local variables.
  !
  integer(kind=intType) :: kk, mm, nn, i, j, l,ii,jj

  real(kind=realType) :: vn, nnx, nny, nnz

  !real(kind=realType), dimension(:,:),   pointer :: gamma1, gamma2
  real(kind=realType), dimension(-2:2,-2:2) :: gamma1,gamma2

  real(kind=realType), dimension(-2:2,-2:2,nw) :: wAdj0, wAdj1
  real(kind=realType), dimension(-2:2,-2:2,nw) :: wAdj2, wAdj3
  real(kind=realType), dimension(-2:2,-2:2)    :: pAdj0, pAdj1
  real(kind=realType), dimension(-2:2,-2:2)    :: pAdj2, pAdj3

  real(kind=realType), dimension(-2:2,-2:2,-2:2)::rlvAdj, revAdj
  real(kind=realType), dimension(-2:2,-2:2)::rlvAdj1, rlvAdj2
  real(kind=realType), dimension(-2:2,-2:2)::revAdj1, revAdj2

  integer(kind=intType) ::iCell, jCell,kCell
  integer(kind=intType) ::isbeg,jsbeg,ksbeg,isend,jsend,ksend
  integer(kind=intType) ::ibbeg,jbbeg,kbbeg,ibend,jbend,kbend
  integer(kind=intType) ::icbeg,jcbeg,kcbeg,icend,jcend,kcend
  integer(kind=intType) :: iOffset, jOffset, kOffset
  logical :: computeBC

  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !

  ! Loop over the boundary condition subfaces of this block.

  bocos: do nn=1,nBocos

     call checkOverlapNKPC(nn,icell,jcell,kcell,isbeg,jsbeg,&
          ksbeg,isend,jsend,ksend,ibbeg,jbbeg,kbbeg,ibend,jbend,kbend,&
          computeBC)

     if (computeBC) then

        ! Check for symmetry boundary condition.

        symmetry: if(BCType(nn) == symm) then

           !Copy the states and other parameters to subfaces
           call extractBCStatesNKPC(nn,wAdj,pAdj,wAdj0, wAdj1, wAdj2,wAdj3,&
                pAdj0,pAdj1, pAdj2,pAdj3,&
                rlvAdj, revAdj,rlvAdj1, rlvAdj2,revAdj1, revAdj2,iOffset,&
                jOffset, kOffset,iCell, jCell,kCell,&
                isbeg,jsbeg,ksbeg,isend,jsend,ksend,ibbeg,jbbeg,kbbeg,ibend,&
                jbend,kbend,icbeg,jcbeg,icend,jcend,secondHalo,nnn,level,sps,sps2)

           select case (BCFaceID(nn))
           case (iMin)
              gamma1 = gamma(1, jcell-2:jcell+2,kcell-2:kcell+2)
              gamma2 = gamma(2, jcell-2:jcell+2,kcell-2:kcell+2)
           case (iMax)
              gamma1 = gamma(ie,jcell-2:jcell+2,kcell-2:kcell+2)
              gamma2 = gamma(il,jcell-2:jcell+2,kcell-2:kcell+2)
           case (jMin)
              gamma1 = gamma(icell-2:icell+2,1, kcell-2:kcell+2)
              gamma2 = gamma(icell-2:icell+2,2, kcell-2:kcell+2)
           case (jMax)
              gamma1 = gamma(icell-2:icell+2,je,kcell-2:kcell+2)
              gamma2 = gamma(icell-2:icell+2,jl,kcell-2:kcell+2)
           case (kMin)
              gamma1 = gamma(icell-2:icell+2,jcell-2:jcell+2,1 )
              gamma2 = gamma(icell-2:icell+2,jcell-2:jcell+2,2 )
           case (kMax)
              gamma1 = gamma(icell-2:icell+2,jcell-2:jcell+2,ke)
              gamma2 = gamma(icell-2:icell+2,jcell-2:jcell+2,kl)
           end select

           ! Loop over the generic subface to set the state in the
           ! halo cells.

           do j=jcBeg,jcEnd
              do i=icBeg,icEnd
                 ii = i - iOffset
                 jj = j - jOffset

                 ! Store the three components of the unit normal a
                 ! bit easier.

                 nnx = normAdj(nn,ii,jj,1,sps2)!BCData(nn)%norm(i,j,1)
                 nny = normAdj(nn,ii,jj,2,sps2)!BCData(nn)%norm(i,j,2)
                 nnz = normAdj(nn,ii,jj,3,sps2)!BCData(nn)%norm(i,j,3)

                 ! Determine twice the normal velocity component,
                 ! which must be substracted from the donor velocity
                 ! to obtain the halo velocity.

                 vn = two*(wAdj2(ii,jj,ivx)*nnx + wAdj2(ii,jj,ivy)*nny &
                      +      wAdj2(ii,jj,ivz)*nnz)

                 ! Determine the flow variables in the halo cell.

                 wAdj1(ii,jj,irho) = wAdj2(ii,jj,irho)

                 wAdj1(ii,jj,ivx) = wAdj2(ii,jj,ivx) - vn*nnx
                 wAdj1(ii,jj,ivy) = wAdj2(ii,jj,ivy) - vn*nny
                 wAdj1(ii,jj,ivz) = wAdj2(ii,jj,ivz) - vn*nnz

                 wAdj1(ii,jj,irhoE) = wAdj2(ii,jj,irhoE)

                 ! Simply copy the turbulent variables.

                 do l=nt1MG,nt2MG
                    wAdj1(ii,jj,l) = wAdj2(ii,jj,l)
                 enddo

                 ! Set the pressure and gamma and possibly the
                 ! laminar and eddy viscosity in the halo.

                 gamma1(ii,jj) = gamma2(ii,jj)
                 pAdj1(ii,jj)    = pAdj2(ii,jj)
                 if( viscous )   rlvAdj1(ii,jj) = rlvAdj2(ii,jj)
                 if( eddyModel ) revAdj1(ii,jj) = revAdj2(ii,jj)

                 if (secondHalo) then
                    ! Determine twice the normal velocity component,
                    ! which must be substracted from the donor velocity
                    ! to obtain the halo velocity.

                    vn = two*(wAdj3(ii,jj,ivx)*nnx + wAdj3(ii,jj,ivy)*nny &
                         +      wAdj3(ii,jj,ivz)*nnz)

                    ! Determine the flow variables in the halo cell.

                    wAdj0(ii,jj,irho) = wAdj3(ii,jj,irho)

                    wAdj0(ii,jj,ivx) = wAdj3(ii,jj,ivx) - vn*nnx
                    wAdj0(ii,jj,ivy) = wAdj3(ii,jj,ivy) - vn*nny
                    wAdj0(ii,jj,ivz) = wAdj3(ii,jj,ivz) - vn*nnz

                    wAdj0(ii,jj,irhoE) = wAdj3(ii,jj,irhoE)

                    ! Simply copy the turbulent variables.

                    do l=nt1MG,nt2MG
                       wAdj1(ii,jj,l) = wAdj3(ii,jj,l)
                    enddo

                    ! Set the pressure and gamma and possibly the
                    ! laminar and eddy viscosity in the halo.

                    gamma1(ii,jj) = gamma2(ii,jj)
                    pAdj0(ii,jj)    = pAdj3(ii,jj)
                 endif

              enddo
           enddo

           call replaceBCStatesNKPC(nn,  wAdj0,wAdj1, wAdj2, wAdj3,&
                pAdj0,pAdj1, pAdj2, pAdj3,rlvAdj1, rlvAdj2,revAdj1, revAdj2,&
                iCell, jCell,kCell,&
                wAdj,pAdj,rlvAdj,revAdj,secondHalo,nnn,level,sps,sps2)
        endif symmetry
     endif
  enddo bocos
end subroutine bcSymmNKPC
