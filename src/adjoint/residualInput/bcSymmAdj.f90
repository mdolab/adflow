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
subroutine bcSymmAdj(wAdj,pAdj,normAdj,iCell,jCell,kCell,secondHalo)
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
  real(kind=realType), dimension(-2:2,-2:2,-2:2,nw), &
       intent(in) :: wAdj
  real(kind=realType), dimension(-2:2,-2:2,-2:2),intent(in) :: pAdj
  real(kind=realType), dimension(nBocos,-2:2,-2:2,3), intent(in) :: normAdj

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
  !print *,'in bcSymmAdj'

!Accounted for inside bocos loop
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

     bocos: do nn=1,nBocos

        call checkOverlapAdj(nn,icell,jcell,kcell,isbeg,jsbeg,&
            ksbeg,isend,jsend,ksend,ibbeg,jbbeg,kbbeg,ibend,jbend,kbend,&
            computeBC)
        !(nn,icell,jcell,kcell,isbeg,jsbeg,&
        !     ksbeg,isend,jsend,ksend,ibbeg,jbbeg,kbbeg,ibend,jbend,kbend,&
        !     computeBC)


        if (computeBC) then

           ! Check for symmetry boundary condition.

           symmetry: if(BCType(nn) == symm) then

!!$             ! Nullify the pointers, because some compilers require that.
!!$
!!$             nullify(ww1, ww2, pp1, pp2, rlv1, rlv2, rev1, rev2)

!!$              ! Set the pointers to the correct subface.
!!$              call setBcPointers(nn, ww1, ww2, pp1, pp2, rlv1, rlv2, &
!!$                                rev1, rev2, mm)

              !Copy the states and other parameters to subfaces
              call extractBCStatesAdj(nn,wAdj,pAdj,wAdj0, wAdj1, wAdj2,wAdj3,&
                   pAdj0,pAdj1, pAdj2,pAdj3,&
                   rlvAdj, revAdj,rlvAdj1, rlvAdj2,revAdj1, revAdj2,iOffset,&
                   jOffset, kOffset,iCell, jCell,kCell,&
                   isbeg,jsbeg,ksbeg,isend,jsend,ksend,ibbeg,jbbeg,kbbeg,ibend,&
                   jbend,kbend,icbeg,jcbeg,icend,jcend,secondHalo)
              !(nn,wAdj,pAdj, wAdj1, wAdj2, pAdj1, pAdj2,&
              !     rlvAdj, revAdj,rlvAdj1, rlvAdj2,revAdj1, revAdj2,iOffset,&
              !     jOffset, kOffset,iCell, jCell,kCell,&
              !     isbeg,jsbeg,ksbeg,isend,jsend,ksend,ibbeg,jbbeg,kbbeg,ibend,&
              !     jbend,kbend,icbeg,jcbeg,icend,jcend)



              ! Set the additional pointers for gamma1 and gamma2.

!              select case (BCFaceID(nn))
!              case (iMin)
!                 gamma1 => gamma(1, 1:,1:); gamma2 => gamma(2, 1:,1:)
!              case (iMax)
!                 gamma1 => gamma(ie,1:,1:); gamma2 => gamma(il,1:,1:)
!              case (jMin)
!                 gamma1 => gamma(1:,1, 1:); gamma2 => gamma(1:,2, 1:)
!              case (jMax)
!                 gamma1 => gamma(1:,je,1:); gamma2 => gamma(1:,jl,1:)
!              case (kMin)
!                 gamma1 => gamma(1:,1:,1 ); gamma2 => gamma(1:,1:,2 )
!              case (kMax)
!                 gamma1 => gamma(1:,1:,ke); gamma2 => gamma(1:,1:,kl)
!              end select

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

!!$             do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
!!$               do i=BCData(nn)%icBeg, BCData(nn)%icEnd
              do j=jcBeg,jcEnd
                 do i=icBeg,icEnd
                    ii = i - iOffset
                    jj = j - jOffset

                    !print *,'i,j',i,j,jcBeg,jcEnd,ii,jj
                    ! Store the three components of the unit normal a
                    ! bit easier.

                    nnx = normAdj(nn,ii,jj,1)!BCData(nn)%norm(i,j,1)
                    nny = normAdj(nn,ii,jj,2)!BCData(nn)%norm(i,j,2)
                    nnz = normAdj(nn,ii,jj,3)!BCData(nn)%norm(i,j,3)

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
!!$                       if( viscous )   rlvAdj0(ii,jj) = rlvAdj3(ii,jj)
!!$                       if( eddyModel ) revAdj0(ii,jj) = revAdj3(ii,jj)
                    endif
                    
                 enddo
              enddo

              call replaceBCStatesAdj(nn,  wAdj0,wAdj1, wAdj2, wAdj3,&
                   pAdj0,pAdj1, pAdj2, pAdj3,rlvAdj1, rlvAdj2,revAdj1, revAdj2,&
                   iCell, jCell,kCell,&
                   wAdj,pAdj,rlvAdj,revAdj,secondHalo)

           endif symmetry
        endif
     enddo bocos
!  enddo nHalo

end subroutine bcSymmAdj
