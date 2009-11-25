!
!      ******************************************************************
!      *                                                                *
!      * File:          bcFarfieldAdj.f90                               *
!      * Author:        Edwin van der Weide                             *
!      *                Seongim Choi,C.A.(Sandy) Mader                  *
!      * Starting date: 03-21-2006                                      *
!      * Last modified: 04-23-2008                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine bcFarfieldAdj(secondHalo,wInfAdj,pInfCorrAdj, wAdj,pAdj,      &
     siAdj, sjAdj, skAdj, normAdj,rFaceAdj,iCell,jCell,kCell,nn,level,sps,sps2)

  !
  !      ******************************************************************
  !      *                                                                *
  !      * bcFarfieldAdj applies the farfield boundary condition to       *
  !      * subface nn of the block to which the pointers in blockPointers *
  !      * currently point.                                               *
  !      *                                                                *
  !      ******************************************************************
  !
  use blockPointers, only : BCData,nBocos,BCType,bcfaceid,gamma,il,jl,kl,w,p
  use constants         ! irho,ivx,ivy,ivz
  use flowVarRefState   ! gammaInf, wInf, pInfCorr
  use BCTypes
  use iteration
  use inputTimeSpectral !nIntervalTimespectral
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType) ::nn,level,sps,sps2! it's not needed anymore w/ normAdj
  !       integer(kind=intType), intent(in) :: icBeg, icEnd, jcBeg, jcEnd
  !       integer(kind=intType), intent(in) :: iOffset, jOffset

  integer(kind=intType) ::iCell, jCell,kCell
  integer(kind=intType) ::isbeg,jsbeg,ksbeg,isend,jsend,ksend
  integer(kind=intType) ::ibbeg,jbbeg,kbbeg,ibend,jbend,kbend
  integer(kind=intType) ::icbeg,jcbeg,kcbeg,icend,jcend,kcend
  integer(kind=intType) :: iOffset, jOffset, kOffset

  real(kind=realType), dimension(-2:2,-2:2,-2:2)::rlvAdj, revAdj
  real(kind=realType), dimension(-2:2,-2:2)::rlvAdj1, rlvAdj2
  real(kind=realType), dimension(-2:2,-2:2)::revAdj1, revAdj2

  real(kind=realType), dimension(nw),intent(in)::wInfAdj


!  real(kind=realType), dimension(-2:2,-2:2,-2:2,3), intent(in) :: &
!       siAdj, sjAdj, skAdj
  real(kind=realType), dimension(-3:2,-3:2,-3:2,3,nTimeIntervalsSpectral), intent(in) :: &
       siAdj, sjAdj, skAdj
  real(kind=realType), dimension(nBocos,-2:2,-2:2,3,nTimeIntervalsSpectral), intent(in) :: normAdj
  real(kind=realType), dimension(nBocos,-2:2,-2:2,nTimeIntervalsSpectral), intent(in) :: rFaceAdj
  real(kind=realType), dimension(-2:2,-2:2,-2:2,nw,nTimeIntervalsSpectral), &
       intent(in) :: wAdj
  real(kind=realType), dimension(-2:2,-2:2,-2:2,nTimeIntervalsSpectral),intent(in) :: pAdj
  real(kind=realType)::pInfCorrAdj

  logical, intent(in) :: secondHalo

  real(kind=realType), dimension(-2:2,-2:2,nw) :: wAdj0, wAdj1
  real(kind=realType), dimension(-2:2,-2:2,nw) :: wAdj2, wAdj3
  real(kind=realType), dimension(-2:2,-2:2)    :: pAdj0, pAdj1
  real(kind=realType), dimension(-2:2,-2:2)    :: pAdj2, pAdj3

  !real(kind=realType), dimension(nBocos,-2:2,-2:2,3), intent(in) :: normAdj
  !
  !      Local variables.
  !
  integer(kind=intType) :: i, j, l, ii, jj,nnbcs

  real(kind=realType) :: nnx, nny, nnz
  real(kind=realType) :: gm1, ovgm1, gm53, factK, ac1, ac2
  real(kind=realType) :: r0, u0, v0, w0, qn0, vn0, c0, s0
  real(kind=realType) :: re, ue, ve, we, qne, ce
  real(kind=realType) :: qnf, cf, uf, vf, wf, sf, cc, qq
  real(kind=realType) :: rface



  logical :: computeBC

!!$!File Parameters remove for AD
!!$      integer :: unitxAD = 15,ierror
!!$      integer ::iii,iiii,jjj,jjjj,kkk,kkkk,nnnn,istart2,jstart2,kstart2,iend2,jend2,kend2,n
!!$      character(len = 16)::outfile
!!$      
!!$      outfile = "xAD.txt"
!!$      
!!$      open (UNIT=unitxAD,File=outfile,status='old',position='append',action='write',iostat=ierror)
!!$      if(ierror /= 0)                        &
!!$           call terminate("verifyResiduals", &
!!$           "Something wrong when &
!!$           &calling open")

!
!      Interfaces
!



  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Some constants needed to compute the riemann invariants.

  gm1   = gammaInf - one
  ovgm1 = one/gm1
  gm53  =  gammaInf - five*third
  factK = -ovgm1*gm53

  ! Compute the three velocity components, the speed of sound and
  ! the entropy of the free stream.

  r0  = one/wInfAdj(irho)
  u0  = wInfAdj(ivx)
  v0  = wInfAdj(ivy)
  w0  = wInfAdj(ivz)
  c0  = sqrt(gammaInf*pInfCorrAdj*r0)
  s0  = wInfAdj(irho)**gammaInf/pInfCorrAdj

  ! Loop over the boundary condition subfaces of this block.

  bocos: do nnbcs=1,nBocos
     !print *,'nn',nn,nnbcs
     call checkOverlapAdj(nnbcs,icell,jcell,kcell,isbeg,jsbeg,&
          ksbeg,isend,jsend,ksend,ibbeg,jbbeg,kbbeg,ibend,jbend,kbend,&
          computeBC)


     if (computeBC) then

        ! Check for farfield boundary conditions.
        

        testFarfield: if(BCType(nnbcs) == FarField) then

           
           call extractBCStatesAdj(nnbcs,wAdj,pAdj,wAdj0, wAdj1, wAdj2,wAdj3,&
            pAdj0,pAdj1, pAdj2,pAdj3,&
            rlvAdj, revAdj,rlvAdj1, rlvAdj2,revAdj1, revAdj2,iOffset,&
            jOffset, kOffset,iCell, jCell,kCell,&
            isbeg,jsbeg,ksbeg,isend,jsend,ksend,ibbeg,jbbeg,kbbeg,ibend,&
            jbend,kbend,icbeg,jcbeg,icend,jcend,secondHalo,nn,level,sps,sps2)
           !print *,'extract',secondhalo,icell,jcell,kcell

           ! Loop over the generic subface to set the state in the
           ! halo cells.

           do j=jcBeg, jcEnd
              do i=icBeg, icEnd
          
                 ii = i - iOffset
                 jj = j - jOffset
             
                 rface = rFaceAdj(nnbcs,ii,jj,sps2)!BCData(nn)%rface(i,j)
                 
                 ! Store the three components of the unit normal a
                 ! bit easier.

                 nnx = normAdj(nnbcs,ii,jj,1,sps2)
                 nny = normAdj(nnbcs,ii,jj,2,sps2)
                 nnz = normAdj(nnbcs,ii,jj,3,sps2)
!!$                 if(i==1 .and. j ==1)then
!!$                    print *,'normxAdj',nnx,nnbcs!,i,j,ii,jj
!!$                    print *,'normyAdj',nny,nnbcs
!!$                    print *,'normzAdj',nnz,nnbcs
!!$                    print *,'rfaceAdj',rface,nnbcs
!!$                 endif

!!$                 write(unitxAD,11)int((icbeg+icend)/2),int((jcbeg+jcend)/2), i,j,nnbcs,normAdj(nnbcs,ii,jj,1,sps2), rFaceAdj(nnbcs,ii,jj,sps2) 
!!$11               format(1x,'wadj',5I8,2f20.14)       
                 ! Compute the normal velocity of the free stream and
                 ! substract the normal velocity of the mesh.

                 qn0 = u0*nnx + v0*nny + w0*nnz
                 vn0 = qn0 - rface

                 ! Compute the three velocity components, the normal
                 ! velocity and the speed of sound of the current state
                 ! in the internal cell.

                 re  = one/wAdj2(ii,jj,irho)
                 ue  = wAdj2(ii,jj,ivx)
                 ve  = wAdj2(ii,jj,ivy)
                 we  = wAdj2(ii,jj,ivz)
                 qne = ue*nnx + ve*nny + we*nnz
                 ce  = sqrt(gammaInf*pAdj2(ii,jj)*re)

                 ! Compute the new values of the riemann invariants in
                 ! the halo cell. Either the value in the internal cell
                 ! is taken (positive sign of the corresponding
                 ! eigenvalue) or the free stream value is taken
                 ! (otherwise).

                 if(vn0 > -c0) then       ! Outflow or subsonic inflow.
                    ac1 = qne + two*ovgm1*ce
                 else                     ! Supersonic inflow.
                    ac1 = qn0 + two*ovgm1*c0
                 endif

                 if(vn0 > c0) then        ! Supersonic outflow.
                    ac2 = qne - two*ovgm1*ce
                 else                     ! Inflow or subsonic outflow.
                    ac2 = qn0 - two*ovgm1*c0
                 endif

                 qnf = half*  (ac1 + ac2)
                 cf  = fourth*(ac1 - ac2)*gm1
                
                 if(vn0 > zero) then            ! Outflow.

                    uf = ue + (qnf - qne)*nnx
                    vf = ve + (qnf - qne)*nny
                    wf = we + (qnf - qne)*nnz
                    sf = wAdj2(ii,jj,irho)**gammaInf/pAdj2(ii,jj)


                    do l=nt1MG,nt2MG
                       wAdj1(ii,jj,l) = wAdj2(ii,jj,l)
                    enddo

                 else                           ! Inflow

                    uf = u0 + (qnf - qn0)*nnx
                    vf = v0 + (qnf - qn0)*nny
                    wf = w0 + (qnf - qn0)*nnz
                    sf = s0
                  
                    do l=nt1MG,nt2MG
                       wAdj1(ii,jj,l) = wInfAdj(l)
                    enddo

                 endif

                 ! Compute the density, velocity and pressure in the
                 ! halo cell.

                 cc = cf*cf/gammaInf
                 qq = uf*uf + vf*vf + wf*wf
               
                 wAdj1(ii,jj,irho) = (sf*cc)**ovgm1
                 wAdj1(ii,jj,ivx)  = uf
                 wAdj1(ii,jj,ivy)  = vf
                 wAdj1(ii,jj,ivz)  = wf
                 pAdj1(ii,jj)      = wAdj1(ii,jj,irho)*cc
                 
                 !print*,'wadj1',i,j,wAdj1(ii,jj,irho),wAdj1(ii,jj,ivx),wAdj1(ii,jj,ivy),wAdj1(ii,jj,ivz),wAdj1(ii,jj,irhoE),pAdj1(ii,jj) 
                 
                 ! Compute the total energy.

                 wAdj1(ii,jj,irhoE) = ovgm1*pAdj1(ii,jj)     &
                      + half*wAdj1(ii,jj,irho) &
                      *      (uf**2 + vf**2 + wf**2)

                 if( kPresent )                            &
                      wAdj1(ii,jj,irhoE) = wAdj1(ii,jj,irhoE) &
                      - factK*wAdj1(ii,jj,irho) &
                      *       wAdj1(ii,jj,itu1)

                 !
                 !        Input the viscous effects - rlv1(), and rev1()
                 !

              enddo
           enddo

           ! Extrapolate the state vectors in case a second halo
           ! is needed.
          ! print *,'second halo adj', secondHalo
           if( secondHalo )                                             &
                call extrapolate2ndHaloAdj(nnbcs,icBeg, icEnd, jcBeg, jcEnd,  &
                iOffset, jOffset, wAdj0, wAdj1, &
                wAdj2, pAdj0, pAdj1, pAdj2)
           !print *,'replaceadj', secondHalo
           call replaceBCStatesAdj(nnbcs,  wAdj0,wAdj1, wAdj2, wAdj3,&
                pAdj0,pAdj1, pAdj2, pAdj3,rlvAdj1, rlvAdj2,revAdj1, revAdj2,&
                iCell, jCell,kCell,&
                wAdj,pAdj,rlvAdj,revAdj,secondHalo,nn,level,sps,sps2)
           
        endif testFarfield
    
     endif

 enddo bocos
!print *,'nnendbcfar',nn,nnbcs
!close (UNIT=unitxAD)
       end subroutine bcFarfieldAdj
