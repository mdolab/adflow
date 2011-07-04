!
!      ******************************************************************
!      *                                                                *
!      * File:          farField Drag.f90                               *
!      * Author:        Gaetan Kenway,C.A.(Sandy) Mader                 *
!      * Starting date: 06-30-2011                                      *
!      * Last modified: 06-30-2011                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine farFieldDrag()
  !
  !      ******************************************************************
  !      *                                                                *
  !      * farFieldDrag compuetes the total drag on the body using        *
  !      * a far-field method                                             *
  !      *                                                                *
  !      ******************************************************************
  !
  use blockPointers
  use BCTypes
  use flowVarRefState
  use inputPhysics
  use inputTimeSpectral
  use communication
  implicit none

  ! Woring Variables

  integer(kind=intType) :: nn,level,sps
  integer(kind=intType) :: i,j,k,l,m,n,ierr,liftindex

  ! Temporary Real Arrays:
  real(kind=realType), dimension(:,:,:,:), allocatable :: V_wind,fvw
  real(kind=realType), dimension(:,:,:), allocatable :: ds,dH,du_ir,res
  real(kind=realType) :: gm1,fact
  real(kind=realType) :: alpha,beta

  ! Expansion Coefficients:
  real(kind=realType) :: ffp1,ffs1,ffh1,ffp2,ffs2,ffh2,ffps2,ffph2,ffsh2
  real(kind=realType) :: uouinf
  ! Ratios used in expansion:
  real(kind=realType) :: dPoP,dSoR,dHou2

  ! Drag Values:
  real(kind=realType) :: drag_local,drag

  ! Temp Variables for divergence calc
  real(kind=realType) :: qsp, qsm, rqsp, rqsm, porVel, porFlux
  real(kind=realType) :: pa, fs, sFace, vnp, vnm
  real(kind=realType) :: wx, wy, wz, rvol
  real(kind=realType) :: vp1,vm1,v(3),left,right,face,v1(3),v2
  ! Define coefficients:
  !mach needs to be mach+machgrid I think...
  
  gm1 = gammaConstant - 1.0_realType

  ffp1 = -1/(gammaConstant*Mach**2)
  ffs1 = -1/(gammaConstant*Mach**2)
  ffh1 = 1.0

  ffp2 = -(1 + gammaConstant * Mach**2)/(2*gammaConstant**2*Mach**4)
  ffs2 = -(1 + gm1*Mach**2)/(2*gammaConstant**2*Mach**4)
  ffh2 = -1.0/2.0

  ffps2 = -(1 + gm1*Mach**2)/(gammaConstant**2*Mach**4)
  ffph2 = 1/(gammaConstant * Mach**2)
  ffsh2 = 1/(gammaConstant * Mach**2)
!    print *,'coefficients',ffp1,ffs1,ffh1
!    print *,'coefficients2',ffp2,ffs2,ffh2
!    print *,'coefficients3',ffps2,ffph2,ffsh2
!    print *, 'uinf:',uinf
!    print *,'Pinf:',Pinf
!    print *,'rhoinf:',rhoinf
!    print *,'Gam:',GammaConstant
!    print *,'RGas:',RGas
   call getDirAngle(velDirFreeStream,liftDirection,liftIndex,alpha,beta)

!   print *,'alpha,beta:',alpha,beta

  level = 1
  
  drag_local = 0.0
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointers(nn,level,sps)

        ! Allocate memory for this block:
        allocate(V_wind(0:ib,0:jb,0:kb,3))
        allocate(ds(0:ib,0:jb,0:kb), &
                 dH(0:ib,0:jb,0:kb), &
                 du_ir(0:ib,0:jb,0:kb),&
                 fvw(0:ib,0:jb,0:kb,3),&
                 res(0:ib,0:jb,0:kb))
        V_wind(:,:,:,:) = 0.0
        ds(:,:,:) = 0.0
        dH(:,:,:) = 0.0
        du_ir(:,:,:) = 0.0
        fvw(:,:,:,:) = 0.0
        res(:,:,:) = 0.0
        ! Loop over owned cells + 1 level of halos:

        if (myid == 0 .and. nn==1) then
           print *,'nx,ny,nz:',nx,ny,nz
           print *,'il,jl,kl:',il,jl,kl
           print *,'ie,je,ke:',ie,je,ke
           print *,'ib,jb,kb:',ib,jb,kb
        end if


        do k=1,ke
           do j=1,je
              do i=1,ie
                 !make w relative velocity????w-sface?
                 ! Compute the three components of the wind-oriented velocities:
                 call getWindAxis(w(i,j,k,ivx:ivz),V_wind(i,j,k,:),alpha,liftIndex)

                 !The variation of Entropy wrt the free stream:
                 ds(i,j,k) = (RGas/gm1)*log((P(i,j,k)/Pinf)*(rhoInf/w(i,j,k,iRho))**gammaConstant)

                 ! The variation of Stagnation Enthalpy relative to free stream:
                 dH(i,j,k) = (gammaConstant/gm1)*(P(i,j,k)/w(i,j,k,iRho)-Pinf/rhoinf) + &
                      0.5*(V_wind(i,j,k,1)**2 + V_wind(i,j,k,2)**2 + V_wind(i,j,k,3)**2 - Uinf**2)

                 dPoP = (P(i,j,k)-Pinf)/Pinf
                 dSoR = ds(i,j,k)/RGas
                 dHou2 = dH(i,j,k)/Uinf**2

                 du_ir(i,j,k)= uinf*sqrt(1 + 2*dH(i,j,k)/uinf**2 - 2/(gm1*Mach**2)*((dPoP + 1)**(gm1/gammaConstant)*exp(dSoR*(gm1/gammaconstant))-1))-uinf
                 !du_ir(i,j,k)= uinf*sqrt(1 + 2*dH(i,j,k)/uinf**2 - 2/(gm1*Mach**2)*(exp(dSoR)**(gm1/gammaConstant)-1))-uinf
                 
!                  du_ir(i,j,k) = uInf*(ffp1*dPoP + ffs1*dSoR + ffH1*dHou2 + &
!                                       ffp2*dPoP + ffs2*dSoR + ffH2*dHou2 + &
!                                       ffps2*dPoP*dSoR + ffph2*dPoP*dHou2 + ffsh2*dSoR*dHou2)
                 


                 !if (myid == 0) then
                    !print *,ds(i,j,k)
                    !print *,dH(i,j,k)
                  !  print *,du_ir(i,j,k)
                 !end if
                 fvw(i,j,k,:) = -w(i,j,k,iRho) * du_ir(i,j,k) * V_wind(i,j,k,:)
                 !fvw(i,j,k,:) = w(i,j,k,ivx:ivz)
                 !fvw(i,j,k,:) = 1.0
!                  if (isnan(sum(fw(i,j,k,:)))) then
!                     fw(i,j,k,:) = 0.0
!                  end if
              end do
           end do
        end do

!         do k=2,kl
!            do j=1,je
!               do i=2,il
!                  !make w relative velocity????w-sface?
!                  ! Compute the three components of the wind-oriented velocities:
!                  call getWindAxis(w(i,j,k,ivx:ivz),V_wind(i,j,k,:),alpha,liftIndex)

!                  !The variation of Entropy wrt the free stream:
!                  ds(i,j,k) = (RGas/gm1)*log((P(i,j,k)/Pinf)*(rhoInf/w(i,j,k,iRho))**gammaConstant)

!                  ! The variation of Stagnation Enthalpy relative to free stream:
!                  dH(i,j,k) = (gammaConstant/gm1)*(P(i,j,k)/w(i,j,k,iRho)-Pinf/rhoinf) + &
!                       0.5*(V_wind(i,j,k,1)**2 + V_wind(i,j,k,2)**2 + V_wind(i,j,k,3)**2 - Uinf**2)

!                  dPoP = (P(i,j,k)-Pinf)/Pinf
!                  dSoR = ds(i,j,k)/RGas
!                  dHou2 = dH(i,j,k)/Uinf**2

!                  !du_ir(i,j,k)= uinf*sqrt(1 + 2*dH(i,j,k)/uinf**2 - 2/(gm1*Mach**2)*((dPoP + 1)**(gm1/gammaConstant)*exp(dSoR)**(gm1/gammaconstant)-1))-uinf
!                  du_ir(i,j,k)= uinf*sqrt(1 + 2*dH(i,j,k)/uinf**2 - 2/(gm1*Mach**2)*(exp(dSoR)**(gm1/gammaConstant)-1))-uinf

!                  if (isnan(du_ir(i,j,k))) then
!                     print *,myid,i,j,k
!                     print *,il,jl,kl
!                     du_ir(i,j,k) = 0.0
!                  end if



!                  fvw(i,j,k,:) = -w(i,j,k,iRho) * du_ir(i,j,k) * V_wind(i,j,k,:)


!               end do
!            end do
!         end do

!         do k=1,ke
!            do j=2,jl
!               do i=2,il
!                  !make w relative velocity????w-sface?
!                  ! Compute the three components of the wind-oriented velocities:
!                  call getWindAxis(w(i,j,k,ivx:ivz),V_wind(i,j,k,:),alpha,liftIndex)

!                  !The variation of Entropy wrt the free stream:
!                  ds(i,j,k) = (RGas/gm1)*log((P(i,j,k)/Pinf)*(rhoInf/w(i,j,k,iRho))**gammaConstant)

!                  ! The variation of Stagnation Enthalpy relative to free stream:
!                  dH(i,j,k) = (gammaConstant/gm1)*(P(i,j,k)/w(i,j,k,iRho)-Pinf/rhoinf) + &
!                       0.5*(V_wind(i,j,k,1)**2 + V_wind(i,j,k,2)**2 + V_wind(i,j,k,3)**2 - Uinf**2)

!                  dPoP = (P(i,j,k)-Pinf)/Pinf
!                  dSoR = ds(i,j,k)/RGas
!                  dHou2 = dH(i,j,k)/Uinf**2

!                  !du_ir(i,j,k)= uinf*sqrt(1 + 2*dH(i,j,k)/uinf**2 - 2/(gm1*Mach**2)*((dPoP + 1)**(gm1/gammaConstant)*exp(dSoR)**(gm1/gammaconstant)-1))-uinf
!                  du_ir(i,j,k)= uinf*sqrt(1 + 2*dH(i,j,k)/uinf**2 - 2/(gm1*Mach**2)*(exp(dSoR)**(gm1/gammaConstant)-1))-uinf
 
!                  if (isnan(du_ir(i,j,k))) then
!                     print *,myid,i,j,k
!                     print *,il,jl,kl
!                     du_ir(i,j,k) = 0.0
!                  end if
                 
!                 fvw(i,j,k,:) = -w(i,j,k,iRho) * du_ir(i,j,k) * V_wind(i,j,k,:)
                
!               end do
!            end do
!         end do

!        print *,'myid sum:',myid,sum(fvw(:,:,:,:))




!                  du_ir(i,j,k) = uInf*(ffp1*dPoP + ffs1*dSoR + ffH1*dHou2 + &
!                                       ffp2*dPoP + ffs2*dSoR + ffH2*dHou2 + &
!                                       ffps2*dPoP*dSoR + ffph2*dPoP*dHou2 + ffsh2*dSoR*dHou2)


!                  if (myid == 0 .and. nn==1 .and. i==2 .and. j==2 .and. k==2) then
!                     print *,'Check ds calc:'
!                     print *,'RGas=',RGas
!                     print *,'gm1=',gm1
!                     print *,'P=',P(i,j,k)
!                     print *,'Pinf=',Pinf
!                     print *,'rhoInf=',rhoInf
!                     print *,'rho=',w(i,j,k,iRho)
!                     print *,'ds=',ds(i,j,k)
!                     print *,'Vwind[0]=',V_wind(i,j,k,1)
!                     print *,'Vwind[1]=',V_wind(i,j,k,2)
!                     print *,'Vwind[2]=',V_wind(i,j,k,3)
!                     print *,'dH=',dH(i,j,k)
!                     print *,'du1=',uouinf
!                     print *,'du2=',du_ir(i,j,k)
!                  end if

!                  if (myid == 0 .and. nn==1 .and. i==2 .and. j==2 .and. k==2) then
!                     print *,'P=',P(i,j,k)
!                     print *,'RGas=',RGas
!                     print *,'Pinf=',Pinf
!                     print *,'rhoInf=',rhoInf
!                     print *,'rho=',w(i,j,k,iRho)
!                     print *,'ds=',ds(i,j,k)
!                     print *,'V_wind(1)=',V_wind(i,j,k,1)
!                     print *,'V_wind(2)=',V_wind(i,j,k,2)
!                     print *,'V_wind(3)=',V_wind(i,j,k,3)
!                     print *,'Uinf=',uinf
!                     print *,'dH=',dH(i,j,k)
!                     print *,'one=',uouinf
!                     print *,'two=',du_ir(i,j,k)
!                  end if

                 !print *,uouinf,du_ir(i,j,k)
                 
                 ! We should now be in a position to integrate the
                 ! irreversible drag over the volume:
                 !         /
                 ! D_irr = | div ( -rho * du_ir * V ) dV
                 !         /
                 !          V

                 ! produce the fvw vector:
                 !fvw(i,j,k,:) = -w(i,j,k,iRho) * du_ir(i,j,k) * V_wind(i,j,k,:)
                 !fvw(i,j,k,:) = 1.0
                 !print *,myid,nn,du_ir(i,j,k)
!                  if (myid == 0) then
!                     print *,i,j,k
!                     print *,fvw(i,j,k,:)
!                  end if
!               end do
!            end do
!         end do


        do k=2,kl
           n = k -1
           do j=2,jl
              m = j -1
              do i=2,il
                 l = i -1

                 ! I Left
                 v = (fvw(i,j,k,:)+fvw(i-1,j,k,:))*0.5
                 face  = v(1)*sI(l,j,k,1) + v(2)*sI(l,j,k,2) + v(3)*sI(l,j,k,3)
                 res(i,j,k) = res(i,j,k) - face*porI(l,j,k)

                 ! I Right
                 v = (fvw(i+1,j,k,:)+fvw(i,j,k,:))*0.5
                 face  = v(1)*sI(i,j,k,1) + v(2)*sI(i,j,k,2) + v(3)*sI(i,j,k,3)
                 res(i,j,k) = res(i,j,k) + face*porI(i,j,k)

                 ! J Left
                 v = (fvw(i,j,k,:)+fvw(i,j-1,k,:))*0.5
                 face  = v(1)*sJ(i,m,k,1) + v(2)*sJ(i,m,k,2) + v(3)*sJ(i,m,k,3)
                 res(i,j,k) = res(i,j,k) - face*porJ(i,m,k)

                 ! J Right
                 v = (fvw(i,j+1,k,:)+fvw(i,j,k,:))*0.5
                 face  = v(1)*sJ(i,j,k,1) + v(2)*sJ(i,j,k,2) + v(3)*sJ(i,j,k,3)
                 res(i,j,k) = res(i,j,k) + face*porJ(i,j,k)

                 ! K Left
                 v = (fvw(i,j,k,:)+fvw(i,j,k-1,:))*0.5
                 face  = v(1)*sK(i,j,n,1) + v(2)*sK(i,j,n,2) + v(3)*sK(i,j,n,3)
                 res(i,j,k) = res(i,j,k) - face*porK(i,j,n)

                 ! K Right
                 v = (fvw(i,j,k+1,:)+fvw(i,j,k,:))*0.5
                 face  = v(1)*sK(i,j,k,1) + v(2)*sK(i,j,k,2) + v(3)*sK(i,j,k,3)
                 res(i,j,k) = res(i,j,k) + face*porK(i,j,k)

              end do
           end do
        end do
        ! Now we have the divergence of fvw computed in res, we simply sum up the contributions:

        !print *,'b,d:',nbkglobal,sum(fvw)
        do k=2,kl
           do j=2,jl
              do i=2,il
                 drag_local = drag_local + res(i,j,k)
              end do
           end do
        end do
      
        deallocate(V_wind,ds,dH,du_ir,fvw,res)

     end do
  end do

  ! Reduce the drag to root proc:
  call mpi_reduce(drag_local,drag,1,sumb_real,mpi_sum,0,SUmb_comm_world,ierr)
  fact = two/(gammaInf*pInf*MachCoef*MachCoef &
            *surfaceRef*LRef*LRef)
  fact = two/(surfaceRef)
  if (myid == 0) then
     print *,'fact',fact
  end if
  !print *,'factff',fact,two,gammaInf,pInf,pRef,MachCoef,MachCoef,surfaceRef,LRef,LRef
  if (myid == 0) then
     print *,'Irreversable drag:',drag
     print *,'Irreversable CD',drag*fact
  end if

end subroutine farFieldDrag

subroutine getWindAxis(V1,V2,alpha,liftIndex)

  ! Return vector V1 specified in body axis frame in wind axis frame. Only works for alpha (and not beta)
  use constants
  
  real(kind=realType) :: V1(3),V2(3),alpha
  integer(kind=intType)::liftIndex
  !print *,'getting wind axis',v1,v2,alpha,liftindex
  if (liftIndex == 2) then
     
     call VectorRotation(V2(1),V2(2),V2(3),3, alpha, V1(1),V1(2),v1(3))

  else if (liftIndex==3) then

     call VectorRotation(V2(1),V2(2),V2(3),2,-alpha, V1(1),V1(2),v1(3))
  end if

end subroutine getWindAxis


