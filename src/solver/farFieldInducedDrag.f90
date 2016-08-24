! #############################
! THIS ROUTINE DOES NOT WORK
! #############################
!
subroutine farFieldInducedDrag(value)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * farFieldInducedDrag computes the induced drag at the farfield  *
  !      * a far-field method                                             *
  !      *                                                                *
  !      ******************************************************************
  !
  use blockPointers
  use flowVarRefState
  use inputPhysics
  use inputTimeSpectral
  use communication
  implicit none


  real(kind=realType) :: value
  ! Boring Variables

!   integer(kind=intType) :: nn,level,sps,mm,m,l
!   integer(kind=intType) :: i,j,k,ierr,liftindex
!   integer(kind=intType) :: groundLevel
!   ! Expansion Coefficients:
!   real(kind=realType) :: ffp1,ffs1,ffh1,ffp2,ffs2,ffh2,ffps2,ffph2,ffsh2

!   ! Ratios used in expansion:
!   real(kind=realType) :: dPoP,dSoR,dHou2,popInf,uouinf,uouinf2
  
!   ! Drag Values:
!   !real(kind=realType),dimension(3) :: drag_local,drag
!   real(kind=realType) :: drag_local,drag
!   real(kind=realType), dimension(3) :: V_wind,fi,f
!   real(kind=realType), dimension(3) :: norm_wind,ss_wind
!   real(kind=realType), dimension(3) :: vAvg,vrot
!   real(kind=realType) :: rhoAvg,pAvg
!   real(kind=realType) :: ds,dH,du_ir,dustar
!   real(kind=realType) :: gm1,fact,ovfact,a,vref,area
!   real(kind=realType) :: alpha,beta
!   real(kind=realType), dimension(3,3) :: vTen,pTen
!   real(kind=realType), dimension(3) :: force,force_local
!   real(kind=realType)::CD,CL
!   !pointers
!   real(kind=realType), dimension(:,:),   pointer :: pp2, pp1
!   real(kind=realType), dimension(:,:),   pointer :: rho2, rho1
!   real(kind=realType), dimension(:,:,:),   pointer :: vv2, vv1
!   real(kind=realType), dimension(:,:,:), pointer :: ss
!   real(kind=realType), dimension(:,:,:), pointer :: norm
! !Begin Execution
!   ! Define coefficients:
!   gm1 = gammaConstant - 1.0_realType
 
!   !mach needs to be mach+machgrid I think...
  
!   ffp1 = -1/(gammaConstant*Mach**2)
!   ffs1 = -1/(gammaConstant*Mach**2)
!   ffh1 = 1.0

!   ffp2 = -(1 + gammaConstant * Mach**2)/(2*gammaConstant**2*Mach**4)
!   ffs2 = -(1 + gm1*Mach**2)/(2*gammaConstant**2*Mach**4)
!   ffh2 = -1.0/2.0

!   ffps2 = -(1 + gm1*Mach**2)/(gammaConstant**2*Mach**4)
!   ffph2 = 1/(gammaConstant * Mach**2)
!   ffsh2 = 1/(gammaConstant * Mach**2)

!   call getDirAngle(velDirFreeStream,liftDirection,liftIndex,alpha,beta)
  
!   drag = 0.0
!   drag_local = 0.0
!   force(:) = 0.0

!   groundlevel = 1
!   spectralLoop: do sps=1,nTimeIntervalsSpectral
     
!      ! Loop over the blocks.
     
!      domains: do nn=1,nDom
        
!         ! Set the pointers for this block.
        
!         call setPointers(nn, groundLevel, sps)
        
        
!         ! Loop over the boundary subfaces of this block.
        
!         bocos: do mm=1,nBocos
! !
! !        ****************************************************************
! !        *                                                              *
! !        * Integrate the induced contribution over the farfield.        *
! !        *                                                              *
! !        ****************************************************************
! !
        
!          farfieldForce: if(BCType(mm) == Farfield) then!symm)then!eulerwall)then!Farfield) then!

!            ! Subface is a farfield boundary.

!            ! Set a bunch of pointers depending on the face id to make
!            ! a generic treatment possible. The routine setBcPointers
!            ! is not used, because quite a few other ones are needed.

!            select case (BCFaceID(mm))

!              case (iMin)
!                  !print *,'imin'
!                pp2  => p(2,1:,1:);      pp1  => p(1,1:,1:)
!                rho2 => w(2,1:,1:,irho); rho1 => w(1,1:,1:,irho)
!                vv2 => w(2,1:,1:,ivx:ivz); vv1 => w(1,1:,1:,ivx:ivz)
!                ss   => si(1,:,:,:)
!                fact = -one

             
!              !===========================================================

!              case (iMax)
!                 !print *,'imax'
!                pp2  => p(il,1:,1:);      pp1  => p(ie,1:,1:)
!                rho2 => w(il,1:,1:,irho); rho1 => w(ie,1:,1:,irho)
!                vv2 => w(il,1:,1:,ivx:ivz); vv1 => w(ie,1:,1:,ivx:ivz)
!                ss   => si(il,:,:,:)
!                fact = one

              
!              !===========================================================

!              case (jMin)
!                   !print *,'jmin'
!                pp2  => p(1:,2,1:);      pp1  => p(1:,1,1:)
!                rho2 => w(1:,2,1:,irho); rho1 => w(1:,1,1:,irho)
!                vv2  => w(1:,2,1:,ivx:ivz); vv1 => w(1:,1,1:,ivx:ivz)
!                ss   => sj(:,1,:,:)
!                fact = -one


!              !===========================================================

!              case (jMax)
!                 !print *,'jmax'
!                pp2  => p(1:,jl,1:);      pp1  => p(1:,je,1:)
!                rho2 => w(1:,jl,1:,irho); rho1 => w(1:,je,1:,irho)
!                vv2  => w(1:,jl,1:,ivx:ivz); vv1 => w(1:,je,1:,ivx:ivz)
!                ss   => sj(:,jl,:,:)
!                fact = one

               
!              !===========================================================

!              case (kMin)
!                 !   print *,'kmin'
!                pp2  => p(1:,1:,2);      pp1  => p(1:,1:,1)
!                rho2 => w(1:,1:,2,irho); rho1 => w(1:,1:,1,irho)
!                vv2  => w(1:,1:,2,ivx:ivz); vv1 => w(1:,1:,1,ivx:ivz)
!                ss   => sk(:,:,1,:)
!                fact = -one

              
!              !===========================================================

!              case (kMax)
!                 ! print *,'kmax'
!                pp2  => p(1:,1:,kl);      pp1  => p(1:,1:,ke)
!                rho2 => w(1:,1:,kl,irho); rho1 => w(1:,1:,ke,irho)
!                vv2  => w(1:,1:,kl,ivx:ivz); vv1 => w(1:,1:,ke,ivx:ivz)
!                ss   => sk(:,:,kl,:)
!                fact = one

              
!            end select


!            !set the pointer to the normals for this face
!            norm => BCData(mm)%norm
           
!            ! Loop over the quadrilateral faces of the subface. Note
!            ! that the nodal range of BCData must be used and not the
!            ! cell range, because the latter may include the halo's in i
!            ! and j-direction. The offset +1 is there, because inBeg and
!            ! jnBeg refer to nodal ranges and not to cell ranges.
           
!            do j=(BCData(mm)%jnBeg+1),BCData(mm)%jnEnd
!              do i=(BCData(mm)%inBeg+1),BCData(mm)%inEnd

!                 !Calculate the velocities at the farfield surface by averageing the
!                 ! cell inside and outside.
!                 vAvg(:) = half*(vv1(i,j,:)+vv2(i,j,:))
!                 !print *,'vavg',vavg,sqrt(vavg(1)**2+vavg(2)**2+vavg(3)**2),sqrt(gammaInf*pRef/rhoRef),vv1(i,j,:),vv2(i,j,:)
!                 pAvg = half*(pp1(i,j)+pp2(i,j))
!                 rhoAvg = half*(rho1(i,j)+rho2(i,j))
                
!                 !perform the basic rho VV+P integration over the far field

!                 vRef =sqrt(gammaInf*pRef/rhoRef)/sqrt(gammaInf*pInf/rhoInf)
!                 !print *,'vavg',vavg,vavg*vref
!                 !create the outer produce tensor and the pressure tensor
!                 do l = 1,3
!                    do m = 1,3
!                       vTen(l,m)  = fact*rhoAvg*rhoRef*vAvg(l) * vAvg(m)*vRef**2
!                       !vTen(l,m)  = fact*rhoAvg*vAvg(l) * vAvg(m)
!                    end do
!                 end do
!                 !print *,'vten',vten
!                 !stop
!                 !print *,'refs',pref, rhoref*vref**2
!                 do l = 1,3
!                    pTen(l,l) = 0.0!fact*(pAvg-pInf)*pRef
!                    !pTen(l,l) = fact*(pAvg-pInf)
!                 end do
!                 !print *,'pten',pten
!                 !add in tauTen for viscous

!                 !Sum the tensors
!                 vTen = vTen+pTen!+tauTen
!                 !print *,'vten',vten
!                 force_local = force_local- matmul(vTen,ss(i,j,:))
!                 !print *,'forcelocal',force_local
!                 !De normalize the quanities
!                 vAvg = vAvg*vRef
!                 pAvg = pAvg*pRef
!                 rhoAvg = rhoAvg*rhoRef

!                 !make w relative velocity????w-sface?
!                 ! Compute the three components of the wind-oriented velocities:
!                 call getWindAxis(vAvg(:),V_wind(:),alpha,liftIndex)
!                 call getWindAxis(vAvg(:)-(/1,0,0/)*(Uinf*vRef),Vrot(:),alpha,liftIndex)
!                 !print *,'v',vAvg,'vw',v_wind,alpha
             
!                 !The variation of Entropy wrt the free stream:
!                 ds = (RGasDim/gm1)*log((Pavg/(Pinf*pRef))*((rhoInf*rhoRef)/rhoAvg)**gammaConstant)
!                 !print *,'RGas',rgas,gm1
!                 !print *,'ds',ds!,(RGas/gm1),(Pavg/Pinf),(rhoInf/rhoAvg),gammaConstant
!                 !print *,'refs',pref/rhoref, vref**2
!                 ! The variation of Stagnation Enthalpy relative to free stream:
!                 dH = (gammaConstant/gm1)*(PAvg/rhoAvg-(Pinf*pref)/(rhoinf*rhoRef)) + &
!                      0.5*(( V_wind(1)**2 + V_wind(2)**2 + V_wind(3)**2) - (Uinf*vRef)**2)
!                 !print *,'vels', 0.5*(( V_wind(1)**2 + V_wind(2)**2 + V_wind(3)**2) - Uinf**2),(PAvg/rhoAvg-Pinf/rhoinf), (gammaConstant/gm1)
                
!                 popInf  = exp(-ds/rgasdim)*(1+((gm1* Mach**2)/2.0)*(1-(V_wind(1)**2 +&
!                      V_wind(2)**2 + V_wind(3)**2)/(Uinf*vRef)**2+&
!                      2*(dh/(Uinf*vRef)**2)))**(gammaConstant/gm1)
!                 !print *,'t1',(gm1* Mach**2)/2.0
!                 !print *,'t2',(V_wind(1)**2 + V_wind(2)**2 + V_wind(3)**2)/(Uinf*vRef)**2
!                 !print *,'t3',2*(dh/(Uinf*vRef)**2)
!                 !print *,'test', exp(-ds/rgasdim),((1+(gm1* Mach**2)/2.0)*(1-(V_wind(1)**2 + V_wind(2)**2 + V_wind(3)**2)/(Uinf*vRef)**2)+2*(dh/(Uinf*vRef)**2)),(gammaConstant/gm1)
!                 !print *,'test2', (V_wind(2)**2 + V_wind(3)**2)/(Uinf*vRef)**2
!                 !print *,'popInf',popInf,pAvg/(pInf*pref),(pAvg/(pInf*pref))/popInf

!                 uouInf = sqrt(1+2*(dh/(Uinf*vRef)**2)-(2.0/(gm1* Mach**2))*&
!                      (popInf**(gm1/gammaconstant)*exp(ds/rgasDim)**(gm1/gammaconstant)-1))
                
!                 !print *,'dh',dh
                
!                 dPoP = (Pavg-Pinf*pRef)/(Pinf*pRef)!(Pavg-Pinf)/Pinf
!                 !print *,'dpop',dpop
!                 dSoR = ds/RGasDim
!                 !print *,'dSoR',dsor
!                 dHou2 = dH/(Uinf*vref)**2
!                 !print *,'dhou2',dhou2 ,dH,Uinf**2
!                 uouinf2 = 1+ffp1*dPoP + ffs1*dSoR + ffH1*dHou2 + &
!                      ffp2*dPoP + ffs2*dSoR + ffH2*dHou2 + &
!                      ffps2*dPoP*dSoR + ffph2*dPoP*dHou2 + ffsh2*dSoR*dHou2
!                 !print *,'uouInf',uouInf,V_wind(1)/(Uinf*vRef),uouInf2
!                 !print *,'error',dPoP**3,dSoR**3,dhou2**3,abs(dPoP)**3+abs(dSoR)**3+abs(dhou2)**3
!                 du_ir = uInf*(ffp1*dPoP + ffs1*dSoR + ffH1*dHou2 + &
!                      ffp2*dPoP + ffs2*dSoR + ffH2*dHou2 + &
!                      ffps2*dPoP*dSoR + ffph2*dPoP*dHou2 + ffsh2*dSoR*dHou2)
!                ! print *,'du_ir', du_ir, uInf,ffp1*dPoP ,ffs1*dSoR, ffH1*dHou2,&
!                !      ffp2*dPoP,ffs2*dSoR, ffH2*dHou2, &
!                !      ffps2*dPoP*dSoR,ffph2*dPoP*dHou2, ffsh2*dSoR*dHou2
!                 !stop

!                 !method from Onera Paper
!                 dustar = v_wind(1)-uinf*vRef-(uinf*vRef*(dSoR/(gammaConstant*Mach**2)-dHou2))

!                 fi(1) = (rhoInf*rhoRef/2.0)*(v_wind(2)**2+v_wind(3)**2-(1-Mach**2)*dustar**2)
!                 fi(2) = (rhoInf*rhoRef/2.0)*(-2*dustar*v_wind(2))
!                 fi(3) = (rhoInf*rhoRef/2.0)*(-2*dustar*v_wind(3))
                
!                 fi(:) = 0.0
!                 fi(1) = 1.0
!                 fi(2) = 1.0
!                 fi = v_wind
!                 !fi(3) = 1.0
!                 !rotate vector quantities to the wind fram
                
!                 call getWindAxis(norm(i,j,:),norm_wind(:),alpha,liftIndex)
!                 call getWindAxis(ss(i,j,:),ss_wind(:),alpha,liftIndex)
!                 !fi = -rhoAvg*(v_wind(1)-uInf-du_ir)*v_wind(:)-(pAvg-pInf)*norm(i,j,:)
!                 !f = -rhoAvg*(v_wind(1))*v_wind(:)!fact*-rhoAvg*(v_wind(1)-(uInf*vRef))*v_wind(:)!-fact*(pAvg-(pInf*pRef))*(/1,0,0/)!-(uInf*vRef)
!                 f = rhoAvg*vAvg(:)!v_wind(:)
!                 !f = fact*-rhoAvg*(-(uInf*vRef))*v_wind(:)!v_wind(1)
!                 !print *,'deltau',v_wind(1)-(uInf*vRef),fact,v_wind(1),(uInf*vRef)
!                 !print *,'norm',norm(i,j,:)
!                 !print *,'f',f,'comp',(pAvg-(pInf*pRef))*(/1,0,0/),norm_wind(:),norm_wind(1)**2+norm_wind(2)**2+norm_wind(3)**2!-rhoAvg,(v_wind(1)-(uInf*vRef)),v_wind(:),(pAvg-(pInf*pRef)),norm(i,j,:)
!                 !print *,'fi',fi,-rhoAvg,v_wind(1),uInf,du_ir,'v',v_wind(:),pAvg,pInf,norm(i,j,:)
!                 !print *,'fsum',f(1)*(norm(i,j,1))**2+f(2)*(norm(i,j,2))**2+f(3)*(norm(i,j,3))**2
!                 !area = sqrt(ss(i,j,1)**2+ss(i,j,2)**2+ss(i,j,3)**2)
!                 !print *,'normcheck',fact*ss(i,j,:),'area',norm(i,j,:)*area
!                 !print *,'normcheck',fact*ss_wind(:),'area',norm_wind(:)*area
!                 !do m = 1,3
!                 do k=1,3
!                    !!drag_local(m) = drag_local(m)+fi(k)*ss(i,j,k)!*norm(i,j,k)
!                    !!drag_local = drag_local+f(k)*ss(i,j,k)!*norm(i,j,k)
!                    drag_local = drag_local+f(k)*fact*ss(i,j,k)!fact*!ss_wind(k)
!                    !drag_local = drag_local+fi(k)*fact*ss_wind(k)
!                    !drag_local = drag_local+fi(k)*norm_wind(k)*area!norm(i,j,k)
!                    !print *,'fi',fi(k),norm(i,j,k),ss_wind(k),ss(i,j,k)
                   
!                    !print *,'drag',f(k)*ss(i,j,k),f(k),ss(i,j,k)
!                 end do
!                 !end do
!                 !print *,'areas',sqrt(ss(i,j,1)**2+ss(i,j,2)**2+ss(i,j,3)**2),sqrt(ss_wind(1)**2+ss_wind(2)**2+ss_wind(3)**2)
              
!             enddo
!          enddo


!       endif farfieldForce

!    enddo bocos
! end do domains
! end do spectralLoop

!   ! Reduce the drag to root proc:
!   call mpi_reduce(drag_local,drag,1,sumb_real,mpi_sum,0,SUmb_comm_world,ierr)
  
!   fact = two/(gammaInf*pInf*MachCoef*MachCoef &
!             *surfaceRef*LRef*LRef)
!   a = sqrt(gammaInf*pRef/rhoRef)
!   ovfact = rhoRef*vRef!0.5*rhoref*(machCoef*a)**2*surfaceRef
!   if (myid == 0) then
!      print *,'fact',fact,ovfact
!   end if
!   !reduce force computation
!   call mpi_reduce(force_local,force,3,sumb_real,mpi_sum,0,SUmb_comm_world,ierr)
  
!   !print *,'fact',fact
!   !print *,'factff',fact,two,gammaInf,pInf,pRef,MachCoef,MachCoef,surfaceRef,LRef,LRef
!   CL = force(1)*liftDirection(1) &
!        + force(2)*liftDirection(2) &
!        + force(3)*liftDirection(3)

!   CD = force(1)*dragDirection(1) &
!        + force(2)*dragDirection(2) &
!        + force(3)*dragDirection(3)
  
!   if (myid == 0) then
!      print *,'Induced drag:',drag
!      print *,'Induced CD',drag*fact,drag/ovfact
!      print *,'force',force,force*fact,force/ovfact
!      print *,'force CL',CL,CL*fact,CL/ovfact
!      print *,'force CD',CD,CD*fact,CD/ovfact
!   end if

!   value = CD/ovfact
end subroutine farFieldInducedDrag
