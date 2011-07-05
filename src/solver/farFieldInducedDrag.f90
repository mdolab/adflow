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
subroutine farFieldInducedDrag()
  !
  !      ******************************************************************
  !      *                                                                *
  !      * farFieldInducedDrag computes the induced drag at the farfield  *
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


  ! Boring Variables

  integer(kind=intType) :: nn,level,sps,mm,m,l
  integer(kind=intType) :: i,j,k,ierr,liftindex
  integer(kind=intType) :: groundLevel
  ! Expansion Coefficients:
  real(kind=realType) :: ffp1,ffs1,ffh1,ffp2,ffs2,ffh2,ffps2,ffph2,ffsh2

  ! Ratios used in expansion:
  real(kind=realType) :: dPoP,dSoR,dHou2
  
  ! Drag Values:
  real(kind=realType),dimension(3) :: drag_local,drag
  real(kind=realType), dimension(3) :: V_wind,fi
  real(kind=realType), dimension(3) :: vAvg
  real(kind=realType) :: rhoAvg,pAvg
  real(kind=realType) :: ds,dH,du_ir
  real(kind=realType) :: gm1,fact,ovfact,a,vref
  real(kind=realType) :: alpha,beta
  real(kind=realType), dimension(3,3) :: vTen,pTen
  real(kind=realType), dimension(3) :: force,force_local
  real(kind=realType)::CD,CL
  !pointers
  real(kind=realType), dimension(:,:),   pointer :: pp2, pp1
  real(kind=realType), dimension(:,:),   pointer :: rho2, rho1
  real(kind=realType), dimension(:,:,:),   pointer :: vv2, vv1
  real(kind=realType), dimension(:,:,:), pointer :: ss
  real(kind=realType), dimension(:,:,:), pointer :: norm
!Begin Execution
  ! Define coefficients:
  gm1 = gammaConstant - 1.0_realType
 
  !mach needs to be mach+machgrid I think...
  
  ffp1 = -1/(gammaConstant*Mach**2)
  ffs1 = -1/(gammaConstant*Mach**2)
  ffh1 = 1.0

  ffp2 = -(1 + gammaConstant * Mach**2)/(2*gammaConstant**2*Mach**4)
  ffs2 = -(1 + gm1*Mach**2)/(2*gammaConstant**2*Mach**4)
  ffh2 = -1.0/2.0

  ffps2 = -(1 + gm1*Mach**2)/(gammaConstant**2*Mach**4)
  ffph2 = 1/(gammaConstant * Mach**2)
  ffsh2 = 1/(gammaConstant * Mach**2)

  call getDirAngle(velDirFreeStream,liftDirection,liftIndex,alpha,beta)
  
  drag(:) = 0.0
  drag_local(:) = 0.0
  force(:) = 0.0

  groundlevel = 1
  spectralLoop: do sps=1,nTimeIntervalsSpectral
     
     ! Loop over the blocks.
     
     domains: do nn=1,nDom
        
        ! Set the pointers for this block.
        
        call setPointers(nn, groundLevel, sps)
        
        
        ! Loop over the boundary subfaces of this block.
        
        bocos: do mm=1,nBocos
!
!        ****************************************************************
!        *                                                              *
!        * Integrate the induced contribution over the farfield.        *
!        *                                                              *
!        ****************************************************************
!
        
         farfieldForce: if(BCType(mm) == Farfield) then!eulerwall)then!Farfield) then!

           ! Subface is a farfield boundary.

           ! Set a bunch of pointers depending on the face id to make
           ! a generic treatment possible. The routine setBcPointers
           ! is not used, because quite a few other ones are needed.

           select case (BCFaceID(mm))

             case (iMin)
                 !print *,'imin'
               pp2  => p(2,1:,1:);      pp1  => p(1,1:,1:)
               rho2 => w(2,1:,1:,irho); rho1 => w(1,1:,1:,irho)
               vv2 => w(2,1:,1:,ivx:ivz); vv1 => w(1,1:,1:,ivx:ivz)
               ss   => si(1,:,:,:)
               fact = -one

             
             !===========================================================

             case (iMax)
                !print *,'imax'
               pp2  => p(il,1:,1:);      pp1  => p(ie,1:,1:)
               rho2 => w(il,1:,1:,irho); rho1 => w(ie,1:,1:,irho)
               vv2 => w(il,1:,1:,ivx:ivz); vv1 => w(ie,1:,1:,ivx:ivz)
               ss   => si(il,:,:,:)
               fact = one

              
             !===========================================================

             case (jMin)
                  !print *,'jmin'
               pp2  => p(1:,2,1:);      pp1  => p(1:,1,1:)
               rho2 => w(1:,2,1:,irho); rho1 => w(1:,1,1:,irho)
               vv2  => w(1:,2,1:,ivx:ivz); vv1 => w(1:,1,1:,ivx:ivz)
               ss   => sj(:,1,:,:)
               fact = -one


             !===========================================================

             case (jMax)
                !print *,'jmax'
               pp2  => p(1:,jl,1:);      pp1  => p(1:,je,1:)
               rho2 => w(1:,jl,1:,irho); rho1 => w(1:,je,1:,irho)
               vv2  => w(1:,jl,1:,ivx:ivz); vv1 => w(1:,je,1:,ivx:ivz)
               ss   => sj(:,jl,:,:)
               fact = one

               
             !===========================================================

             case (kMin)
                !   print *,'kmin'
               pp2  => p(1:,1:,2);      pp1  => p(1:,1:,1)
               rho2 => w(1:,1:,2,irho); rho1 => w(1:,1:,1,irho)
               vv2  => w(1:,1:,2,ivx:ivz); vv1 => w(1:,1:,1,ivx:ivz)
               ss   => sk(:,:,1,:)
               fact = -one

              
             !===========================================================

             case (kMax)
                ! print *,'kmax'
               pp2  => p(1:,1:,kl);      pp1  => p(1:,1:,ke)
               rho2 => w(1:,1:,kl,irho); rho1 => w(1:,1:,ke,irho)
               vv2  => w(1:,1:,kl,ivx:ivz); vv1 => w(1:,1:,ke,ivx:ivz)
               ss   => sk(:,:,kl,:)
               fact = one

              
           end select


           !set the pointer to the normals for this face
           norm => BCData(mm)%norm
           
           ! Loop over the quadrilateral faces of the subface. Note
           ! that the nodal range of BCData must be used and not the
           ! cell range, because the latter may include the halo's in i
           ! and j-direction. The offset +1 is there, because inBeg and
           ! jnBeg refer to nodal ranges and not to cell ranges.
           
           do j=(BCData(mm)%jnBeg+1),BCData(mm)%jnEnd
             do i=(BCData(mm)%inBeg+1),BCData(mm)%inEnd

                !Calculate the velocities at the farfield surface by averageing the
                ! cell inside and outside.
                vAvg(:) = half*(vv1(i,j,:)+vv2(i,j,:))
                !print *,'vavg',vavg,sqrt(vavg(1)**2+vavg(2)**2+vavg(3)**2),sqrt(gammaInf*pRef/rhoRef),vv1(i,j,:),vv2(i,j,:)
                pAvg = half*(pp1(i,j)+pp2(i,j))
                rhoAvg = half*(rho1(i,j)+rho2(i,j))
                
                !perform the basic rho VV+P integration over the far field

                vRef =sqrt(gammaInf*pRef/rhoRef)/sqrt(gammaInf*pInf/rhoInf)
                !print *,'vavg',vavg,vavg*vref
                !create the outer produce tensor and the pressure tensor
                do l = 1,3
                   do m = 1,3
                      !vTen(l,m)  = fact*rhoAvg*rhoRef*vAvg(l) * vAvg(m)*vRef**2
                      vTen(l,m)  = fact*rhoAvg*vAvg(l) * vAvg(m)
                   end do
                end do
                !print *,'vten',vten
                !stop
                !print *,'refs',pref, rhoref*vref**2
                do l = 1,3
                   !pTen(l,l) = fact*(pAvg-pInf)*pRef
                   pTen(l,l) = fact*(pAvg-pInf)
                end do
                
                !add in tauTen for viscous

                !Sum the tensors
                vTen = vTen+pTen!+tauTen
                
                force_local = force_local- matmul(vTen,ss(i,j,:))
                
                
                
                !make w relative velocity????w-sface?
                ! Compute the three components of the wind-oriented velocities:
                call getWindAxis(vAvg(:),V_wind(:),alpha,liftIndex)
                 !print *,'gm1',gm1
             
                !The variation of Entropy wrt the free stream:
                ds = (RGas/gm1)*log((Pavg/Pinf)*(rhoInf/rhoAvg)**gammaConstant)
                !print *,'RGas',rgas,gm1
                !print *,'ds',ds,(RGas/gm1),(Pavg/Pinf),(rhoInf/rhoAvg),gammaConstant
                
                ! The variation of Stagnation Enthalpy relative to free stream:
                dH = (gammaConstant/gm1)*(PAvg/rhoAvg-Pinf/rhoinf) + &
                     0.5*(( V_wind(1)**2 + V_wind(2)**2 + V_wind(3)**2) - Uinf**2)
                
                dPoP = (Pavg-Pinf)/Pinf
                !print *,'dpop',dpop
                dSoR = ds/RGas
                !print *,'dSoR',dsor
                dHou2 = dH/Uinf**2
                !print *,'dhou2',dhou2 ,dH,Uinf**2
                du_ir = uInf*(ffp1*dPoP + ffs1*dSoR + ffH1*dHou2 + &
                     ffp2*dPoP + ffs2*dSoR + ffH2*dHou2 + &
                     ffps2*dPoP*dSoR + ffph2*dPoP*dHou2 + ffsh2*dSoR*dHou2)
               ! print *,'du_ir', du_ir, uInf,ffp1*dPoP ,ffs1*dSoR, ffH1*dHou2,&
               !      ffp2*dPoP,ffs2*dSoR, ffH2*dHou2, &
               !      ffps2*dPoP*dSoR,ffph2*dPoP*dHou2, ffsh2*dSoR*dHou2

                fi = -rhoAvg*(v_wind(1)-uInf-du_ir)*v_wind(:)-(pAvg-pInf)*norm(i,j,:)
                !print *,'fi',fi,-rhoAvg,v_wind(1),uInf,du_ir,'v',v_wind(:),pAvg,pInf,norm(i,j,:)
                do m = 1,3
                   do k=1,3
                      drag_local(m) = drag_local(m)+fi(k)*ss(i,j,k)!*norm(i,j,k)
                   end do
                end do
                
              
            enddo
         enddo


      endif farfieldForce

   enddo bocos
end do domains
end do spectralLoop

  ! Reduce the drag to root proc:
  call mpi_reduce(drag_local,drag,1,sumb_real,mpi_sum,0,SUmb_comm_world,ierr)
  
  fact = two/(gammaInf*pInf*MachCoef*MachCoef &
            *surfaceRef*LRef*LRef)
  a = sqrt(gammaInf*pRef/rhoRef)
  ovfact = 0.5*rhoref*(machCoef*a)**2*surfaceRef
  if (myid == 0) then
     print *,'fact',fact,ovfact
  end if
  !reduce force computation
  call mpi_reduce(force_local,force,3,sumb_real,mpi_sum,0,SUmb_comm_world,ierr)
  
  !print *,'fact',fact
  !print *,'factff',fact,two,gammaInf,pInf,pRef,MachCoef,MachCoef,surfaceRef,LRef,LRef
  CL = force(1)*liftDirection(1) &
       + force(2)*liftDirection(2) &
       + force(3)*liftDirection(3)

  CD = force(1)*dragDirection(1) &
       + force(2)*dragDirection(2) &
       + force(3)*dragDirection(3)
  
  if (myid == 0) then
     print *,'Induced drag:',drag
     print *,'Induced CD',drag*fact
     print *,'force',force,force*fact,force/ovfact
     print *,'force CL',CL,CL*fact,CL/ovfact
     print *,'force CD',CD,CD*fact,CD/ovfact
  end if
end subroutine farFieldInducedDrag
