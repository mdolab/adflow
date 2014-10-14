!      ==================================================================

subroutine computedwDADI
  !
  !      ******************************************************************
  !      *                                                                *
  !      * executeRkStage executes one runge kutta stage. The stage       *
  !      * number, rkStage, is defined in the local module iteration.     *
  !      *                                                                *
  !      ******************************************************************
  !
  use blockPointers
  use constants
  use flowVarRefState
  use inputIteration
  use inputPhysics
  use inputTimeSpectral
  use inputUnsteady
  use iteration
  implicit none
  !
  !      Local parameter.
  !
  !      Local variables.
  !
  integer(kind=intType) :: i, j, k, n
  real(kind=realType) :: gm1, epsval, fac, eps2
  real(kind=realType) :: uvel, vvel, wvel, cijk, cijkinv, c2inv, uvw
  real(kind=realType) :: ri1, ri2, ri3, rj1, rj2, rj3, rk1, rk2, rk3, uu
  real(kind=realType) :: ri, rj, rk, qsi, qsj, qsk, currentCfl
  real(kind=realType) :: xfact,  cInf, cInf2
  real(kind=realType) :: dw1, dw2, dw3, dw4, dw5
  real(kind=realType) :: a1, a2, a3, a4, a5, a6, a7, mut, ge
  real(kind=realType) :: viscTerm1, viscTerm2, viscTerm3
  real(kind=realType) :: metterm
  real(kind=realType) :: unsteadyImpl, mult
  real(kind=realType) :: sqrt2, sqrt2inv, alph, alphinv
  real(kind=realType) :: volhalf, volhalfrho, volfact
  real(kind=realType) :: mutpI, mutmI, mutpJ, mutmJ, mutpK, mutmK, &
       mettermp,mettermm, &
       viscTermI,viscTermJ,viscTermK

  real(kind=realType), dimension(:,:,:), pointer:: qq_i,qq_j,qq_k, &
       cc_i,cc_j,cc_k,spectral_i,spectral_j,spectral_k,dual_dt
  real(kind=realType), dimension(ie,5) :: bbi,cci,ddi,ffi
  real(kind=realType), dimension(je,5) :: bbj,ccj,ddj,ffj
  real(kind=realType), dimension(ke,5) :: bbk,cck,ddk,ffk
  real(kind=realType), dimension(ie) :: mettermi
  real(kind=realType), dimension(je) :: mettermj
  real(kind=realType), dimension(ke) :: mettermk
  real(kind=realType), dimension(5) :: diagPlus,diagMinus
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************

  ! Set the value of the current cfl number,
  ! depending on the situation. On the finest grid in the mg cycle
  ! the second halo is computed, otherwise not.

  currentCfl = cflCoarse
  if (currentLevel == 1) then
     currentCfl = cfl
  end if
  qq_i => dadidata(:,:,:,1)
  qq_j => dadidata(:,:,:,2)
  qq_k => dadidata(:,:,:,3)
  cc_i => dadidata(:,:,:,4)
  cc_j => dadidata(:,:,:,5)
  cc_k => dadidata(:,:,:,6)
  spectral_i => dadidata(:,:,:,7)
  spectral_j => dadidata(:,:,:,8)
  spectral_k => dadidata(:,:,:,9)
  dual_dt => dadidata(:,:,:,10)

  !  havent thought about iblank

  !   these are factors for robustness. 

  epsval = 0.08
  fac    = 1.05
  mut    = zero
  cInf2  = gammaInf*pInf/rhoInf
  cInf   = sqrt(cInf2)


  qsi = zero
  qsj = zero
  qsk = zero
  viscTermi=zero
  viscTermj=zero
  viscTermk=zero
  sqrt2=sqrt(two)
  sqrt2inv=one/sqrt2


  if(equationMode.eq.steady) then
     do k=2,kl
        do j=2,jl
           do i=2,il
              dual_dt(i,j,k)= currentCfl*dtl(i,j,k)*vol(i,j,k)
           enddo
        enddo
     enddo
  else  
     do k=2,kl
        do j=2,jl
           do i=2,il
              unsteadyImpl = coefTime(0)*timeRef/deltaT
              mult = currentCfl*dtl(i,j,k)
              mult = mult/(mult*unsteadyImpl*vol(i,j,k) + one)
              dual_dt(i,j,k)   = mult*vol(i,j,k)
           enddo
        enddo
     enddo
  endif

  !     Set up some arrays
  do k=2,kl
     do j=2,jl
        do i=2,il
           volhalf    = half/vol(i,j,k)
           volhalfrho = volhalf/w(i,j,k,irho)
           uvel  = w(i,j,k,ivx)
           vvel  = w(i,j,k,ivy)
           wvel  = w(i,j,k,ivz)

           cijk  = sqrt(gamma(i,j,k)*p(i,j,k)/w(i,j,k,irho))

           ri1  = volhalf*(si(i,j,k,1) + si(i-1,j,k,1))
           ri2  = volhalf*(si(i,j,k,2) + si(i-1,j,k,2))
           ri3  = volhalf*(si(i,j,k,3) + si(i-1,j,k,3))

           rj1  = volhalf*(sj(i,j,k,1) + sj(i,j-1,k,1))
           rj2  = volhalf*(sj(i,j,k,2) + sj(i,j-1,k,2))
           rj3  = volhalf*(sj(i,j,k,3) + sj(i,j-1,k,3))

           rk1  = volhalf*(sk(i,j,k,1) + sk(i,j,k-1,1))
           rk2  = volhalf*(sk(i,j,k,2) + sk(i,j,k-1,2))
           rk3  = volhalf*(sk(i,j,k,3) + sk(i,j,k-1,3))

           if( addGridVelocities ) then
              qsi = (sFaceI(i-1,j,k) + sFaceI(i,j,k))*volhalf
              qsj = (sFaceJ(i,j-1,k) + sFaceJ(i,j,k))*volhalf
              qsk = (sFaceK(i,j,k-1) + sFaceK(i,j,k))*volhalf
           endif


           qq_i(i,j,k) = ri1*w(i,j,k,ivx) + ri2*w(i,j,k,ivy) +  &
                ri3*w(i,j,k,ivz) - qsi
           qq_j(i,j,k) = rj1*w(i,j,k,ivx) + rj2*w(i,j,k,ivy) +  &
                rj3*w(i,j,k,ivz) - qsj
           qq_k(i,j,k) = rk1*w(i,j,k,ivx) + rk2*w(i,j,k,ivy) +  &
                rk3*w(i,j,k,ivz) - qsk

           ri   = sqrt(ri1*ri1+ri2*ri2+ri3*ri3)
           rj   = sqrt(rj1*rj1+rj2*rj2+rj3*rj3)
           rk   = sqrt(rk1*rk1+rk2*rk2+rk3*rk3)

           cc_i(i,j,k) = cijk*ri
           cc_j(i,j,k) = cijk*rj
           cc_k(i,j,k) = cijk*rk
           if( viscous )  then
              mutpI = rlv(i,j,k)+rlv(i+1,j,k)
              mutmI = rlv(i,j,k)+rlv(i-1,j,k)
              mutpJ = rlv(i,j,k)+rlv(i,j+1,k)
              mutmJ = rlv(i,j,k)+rlv(i,j-1,k)
              mutpK = rlv(i,j,k)+rlv(i,j,k+1)
              mutmK = rlv(i,j,k)+rlv(i,j,k-1)
      	      if( eddyModel ) then
       		 mutpI = mutpI+rev(i,j,k)+rev(i+1,j,k)
       		 mutmI = mutpI+rev(i,j,k)+rev(i-1,j,k)
       		 mutpJ = mutpJ+rev(i,j,k)+rev(i,j+1,k)
       		 mutmJ = mutpJ+rev(i,j,k)+rev(i,j-1,k)
       		 mutpK = mutpK+rev(i,j,k)+rev(i,j,k+1)
       		 mutmK = mutpK+rev(i,j,k)+rev(i,j,k-1)
      	      endif

              volfact  =   two/(vol(i,j,k)+vol(i+1,j,k))
              mettermp =   (si(i,j,k,1)*si(i,j,k,1) &
                   +si(i,j,k,2)*si(i,j,k,2) &
                   +si(i,j,k,3)*si(i,j,k,3))*mutpI*volfact
              volfact  =   two/(vol(i,j,k)+vol(i-1,j,k))
              mettermm =   (si(i-1,j,k,1)*si(i-1,j,k,1) &
                   +si(i-1,j,k,2)*si(i-1,j,k,2) &
                   +si(i-1,j,k,3)*si(i-1,j,k,3))*mutmI*volfact
              viscTermi=   (mettermp+mettermm)*volhalfrho

              volfact  =   two/(vol(i,j,k)+vol(i,j+1,k))
              mettermp =   (sj(i,j,k,1)*sj(i,j,k,1) &
                   +sj(i,j,k,2)*sj(i,j,k,2) &
                   +sj(i,j,k,3)*sj(i,j,k,3))*mutpJ*volfact
              volfact  =   two/(vol(i,j,k)+vol(i,j-1,k))
              mettermm =   (sj(i,j-1,k,1)*sj(i,j-1,k,1) &
                   +sj(i,j-1,k,2)*sj(i,j-1,k,2) &
                   +sj(i,j-1,k,3)*sj(i,j-1,k,3))*mutmJ*volfact
              viscTermj=   (mettermp+mettermm)*volhalfrho

              volfact  =   two/(vol(i,j,k)+vol(i,j,k+1))
              mettermp =   (sk(i,j,k,1)*sk(i,j,k,1) &
                   +sk(i,j,k,2)*sk(i,j,k,2) &
                   +sk(i,j,k,3)*sk(i,j,k,3))*mutpK*volfact
              volfact  =   two/(vol(i,j,k)+vol(i,j,k-1))
              mettermm =   (sk(i,j,k-1,1)*sk(i,j,k-1,1) &
                   +sk(i,j,k-1,2)*sk(i,j,k-1,2) &
                   +sk(i,j,k-1,3)*sk(i,j,k-1,3))*mutmK*volfact
              viscTermk=   (mettermp+mettermm)*volhalfrho

           endif

           eps2 = ri*cInf*epsval
           spectral_i(i,j,k) = (fac*sqrt((abs(qq_i(i,j,k))+cc_i(i,j,k))**2 &
                + eps2**2)+viscTermi)*dual_dt(i,j,k)
           eps2 = rj*cInf*epsval
           spectral_j(i,j,k) = (fac*sqrt((abs(qq_j(i,j,k))+cc_j(i,j,k))**2 &
                + eps2**2)+viscTermj)*dual_dt(i,j,k)
           eps2 = rk*cInf*epsval
           spectral_k(i,j,k) = (fac*sqrt((abs(qq_k(i,j,k))+cc_k(i,j,k))**2 & 
                + eps2**2)+viscTermk)*dual_dt(i,j,k)
           spectral_i(i,j,k)=spectral_i(i,j,k)*zero
           spectral_j(i,j,k)=spectral_j(i,j,k)*zero
           spectral_k(i,j,k)=spectral_k(i,j,k)*zero
        enddo
     enddo
  enddo

  !       Multiply by T_eta^inv
  do k=2,kl
     do j=2,jl
        do i=2,il

           gm1   = gamma(i,j,k)-one
           cijk  = sqrt(gamma(i,j,k)*p(i,j,k)/w(i,j,k,irho))
           c2inv = one/(cijk*cijk)
           xfact = two*cijk
           alphinv  = sqrt2*cijk/w(i,j,k,irho)

           uvel  = w(i,j,k,ivx)
           vvel  = w(i,j,k,ivy)
           wvel  = w(i,j,k,ivz)
           uvw   = half*(uvel*uvel+vvel*vvel+wvel*wvel)

           rj1  = half*(sj(i,j,k,1) + sj(i,j-1,k,1))
           rj2  = half*(sj(i,j,k,2) + sj(i,j-1,k,2))
           rj3  = half*(sj(i,j,k,3) + sj(i,j-1,k,3))
           rj   = sqrt(rj1*rj1+rj2*rj2+rj3*rj3)
           uu   = uvel*rj1+vvel*rj2+wvel*rj3
           rj1  = rj1/rj
           rj2  = rj2/rj
           rj3  = rj3/rj

           dw1   = dw(i,j,k,1)
           dw2   = dw(i,j,k,2)
           dw3   = dw(i,j,k,3)
           dw4   = dw(i,j,k,4)
           dw5   = dw(i,j,k,5)

           a1    =  dw2*uvel+dw3*vvel+dw4*wvel-dw5
           a1    =  a1*gm1*c2inv+dw1*(one-uvw*gm1*c2inv)

           a2    =  (rj2*wvel-rj3*vvel)*dw1+rj3*dw3-rj2*dw4
           a3    =  (rj3*uvel-rj1*wvel)*dw1+rj1*dw4-rj3*dw2
           a4    =  (rj1*vvel-rj2*uvel)*dw1+rj2*dw2-rj1*dw3

           a5    =  uvw*dw1-uvel*dw2-vvel*dw3-wvel*dw4+dw5
           a5    =  a5*gm1*c2inv

           a6    =  uu*dw1/rj-rj1*dw2-rj2*dw3-rj3*dw4

           dw(i,j,k,1)  = a1*rj1+a2/w(i,j,k,irho)
           dw(i,j,k,2)  = a1*rj2+a3/w(i,j,k,irho)
           dw(i,j,k,3)  = a1*rj3+a4/w(i,j,k,irho)
           dw(i,j,k,4)  = (half*a5-a6/xfact)*alphinv
           dw(i,j,k,5)  = (half*a5+a6/xfact)*alphinv

        enddo
     enddo
  enddo

  if(jl.gt.2) then

     !  Inversion in j

     do k=2,kl
        do i=2,il
           do j=1,jl
              if( viscous )   mut = rlv(i,j,k)+rlv(i,j+1,k)
              if( eddyModel ) mut = mut+rev(i,j,k)+rev(i,j+1,k)

              volfact=one/(vol(i,j,k)+vol(i,j+1,k))
              metterm =   (sj(i,j,k,1)*sj(i,j,k,1) &
                   +sj(i,j,k,2)*sj(i,j,k,2) &
                   +sj(i,j,k,3)*sj(i,j,k,3))
              mettermj(j)= metterm*mut*volfact
           enddo

           do j=2,jl

              viscTerm1 = mettermj(j)  /vol(i,j,k)/w(i,j,k,irho)
              viscTerm3 = mettermj(j-1)/vol(i,j,k)/w(i,j,k,irho)
              viscTerm2 = viscTerm1+viscTerm3

              volhalf = half/vol(i,j,k)
              rj1  = volhalf*(sj(i,j,k,1) + sj(i,j-1,k,1))
              rj2  = volhalf*(sj(i,j,k,2) + sj(i,j-1,k,2))
              rj3  = volhalf*(sj(i,j,k,3) + sj(i,j-1,k,3))
              metterm     = rj1*rj1+rj2*rj2+rj3*rj3
              eps2        = epsval*epsval*cInf2*metterm
              diagPlus(1) =half*(qq_j(i,j,k)+fac*sqrt(qq_j(i,j,k)**2+eps2))
              diagPlus(2) =diagPlus(1)
              diagPlus(3) =diagPlus(1)
              diagPlus(4) =half*(qq_j(i,j,k)+cc_j(i,j,k)   &
                   +fac*sqrt((qq_j(i,j,k)+cc_j(i,j,k))**2+eps2))
              diagPlus(5) =half*(qq_j(i,j,k)-cc_j(i,j,k)   &
                   +fac*sqrt((qq_j(i,j,k)-cc_j(i,j,k))**2+eps2))
              diagMinus(1) =half*(qq_j(i,j,k)-fac*sqrt(qq_j(i,j,k)**2+eps2))
              diagMinus(2) =diagMinus(1)
              diagMinus(3) =diagMinus(1)
              diagMinus(4) =half*(qq_j(i,j,k)+cc_j(i,j,k)   &
                   -fac*sqrt((qq_j(i,j,k)+cc_j(i,j,k))**2+eps2))
              diagMinus(5) =half*(qq_j(i,j,k)-cc_j(i,j,k)   &
                   -fac*sqrt((qq_j(i,j,k)-cc_j(i,j,k))**2+eps2))

              do n=1,5
                 bbj(j+1,n)= -viscTerm1-diagPlus(n)
                 ddj(j-1,n)= -viscTerm3+diagMinus(n)
                 ccj(j  ,n)=  viscTerm2+diagPlus(n)-diagMinus(n) 
              enddo
           enddo

           do n=1,5
              bbj(je,n)=zero
              ddj(1 ,n)=zero
              do j=2,jl
                 bbj(j,n)=bbj(j,n)*dual_dt(i,j,k)
                 ddj(j,n)=ddj(j,n)*dual_dt(i,j,k)
                 ccj(j,n)=one+ccj(j,n)*dual_dt(i,j,k)+spectral_i(i,j,k)+spectral_k(i,j,k)
                 ffj(j,n)=dw(i,j,k,n)
              enddo
           enddo

           call tridiagsolve(bbj,ccj,ddj,ffj,jl)

           do n=1,5
              do j=2,jl
                 dw(i,j,k,n)=ffj(j,n)
              enddo
           enddo

        enddo
     enddo

  endif
  !       Multiply by T_xi^inv T_eta


  do k=2,kl 
     do j=2,jl
        do i=2,il
           ri1  = half*(si(i,j,k,1) + si(i-1,j,k,1))
           ri2  = half*(si(i,j,k,2) + si(i-1,j,k,2))
           ri3  = half*(si(i,j,k,3) + si(i-1,j,k,3))
           ri   = sqrt(ri1*ri1+ri2*ri2+ri3*ri3)
           ri1  = ri1/ri
           ri2  = ri2/ri
           ri3  = ri3/ri

           rj1  = half*(sj(i,j,k,1) + sj(i,j-1,k,1))
           rj2  = half*(sj(i,j,k,2) + sj(i,j-1,k,2))
           rj3  = half*(sj(i,j,k,3) + sj(i,j-1,k,3))
           rj   = sqrt(rj1*rj1+rj2*rj2+rj3*rj3)
           rj1  = rj1/rj
           rj2  = rj2/rj
           rj3  = rj3/rj

           dw1   = dw(i,j,k,1)
           dw2   = dw(i,j,k,2)
           dw3   = dw(i,j,k,3)
           dw4   = dw(i,j,k,4)
           dw5   = dw(i,j,k,5)

           a1 = ( ri1*rj1 + ri2*rj2 + ri3*rj3 )
           a2 = ( ri1*rj2 - rj1*ri2 )
           a3 = ( ri3*rj2 - rj3*ri2 )
           a4 = ( ri1*rj3 - rj1*ri3 )
           a5 = (dw4-dw5)*sqrt2inv
           a6 = (dw4+dw5)*half
           a7 = (a3*dw1+a4*dw2-a2*dw3-a5*a1)*sqrt2inv

           dw(i,j,k,1)=  a1*dw1 + a2*dw2 + a4*dw3 + a5*a3
           dw(i,j,k,2)=- a2*dw1 + a1*dw2 - a3*dw3 + a5*a4
           dw(i,j,k,3)=- a4*dw1 + a3*dw2 + a1*dw3 - a5*a2
           dw(i,j,k,4)= -a7+a6
           dw(i,j,k,5)=  a7+a6

        enddo
     enddo
  enddo

  !  Multiply by diagonal

  do k=2,kl
     do j=2,jl
        do i=2,il
           xfact = one+spectral_i(i,j,k)+spectral_j(i,j,k) &
                +spectral_k(i,j,k)	
           dw(i,j,k,1)=dw(i,j,k,1)*xfact
           dw(i,j,k,2)=dw(i,j,k,2)*xfact
           dw(i,j,k,3)=dw(i,j,k,3)*xfact
           dw(i,j,k,4)=dw(i,j,k,4)*xfact
           dw(i,j,k,5)=dw(i,j,k,5)*xfact
        enddo
     enddo
  enddo



  if(il.gt.2) then

     !  Inversion in i

     do k=2,kl
        do j=2,jl
           do i=1,il
              if( viscous )   mut = rlv(i,j,k)+rlv(i+1,j,k)
              if( eddyModel ) mut = mut+rev(i,j,k)+rev(i+1,j,k)

              volfact=one/(vol(i,j,k)+vol(i+1,j,k))
              metterm    =(si(i,j,k,1)*si(i,j,k,1) &
                   +si(i,j,k,2)*si(i,j,k,2) &
                   +si(i,j,k,3)*si(i,j,k,3))
              mettermi(i)=metterm*mut*volfact
           enddo

           do i=2,il

              viscTerm1 = mettermi(i)  /vol(i,j,k)/w(i,j,k,irho)
              viscTerm3 = mettermi(i-1)/vol(i,j,k)/w(i,j,k,irho)
              viscTerm2 = viscTerm1+viscTerm3

              volhalf = half/vol(i,j,k)
              ri1  = volhalf*(si(i,j,k,1) + si(i-1,j,k,1))
              ri2  = volhalf*(si(i,j,k,2) + si(i-1,j,k,2))
              ri3  = volhalf*(si(i,j,k,3) + si(i-1,j,k,3))
              metterm     = ri1*ri1+ri2*ri2+ri3*ri3
              eps2        = epsval*epsval*cInf2*metterm
              diagPlus(1) =half*(qq_i(i,j,k)+fac*sqrt(qq_i(i,j,k)**2+eps2))
              diagPlus(2) =diagPlus(1)
              diagPlus(3) =diagPlus(1)
              diagPlus(4) =half*(qq_i(i,j,k)+cc_i(i,j,k)   &
                   +fac*sqrt((qq_i(i,j,k)+cc_i(i,j,k))**2+eps2))
              diagPlus(5) =half*(qq_i(i,j,k)-cc_i(i,j,k)   &
                   +fac*sqrt((qq_i(i,j,k)-cc_i(i,j,k))**2+eps2))
              diagMinus(1) =half*(qq_i(i,j,k)-fac*sqrt(qq_i(i,j,k)**2+eps2))
              diagMinus(2) =diagMinus(1)
              diagMinus(3) =diagMinus(1)
              diagMinus(4) =half*(qq_i(i,j,k)+cc_i(i,j,k)   &
                   -fac*sqrt((qq_i(i,j,k)+cc_i(i,j,k))**2+eps2))
              diagMinus(5) =half*(qq_i(i,j,k)-cc_i(i,j,k)   &
                   -fac*sqrt((qq_i(i,j,k)-cc_i(i,j,k))**2+eps2))
              do n=1,5
                 bbi(i+1,n)= -viscTerm1-diagPlus(n)
                 ddi(i-1,n)= -viscTerm3+diagMinus(n)
                 cci(i  ,n)=  viscTerm2+diagPlus(n)-diagMinus(n)
              enddo
           enddo

           do n=1,5
              bbi(ie ,n)=zero
              ddi(1 ,n)=zero
              do i=2,il
                 bbi(i,n)=bbi(i,n)*dual_dt(i,j,k)
                 ddi(i,n)=ddi(i,n)*dual_dt(i,j,k)
                 cci(i,n)=one+cci(i,n)*dual_dt(i,j,k)+spectral_j(i,j,k)+spectral_k(i,j,k)
                 ffi(i,n)=dw(i,j,k,n)
              enddo
           enddo


           call tridiagsolve(bbi,cci,ddi,ffi,il)

           do n=1,5
              do i=2,il
                 dw(i,j,k,n)=ffi(i,n)
              enddo
           enddo

        enddo
     enddo

  endif
  !       Multiply by T_zeta^inv T_xi


  do k=2,kl 
     do j=2,jl
        do i=2,il
           ri1  = half*(si(i,j,k,1) + si(i-1,j,k,1))
           ri2  = half*(si(i,j,k,2) + si(i-1,j,k,2))
           ri3  = half*(si(i,j,k,3) + si(i-1,j,k,3))
           ri   = sqrt(ri1*ri1+ri2*ri2+ri3*ri3)
           ri1  = ri1/ri
           ri2  = ri2/ri
           ri3  = ri3/ri

           rk1  = half*(sk(i,j,k,1) + sk(i,j,k-1,1))
           rk2  = half*(sk(i,j,k,2) + sk(i,j,k-1,2))
           rk3  = half*(sk(i,j,k,3) + sk(i,j,k-1,3))
           rk   = sqrt(rk1*rk1+rk2*rk2+rk3*rk3)
           rk1  = rk1/rk
           rk2  = rk2/rk
           rk3  = rk3/rk

           dw1   = dw(i,j,k,1)
           dw2   = dw(i,j,k,2)
           dw3   = dw(i,j,k,3)
           dw4   = dw(i,j,k,4)
           dw5   = dw(i,j,k,5)

           a1 =  ri1*rk1 + ri2*rk2 + ri3*rk3 
           a2 =  rk1*ri2 - ri1*rk2 
           a3 =  rk3*ri2 - ri3*rk2 
           a4 =  rk1*ri3 - ri1*rk3 
           a5 = (dw4-dw5)*sqrt2inv
           a6 = (dw4+dw5)*half
           a7 = (a3*dw1+a4*dw2-a2*dw3-a5*a1)*sqrt2inv

           dw(i,j,k,1)=  a1*dw1 + a2*dw2 + a4*dw3 + a5*a3
           dw(i,j,k,2)=- a2*dw1 + a1*dw2 - a3*dw3 + a5*a4
           dw(i,j,k,3)=- a4*dw1 + a3*dw2 + a1*dw3 - a5*a2
           dw(i,j,k,4)= -a7+a6
           dw(i,j,k,5)=  a7+a6
        enddo
     enddo
  enddo

  !  Multiply by diagonal

  do k=2,kl
     do j=2,jl
        do i=2,il
           xfact = one+spectral_i(i,j,k)+spectral_j(i,j,k) &
                +spectral_k(i,j,k)	
           dw(i,j,k,1)=dw(i,j,k,1)*xfact
           dw(i,j,k,2)=dw(i,j,k,2)*xfact
           dw(i,j,k,3)=dw(i,j,k,3)*xfact
           dw(i,j,k,4)=dw(i,j,k,4)*xfact
           dw(i,j,k,5)=dw(i,j,k,5)*xfact
        enddo
     enddo
  enddo

  if(kl.gt.2) then

     !  Inversion in k

     do j=2,jl
        do i=2,il
           do k=1,kl
              if( viscous )   mut = rlv(i,j,k)+rlv(i,j,k+1)
              if( eddyModel ) mut = mut+rev(i,j,k)+rev(i,j,k+1)

              volfact=1./(vol(i,j,k)+vol(i,j,k+1))
              metterm    = sk(i,j,k,1)*sk(i,j,k,1)  &
                   +sk(i,j,k,2)*sk(i,j,k,2)  &
                   +sk(i,j,k,3)*sk(i,j,k,3)
              mettermk(k)= metterm*mut*volfact
           enddo

           do k=2,kl

              viscTerm1 = mettermk(k)  /vol(i,j,k)/w(i,j,k,irho)
              viscTerm3 = mettermk(k-1)/vol(i,j,k)/w(i,j,k,irho)
              viscTerm2 = viscTerm1+viscTerm3

              volhalf = half/vol(i,j,k)
              rk1  = volhalf*(sk(i,j,k,1) + sj(i,j,k-1,1))
              rk2  = volhalf*(sk(i,j,k,2) + sj(i,j,k-1,2))
              rk3  = volhalf*(sk(i,j,k,3) + sj(i,j,k-1,3))
              metterm     = rk1*rk1+rk2*rk2+rk3*rk3
              eps2        = epsval*epsval*cInf2*metterm
              diagPlus(1) =half*(qq_k(i,j,k)+fac*sqrt(qq_k(i,j,k)**2+eps2))
              diagPlus(2) =diagPlus(1)
              diagPlus(3) =diagPlus(1)
              diagPlus(4) =half*(qq_k(i,j,k)+cc_k(i,j,k)   &
                   +fac*sqrt((qq_k(i,j,k)+cc_k(i,j,k))**2+eps2))
              diagPlus(5) =half*(qq_k(i,j,k)-cc_k(i,j,k)   &
                   +fac*sqrt((qq_k(i,j,k)-cc_k(i,j,k))**2+eps2))
              diagMinus(1) =half*(qq_k(i,j,k)-fac*sqrt(qq_k(i,j,k)**2+eps2))
              diagMinus(2) =diagMinus(1)
              diagMinus(3) =diagMinus(1)
              diagMinus(4) =half*(qq_k(i,j,k)+cc_k(i,j,k)   &
                   -fac*sqrt((qq_k(i,j,k)+cc_k(i,j,k))**2+eps2))
              diagMinus(5) =half*(qq_k(i,j,k)-cc_k(i,j,k)   &
                   -fac*sqrt((qq_k(i,j,k)-cc_k(i,j,k))**2+eps2))

              do n=1,5
                 bbk(k+1,n)= -viscTerm1-diagPlus(n)
                 ddk(k-1,n)= -viscTerm3+diagMinus(n)
                 cck(k  ,n)=  viscTerm2+diagPlus(n)-diagMinus(n)
              enddo
           enddo

           do n=1,5
              bbk(ke,n)=zero
              ddk(1 ,n)=zero
              do k=2,kl
                 bbk(k,n)=bbk(k,n)*dual_dt(i,j,k)
                 ddk(k,n)=ddk(k,n)*dual_dt(i,j,k)
                 cck(k,n)=one+cck(k,n)*dual_dt(i,j,k)+spectral_i(i,j,k)+spectral_j(i,j,k)
                 ffk(k,n)=dw(i,j,k,n)
              enddo
           enddo

           call tridiagsolve(bbk,cck,ddk,ffk,kl)

           do n=1,5
              do k=2,kl
                 dw(i,j,k,n)=ffk(k,n)
              enddo
           enddo

        enddo
     enddo
  endif

  !       Multiply by T_zeta

  do k=2,kl
     do j=2,jl
        do i=2,il
           uvel  = w(i,j,k,ivx)
           vvel  = w(i,j,k,ivy)
           wvel  = w(i,j,k,ivz)

           rk1  = half*(sk(i,j,k,1) + sk(i,j,k-1,1))
           rk2  = half*(sk(i,j,k,2) + sk(i,j,k-1,2))
           rk3  = half*(sk(i,j,k,3) + sk(i,j,k-1,3))

           rk = sqrt(rk1*rk1+rk2*rk2+rk3*rk3)
           uu = uvel*rk1+vvel*rk2+wvel*rk3

           rk1  = rk1/rk
           rk2  = rk2/rk
           rk3  = rk3/rk

           uvw   = half*(uvel*uvel+vvel*vvel+wvel*wvel)
           cijkinv = sqrt(w(i,j,k,irho)/gamma(i,j,k)/p(i,j,k))
           alph  = w(i,j,k,irho)*cijkinv*sqrt2inv
           xfact = two/cijkinv

           ge=gamma(i,j,k)*w(i,j,k,irhoE)/w(i,j,k,irho)-  &
                (gamma(i,j,k)-one)*uvw  

           dw1   = dw(i,j,k,1)
           dw2   = dw(i,j,k,2)
           dw3   = dw(i,j,k,3)
           dw4   = dw(i,j,k,4)*alph
           dw5   = dw(i,j,k,5)*alph

           a1    =  dw1*rk1+dw2*rk2+dw3*rk3+dw4+dw5
           a2    =  half*xfact*(dw4-dw5)
           a3    =  uvw*(rk1*dw1+rk2*dw2+rk3*dw3)

           dw(i,j,k,1)  = a1
           dw(i,j,k,2)  = a1*uvel-w(i,j,k,irho)*(rk3*dw2-rk2*dw3) & 
                +a2*rk1
           dw(i,j,k,3)  = a1*vvel-w(i,j,k,irho)*(rk1*dw3-rk3*dw1) & 
                +a2*rk2
           dw(i,j,k,4)  = a1*wvel-w(i,j,k,irho)*(rk2*dw1-rk1*dw2) & 
                +a2*rk3
           dw(i,j,k,5)  = a3+ w(i,j,k,irho)*       &
                ((vvel*rk3-wvel*rk2)*dw1            &   
                +(wvel*rk1-uvel*rk3)*dw2            &         
                +(uvel*rk2-vvel*rk1)*dw3)           & 
                +(ge+half*xfact*uu/rk)*dw4 &   
                +(ge-half*xfact*uu/rk)*dw5       

        enddo
     enddo
  enddo


  !      For consistency with update.

  do k=2,kl
     do j=2,jl
        do i=2,il
           volfact=-one/vol(i,j,k)
           dw(i,j,k,1)=dw(i,j,k,1)*volfact
           dw(i,j,k,2)=dw(i,j,k,2)*volfact
           dw(i,j,k,3)=dw(i,j,k,3)*volfact
           dw(i,j,k,4)=dw(i,j,k,4)*volfact
           dw(i,j,k,5)=dw(i,j,k,5)*volfact
        enddo
     enddo
  enddo


end subroutine computedwDADI

subroutine tridiagsolve(bb,cc,dd,ff,nn)
  use precision
  implicit none
  !
  !      Subroutine arguments
  !
  integer(kind=intType) :: nn
  real(kind=realType), dimension(nn+1,5) :: bb,cc,dd,ff

  !     local variables

  integer(kind=intType) :: m,n
  real(kind=realType) :: d0,d2

  do n=1,5
     m=2
     d0=1./cc(m,n)
     dd(m,n)=dd(m,n)*d0
     ff(m,n)=ff(m,n)*d0


     do m=3,nn
        d2=bb(m,n)
        d0=1./(cc(m,n)-d2*dd(m-1,n))
        ff(m,n)=(ff(m,n)-d2*ff(m-1,n))*d0
        dd(m,n)=dd(m,n)*d0
     enddo


     do m=nn-1,2,-1
        ff(m,n)=ff(m,n)-dd(m,n)*ff(m+1,n)
     enddo

  enddo


end subroutine tridiagsolve
