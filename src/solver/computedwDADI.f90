!      ==================================================================

subroutine computedwDADI
  !
  !       Executes one DADI Stage                                        
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
  real(kind=realType) :: gm1,epsval,fac,eps2,epscinf2,nhat
  real(kind=realType) :: uvel,vvel,wvel,cijk,cijkinv,c2inv,cijk2
  real(kind=realType) :: ri1,ri2,ri3,rj1,rj2,rj3,rk1,rk2,rk3
  real(kind=realType) :: el1,el2,el3,el,em1,em2,em3,em
  real(kind=realType) :: el1i,el2i,el3i,em1i,em2i,em3i,eli
  real(kind=realType) :: el1j,el2j,el3j,em1j,em2j,em3j,elj
  real(kind=realType) :: el1k,el2k,el3k,em1k,em2k,em3k,elk
  real(kind=realType) :: ri,rj,rk,qsi,qsj,qsk,currentCfl,uvw
  real(kind=realType) :: xfact,xfacti,xfactj,xfactk,cInf,cInf2
  real(kind=realType) :: dw1,dw2,dw3,dw4,dw5,uu,vv,ww,mut
  real(kind=realType) :: a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11
  real(kind=realType) :: unsteadyImpl,mult,volfact,rhoinv
  real(kind=realType) :: uup,uum,cijkp,cijkm

  real(kind=realType), dimension(ie,je,ke) :: qq_i,qq_j,qq_k, &
       cc_i,cc_j,cc_k,dual_dt, &
       viscTerm_I,viscTerm_J,viscTerm_K
  real(kind=realType), dimension(ie,5) :: bbi,cci,ddi,ffi
  real(kind=realType), dimension(je,5) :: bbj,ccj,ddj,ffj
  real(kind=realType), dimension(ke,5) :: bbk,cck,ddk,ffk
  real(kind=realType), dimension(5) :: diagPlus,diagMinus
  !
  !       Begin execution                                                

  ! Set the value of the current cfl number,
  ! depending on the situation. On the finest grid in the mg cycle
  ! the second halo is computed, otherwise not.
  currentCfl = cflCoarse
  if (currentLevel == 1) then
     currentCfl = cfl
  end if

  !  havent thought about iblank

  !   these are factors for robustness. 

  epsval = 0.00
  fac    = 1.00
  eps2   = zero
  mut    = zero
  nhat   = one/sqrt(three)
  cInf2  = gammaInf*pInf/rhoInf
  cInf   = sqrt(cInf2)
  epscinf2 = cinf2*epsval*epsval


  qsi = zero
  qsj = zero
  qsk = zero


  if(equationMode.eq.steady) then
     do k=2,kl
        do j=2,jl
           do i=2,il
              dual_dt(i,j,k)= currentCfl*dtl(i,j,k)
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
              dual_dt(i,j,k)   = mult
           enddo
        enddo
     enddo
  endif

  !     Set up some arrays
  do k=2,kl
     do j=2,jl
        do i=2,il
           uvel  = w(i,j,k,ivx)
           vvel  = w(i,j,k,ivy)
           wvel  = w(i,j,k,ivz)
           rhoinv = one/w(i,j,k,irho)
           volfact= one/vol(i,j,k)

           cijk2 = gamma(i,j,k)*p(i,j,k)*rhoinv

           ri1  = half*(si(i,j,k,1) + si(i-1,j,k,1))
           ri2  = half*(si(i,j,k,2) + si(i-1,j,k,2))
           ri3  = half*(si(i,j,k,3) + si(i-1,j,k,3))

           rj1  = half*(sj(i,j,k,1) + sj(i,j-1,k,1))
           rj2  = half*(sj(i,j,k,2) + sj(i,j-1,k,2))
           rj3  = half*(sj(i,j,k,3) + sj(i,j-1,k,3))

           rk1  = half*(sk(i,j,k,1) + sk(i,j,k-1,1))
           rk2  = half*(sk(i,j,k,2) + sk(i,j,k-1,2))
           rk3  = half*(sk(i,j,k,3) + sk(i,j,k-1,3))

           if( addGridVelocities ) then
              qsi = (sFaceI(i-1,j,k) + sFaceI(i,j,k))*half
              qsj = (sFaceJ(i,j-1,k) + sFaceJ(i,j,k))*half
              qsk = (sFaceK(i,j,k-1) + sFaceK(i,j,k))*half
           endif


           qq_i(i,j,k) = ri1*uvel + ri2*vvel +  ri3*wvel - qsi
           qq_j(i,j,k) = rj1*uvel + rj2*vvel +  rj3*wvel - qsj
           qq_k(i,j,k) = rk1*uvel + rk2*vvel +  rk3*wvel - qsk

           ri  = ri1*ri1+ri2*ri2+ri3*ri3
           rj  = rj1*rj1+rj2*rj2+rj3*rj3
           rk  = rk1*rk1+rk2*rk2+rk3*rk3

           cc_i(i,j,k) = sqrt(cijk2*ri)
           cc_j(i,j,k) = sqrt(cijk2*rj)
           cc_k(i,j,k) = sqrt(cijk2*rk)
           if( viscous )  then
              mut = rlv(i,j,k)
              if( eddyModel ) mut=mut+rev(i,j,k)

              viscTerm_I(i,j,k) =  rI*rhoinv*volfact*mut
              viscTerm_J(i,j,k) =  rJ*rhoinv*volfact*mut
              viscTerm_K(i,j,k) =  rK*rhoinv*volfact*mut
           else
              viscTerm_I(i,j,k) =  zero
              viscTerm_J(i,j,k) =  zero
              viscTerm_K(i,j,k) =  zero
           endif

        enddo
     enddo
  enddo

  !    Multiply by T_eta^inv
  do k=2,kl
     do j=2,jl
        do i=2,il

           gm1   = gamma(i,j,k)-one
           cijk  = sqrt(gamma(i,j,k)*p(i,j,k)/w(i,j,k,irho))
           c2inv = one/(cijk*cijk)
           xfact = half*c2inv

           uvel  = w(i,j,k,ivx)
           vvel  = w(i,j,k,ivy)
           wvel  = w(i,j,k,ivz)

           rj1  = half*(sj(i,j,k,1) + sj(i,j-1,k,1))
           rj2  = half*(sj(i,j,k,2) + sj(i,j-1,k,2))
           rj3  = half*(sj(i,j,k,3) + sj(i,j-1,k,3))
           rj   = sqrt(rj1*rj1+rj2*rj2+rj3*rj3)
           if(rj.gt.zero) then
              rj1  = rj1/rj
              rj2  = rj2/rj
              rj3  = rj3/rj
           else
              rj1  = nhat
              rj2  =-nhat
              rj3  = nhat
              print*,'rj blowing up'
           endif

           el1  = rj2-rj3
           el2  = rj3-rj1
           el3  = rj1-rj2
           el   = sqrt(el1*el1+el2*el2+el3*el3)
           if(el.gt.zero) then
              el1  = el1/el
              el2  = el2/el
              el3  = el3/el
           else
              print*,'el blowing up'
           endif

           em1  = rj2*el3 - rj3*el2
           em2  = rj3*el1 - rj1*el3
           em3  = rj1*el2 - rj2*el1

           uu   = uvel*rj1+vvel*rj2+wvel*rj3
           vv   = uvel*el1+vvel*el2+wvel*el3
           ww   = uvel*em1+vvel*em2+wvel*em3

           dw1   = dw(i,j,k,1)
           dw2   = dw(i,j,k,2)
           dw3   = dw(i,j,k,3)
           dw4   = dw(i,j,k,4)
           dw5   = dw(i,j,k,5)

           uvw   = half*(uvel*uvel+vvel*vvel+wvel*wvel)
           a1    =  dw2*uvel+dw3*vvel+dw4*wvel-dw5
           a1    =  a1*gm1*c2inv+dw1*(one-uvw*gm1*c2inv)

           a2    =  -vv*dw1+el1*dw2+el2*dw3+el3*dw4

           a3    =  -ww*dw1+em1*dw2+em2*dw3+em3*dw4

           a5    =  (uvw*dw1-uvel*dw2-vvel*dw3-wvel*dw4+dw5)*gm1

           a6    =  (uu*dw1-rj1*dw2-rj2*dw3-rj3*dw4)*cijk

           dw(i,j,k,1)  = a1
           dw(i,j,k,2)  = a2
           dw(i,j,k,3)  = a3
           dw(i,j,k,4)  = (a5-a6)*xfact
           dw(i,j,k,5)  = (a5+a6)*xfact
        enddo
     enddo
  enddo


  if(jl.gt.2) then

     !  Inversion in j

     do k=2,kl
        do i=2,il
           do j=2,jl

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
                 bbj(j+1,n)= -diagPlus(n)*vol(i,j+1,k)/vol(i,j,k)
                 ddj(j-1,n)=  diagMinus(n)*vol(i,j-1,k)/vol(i,j,k)
                 !     bbj(j+1,n)= -diagPlus(n)
                 !     ddj(j-1,n)=  diagMinus(n)
                 ccj(j  ,n)=  diagPlus(n)-diagMinus(n) 
              enddo
           enddo

           do n=1,5
              bbj(je,n)=zero
              ddj(1 ,n)=zero
              do j=2,jl
                 bbj(j,n)=(bbj(j,n)-viscTerm_J(i,j,k))*dual_dt(i,j,k)
                 ddj(j,n)=(ddj(j,n)-viscTerm_J(i,j,k))*dual_dt(i,j,k)
                 ccj(j,n)=one+(ccj(j,n)+two*viscTerm_J(i,j,k))*dual_dt(i,j,k)
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

  !   Multiply  by T_xi^inv T_eta

  do k=2,kl 
     do j=2,jl
        do i=2,il

           ri1  = half*(si(i,j,k,1) + si(i-1,j,k,1))
           ri2  = half*(si(i,j,k,2) + si(i-1,j,k,2))
           ri3  = half*(si(i,j,k,3) + si(i-1,j,k,3))
           ri   = sqrt(ri1*ri1+ri2*ri2+ri3*ri3)

           if(ri.gt.zero) then
              ri1  = ri1/ri
              ri2  = ri2/ri
              ri3  = ri3/ri
           else
              ri1  = nhat
              ri2  =-nhat
              ri3  = nhat
              print*,'ri blowing up'
           endif

           el1i = ri2-ri3
           el2i = ri3-ri1
           el3i = ri1-ri2
           eli  = sqrt(el1i*el1i+el2i*el2i+el3i*el3i)
           if(eli.gt.zero) then
              el1i = el1i/eli
              el2i = el2i/eli
              el3i = el3i/eli
           else
              print*,'eli blowing up'
           endif

           em1i = ri2*el3i - ri3*el2i
           em2i = ri3*el1i - ri1*el3i
           em3i = ri1*el2i - ri2*el1i


           rj1  = half*(sj(i,j,k,1) + sj(i,j-1,k,1))
           rj2  = half*(sj(i,j,k,2) + sj(i,j-1,k,2))
           rj3  = half*(sj(i,j,k,3) + sj(i,j-1,k,3))
           rj   = sqrt(rj1*rj1+rj2*rj2+rj3*rj3)
           if(rj.gt.zero) then
              rj1  = rj1/rj
              rj2  = rj2/rj
              rj3  = rj3/rj
           else
              rj1  = nhat
              rj2  =-nhat
              rj3  = nhat
              print*,'rj blowing up'
           endif

           el1j = rj2-rj3
           el2j = rj3-rj1
           el3j = rj1-rj2
           elj  = sqrt(el1j*el1j+el2j*el2j+el3j*el3j)
           if(elj.gt.zero) then
              el1j = el1j/elj
              el2j = el2j/elj
              el3j = el3j/elj
           else
              print*,'elj blowing up'
           endif

           em1j = rj2*el3j - rj3*el2j
           em2j = rj3*el1j - rj1*el3j
           em3j = rj1*el2j - rj2*el1j

           dw2   = dw(i,j,k,2)
           dw3   = dw(i,j,k,3)
           dw4   = dw(i,j,k,4)
           dw5   = dw(i,j,k,5)

           cijk = sqrt(gamma(i,j,k)*p(i,j,k)/w(i,j,k,irho))
           xfact= half/cijk 
           a1    = dw4+dw5
           a2    = dw4-dw5

           a3    = el1j*el1i+el2j*el2i+el3j*el3i
           a4    = em1j*el1i+em2j*el2i+em3j*el3i
           a5    =  rj1*el1i+ rj2*el2i+ rj3*el3i

           a6    = el1j*em1i+el2j*em2i+el3j*em3i
           a7    = em1j*em1i+em2j*em2i+em3j*em3i
           a8    =  rj1*em1i+ rj2*em2i+ rj3*em3i

           a9   = (el1j*ri1+el2j*ri2+el3j*ri3)*xfact
           a10  = (em1j*ri1+em2j*ri2+em3j*ri3)*xfact
           a11  =   rj1*ri1+ rj2*ri2+ rj3*ri3

           dw(i,j,k,2)= a3*dw2+a4*dw3+a5*cijk*a2
           dw(i,j,k,3)= a6*dw2+a7*dw3+a8*cijk*a2
           dw(i,j,k,4)= a9*dw2+a10*dw3+half*(a1+a11*a2)
           dw(i,j,k,5)=-a9*dw2-a10*dw3+half*(a1-a11*a2)

        enddo
     enddo
  enddo

  !  Multiply by diagonal

  do k=2,kl
     do j=2,jl
        do i=2,il
           dw(i,j,k,1)=dw(i,j,k,1)
           dw(i,j,k,2)=dw(i,j,k,2)
           dw(i,j,k,3)=dw(i,j,k,3)
           dw(i,j,k,4)=dw(i,j,k,4)
           dw(i,j,k,5)=dw(i,j,k,5)
        enddo
     enddo
  enddo



  if(il.gt.2) then

     !  Inversion in i

     do k=2,kl
        do j=2,jl
           do i=2,il

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
                 bbi(i+1,n)= -diagPlus(n)*vol(i+1,j,k)/vol(i,j,k)
                 ddi(i-1,n)=  diagMinus(n)*vol(i-1,j,k)/vol(i,j,k)
                 !     bbi(i+1,n)= -diagPlus(n)
                 !     ddi(i-1,n)=  diagMinus(n)
                 cci(i  ,n)=  diagPlus(n)-diagMinus(n)
              enddo
           enddo

           do n=1,5
              bbi(ie,n)=zero
              ddi(1 ,n)=zero
              do i=2,il
                 bbi(i,n)=(bbi(i,n)-viscTerm_I(i,j,k))*dual_dt(i,j,k)
                 ddi(i,n)=(ddi(i,n)-viscTerm_I(i,j,k))*dual_dt(i,j,k)
                 cci(i,n)=one+(cci(i,n)+two*viscTerm_I(i,j,k))*dual_dt(i,j,k) 
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

  !   Multiply  by T_zeta^inv T_xi

  do k=2,kl 
     do j=2,jl
        do i=2,il

           rk1  = half*(sk(i,j,k,1) + sk(i,j,k-1,1))
           rk2  = half*(sk(i,j,k,2) + sk(i,j,k-1,2))
           rk3  = half*(sk(i,j,k,3) + sk(i,j,k-1,3))
           rk   = sqrt(rk1*rk1+rk2*rk2+rk3*rk3)

           if(rk.gt.zero) then
              rk1  = rk1/rk
              rk2  = rk2/rk
              rk3  = rk3/rk
           else
              rk1  = nhat
              rk2  =-nhat
              rk3  = nhat
              print*,'rk blowing up'
           endif

           el1k = rk2-rk3
           el2k = rk3-rk1
           el3k = rk1-rk2
           elk  = sqrt(el1k*el1k+el2k*el2k+el3k*el3k)
           if(elk.gt.zero) then
              el1k = el1k/elk
              el2k = el2k/elk
              el3k = el3k/elk
           else
              print*,'elk blowing up'
           endif

           em1k = rk2*el3k - rk3*el2k
           em2k = rk3*el1k - rk1*el3k
           em3k = rk1*el2k - rk2*el1k


           ri1  = half*(si(i,j,k,1) + si(i-1,j,k,1))
           ri2  = half*(si(i,j,k,2) + si(i-1,j,k,2))
           ri3  = half*(si(i,j,k,3) + si(i-1,j,k,3))
           ri   = sqrt(ri1*ri1+ri2*ri2+ri3*ri3)
           if(ri.gt.zero) then
              ri1  = ri1/ri
              ri2  = ri2/ri
              ri3  = ri3/ri
           else
              ri1  = nhat
              ri2  =-nhat
              ri3  = nhat
              print*,'ri blowing up'
           endif

           el1i = ri2-ri3
           el2i = ri3-ri1
           el3i = ri1-ri2
           eli  = sqrt(el1i*el1i+el2i*el2i+el3i*el3i)
           if(eli.gt.zero) then
              el1i = el1i/eli
              el2i = el2i/eli
              el3i = el3i/eli
           else
              print*,'eli blowing up'
           endif

           em1i = ri2*el3i - ri3*el2i
           em2i = ri3*el1i - ri1*el3i
           em3i = ri1*el2i - ri2*el1i

           dw2   = dw(i,j,k,2)
           dw3   = dw(i,j,k,3)
           dw4   = dw(i,j,k,4)
           dw5   = dw(i,j,k,5)

           cijk = sqrt(gamma(i,j,k)*p(i,j,k)/w(i,j,k,irho))
           xfact= half/cijk 
           a1    = dw4+dw5
           a2    = dw4-dw5

           a3    = el1i*el1k+el2i*el2k+el3i*el3k
           a4    = em1i*el1k+em2i*el2k+em3i*el3k
           a5    =  ri1*el1k+ ri2*el2k+ ri3*el3k

           a6    = el1i*em1k+el2i*em2k+el3i*em3k
           a7    = em1i*em1k+em2i*em2k+em3i*em3k
           a8    =  ri1*em1k+ ri2*em2k+ ri3*em3k

           a9   = (el1i*rk1+el2i*rk2+el3i*rk3)*xfact
           a10  = (em1i*rk1+em2i*rk2+em3i*rk3)*xfact
           a11  =   ri1*rk1+ ri2*rk2+ ri3*rk3

           dw(i,j,k,2)= a3*dw2+a4*dw3+a5*cijk*a2
           dw(i,j,k,3)= a6*dw2+a7*dw3+a8*cijk*a2
           dw(i,j,k,4)= a9*dw2+a10*dw3+half*(a1+a11*a2)
           dw(i,j,k,5)=-a9*dw2-a10*dw3+half*(a1-a11*a2)

        enddo
     enddo
  enddo

  !  Multiply by diagonal

  do k=2,kl
     do j=2,jl
        do i=2,il
           dw(i,j,k,1)=dw(i,j,k,1)
           dw(i,j,k,2)=dw(i,j,k,2)
           dw(i,j,k,3)=dw(i,j,k,3)
           dw(i,j,k,4)=dw(i,j,k,4)
           dw(i,j,k,5)=dw(i,j,k,5)
        enddo
     enddo
  enddo

  if(kl.gt.2) then

     !  Inversion in k

     do j=2,jl
        do i=2,il
           do k=2,kl

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
                 bbk(k+1,n)= -diagPlus(n)*vol(i,j,k+1)/vol(i,j,k)
                 ddk(k-1,n)=  diagMinus(n)*vol(i,j,k-1)/vol(i,j,k)
                 !     bbk(k+1,n)= -diagPlus(n)
                 !     ddk(k-1,n)=  diagMinus(n)
                 cck(k  ,n)=  diagPlus(n)-diagMinus(n)
              enddo
           enddo

           do n=1,5
              bbk(ke,n)=zero
              ddk(1 ,n)=zero
              do k=2,kl
                 bbk(k,n)=(bbk(k,n)-viscTerm_K(i,j,k))*dual_dt(i,j,k)
                 ddk(k,n)=(ddk(k,n)-viscTerm_K(i,j,k))*dual_dt(i,j,k)
                 cck(k,n)=one+(cck(k,n)+two*viscTerm_K(i,j,k))*dual_dt(i,j,k)
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

  ! Multiply by T_zeta

  do k=2,kl
     do j=2,jl
        do i=2,il
           uvel  = w(i,j,k,ivx)
           vvel  = w(i,j,k,ivy)
           wvel  = w(i,j,k,ivz)

           rk1  = half*(sk(i,j,k,1) + sk(i,j,k-1,1))
           rk2  = half*(sk(i,j,k,2) + sk(i,j,k-1,2))
           rk3  = half*(sk(i,j,k,3) + sk(i,j,k-1,3))
           rk   = sqrt(rk1*rk1+rk2*rk2+rk3*rk3)
           if(rk.gt.zero) then
              rk1  = rk1/rk
              rk2  = rk2/rk
              rk3  = rk3/rk
           else
              rk1  = nhat
              rk2  =-nhat
              rk3  = nhat
              print*,'rk blowing up'
           endif

           el1  = rk2-rk3
           el2  = rk3-rk1
           el3  = rk1-rk2
           el   = sqrt(el1*el1+el2*el2+el3*el3)
           if(el.gt.zero) then
              el1  = el1/el
              el2  = el2/el
              el3  = el3/el
           else
              print*,'el blowing up'
           endif

           em1  = rk2*el3 - rk3*el2
           em2  = rk3*el1 - rk1*el3
           em3  = rk1*el2 - rk2*el1

           uu   = uvel*rk1+vvel*rk2+wvel*rk3
           vv   = uvel*el1+vvel*el2+wvel*el3
           ww   = uvel*em1+vvel*em2+wvel*em3


           uvw   = half*(uvel*uvel+vvel*vvel+wvel*wvel)
           cijk  = sqrt(gamma(i,j,k)*p(i,j,k)/w(i,j,k,irho))
           gm1   = gamma(i,j,k)-one
           xfact = cijk*cijk/gm1

           dw1   = dw(i,j,k,1)
           dw2   = dw(i,j,k,2)
           dw3   = dw(i,j,k,3)
           dw4   = dw(i,j,k,4)
           dw5   = dw(i,j,k,5)

           a1    =  dw4+dw5
           a2    =  dw4-dw5

           dw(i,j,k,1)  = a1+dw1
           dw(i,j,k,2)  = uvel*dw1+el1*dw2+em1*dw3+rk1*cijk*a2+uvel*a1
           dw(i,j,k,3)  = vvel*dw1+el2*dw2+em2*dw3+rk2*cijk*a2+vvel*a1
           dw(i,j,k,4)  = wvel*dw1+el3*dw2+em3*dw3+rk3*cijk*a2+wvel*a1
           dw(i,j,k,5)  = uvw*dw1 +vv*dw2+ww*dw3+cijk*uu*a2+(uvw+xfact)*a1

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
  real(kind=realType) :: d0,d1,d2

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
