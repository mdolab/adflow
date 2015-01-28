! This is a special function that is sued to alloc derivative values
! in blockpointers_d for use with the AD code.

subroutine alloc_derivative_values(level)

  use blockPointers
  use inputtimespectral
  use flowvarrefstate
  use inputPhysics
  use BCTypes
  use cgnsGrid 
  use paramTurb
  use turbMod
  use inputADjoint
  use inputDiscretization
  use communication
  use wallDistanceData
#ifndef USE_COMPLEX
  use bcroutines_b
#endif
  implicit none

  ! Input parameters
  integer(kind=intType) :: level

  ! Local variables
  integer(kind=intType) :: sps,ierr,i,j,k,l, mm, nState, nn
  integer(kind=intType) :: iBeg, jBeg, iStop, jStop, isizemax, jsizemax
  integer(kind=intType) :: massShape(2), max_face_size

  real(kind=realType) :: alpha, beta, force(3, nTimeINtervalsSpectral), moment(3, nTimeIntervalsSpectral), sepSensor, Cavitation
  integer(kind=intType) :: liftIndex
  
  ! Setup number of state variable based on turbulence assumption
  if ( frozenTurbulence ) then
     nState = nwf
  else
     nState = nw
  endif
  
  ! This routine will not use the extra variables to block_res or the
  ! extra outputs, so we must zero them here

  call getDirAngle(velDirFreestream, liftDirection, liftIndex, alpha, beta)

  ! First create a flowdoms-like structure that is nominally the same
  ! shape as flowDoms. However we will only ALLOCATE values for block
  ! nn

  allocate(flowDomsd(nDom,1,nTimeIntervalsSpectral),stat=ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! winfd hasn't be allocated so we'll do it here
  allocate(winfd(10),stat=ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! If we are not using RANS with walDistance create a dummy xSurfVec
  ! since one does not yet exist
  if (.not. wallDistanceNeeded .or. .not. useApproxWallDistance) then 
     call VecCreateMPI(SUMB_COMM_WORLD, 1, PETSC_DETERMINE, xSurfVec(1), ierr)
  end if

  ! Duplicate the PETSc Xsurf Vec, but only on the first level:
  call VecDuplicate(xSurfVec(1), xSurfVecd, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointers(nn,level,sps)

        ! Allocate d2wall if not already done so
        if (.not. associated(flowDoms(nn, 1, sps)%d2wall)) then 
           allocate(flowDoms(nn, 1, sps)%d2wall(2:il, 2:jl, 2:kl))
           call EChk(ierr,__FILE__,__LINE__)
        end if

        ! Allocate shockSensor in flowDoms *NOT* flowDomsd....and
        ! compute the value depending on equations/dissipation
        ! type. Note we are just doing all cells including corners
        ! halos..those values are not used anyway. 
        allocate(flowDoms(nn,1,sps)%shockSensor(0:ib,0:jb,0:kb))
        shockSensor => flowDoms(nn,1,sps)%shockSensor
        if (equations == EulerEquations .or. spaceDiscr == dissMatrix) then 
           !shockSensor is Pressure
           do k=0,kb
              do j=0,jb
                 do i=0,ib
                    shockSensor(i,j,k) = P(i,j,k)
                 end do
              end do
           end do
        else
           ! Enthalpy is used instead
           do k=0,kb
              do j=2,jl
                 do i=2,il
                    shockSensor(i,j,k) = p(i,j,k)/(w(i,j,k,irho)**gamma(i,j,k))
                 enddo
              enddo
           enddo
           
           do k=2,kl
              do j=2,jl
                 shockSensor(0, j,k) = p(0, j,k)/(w(0, j,k,irho)**gamma(0, j,k))
                 shockSensor(1, j,k) = p(1, j,k)/(w(1, j,k,irho)**gamma(1, j,k))
                 shockSensor(ie,j,k) = p(ie,j,k)/(w(ie,j,k,irho)**gamma(ie,j,k))
                 shockSensor(ib,j,k) = p(ib,j,k)/(w(ib,j,k,irho)**gamma(ib,j,k))
              enddo
           enddo
           
           do k=2,kl
              do i=2,il
                 shockSensor(i,0, k) = p(i,0, k)/(w(i,0, k,irho)**gamma(i,0, k))
                 shockSensor(i,1, k) = p(i,1, k)/(w(i,1, k,irho)**gamma(i,1, k))
                 shockSensor(i,je,k) = p(i,je,k)/(w(i,je,k,irho)**gamma(i,je,k))
                 shockSensor(i,jb,k) = p(i,jb,k)/(w(i,jb,k,irho)**gamma(i,jb,k))
              enddo
           enddo
        end if

        allocate(flowDomsd(nn,1,sps)%x(0:ie,0:je,0:ke,3), stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        allocate(flowDomsd(nn,1,sps)%vol(0:ib,0:jb,0:kb),&
             flowDomsd(nn,1,sps)%si(0:ie,1:je,1:ke,3), &
             flowDomsd(nn,1,sps)%sj(1:ie,0:je,1:ke,3), &
             flowDomsd(nn,1,sps)%sk(1:ie,1:je,0:ke,3), stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        allocate(&
             flowDomsd(nn,1,sps)%rotMatrixI(il,2:jl,2:kl,3,3), &
             flowDomsd(nn,1,sps)%rotMatrixJ(2:il,jl,2:kl,3,3), &
             flowDomsd(nn,1,sps)%rotMatrixK(2:il,2:jl,kl,3,3),stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        allocate(flowDomsd(nn,1,sps)%s(ie,je,ke,3),      &
             flowDomsd(nn,1,sps)%sFaceI(0:ie,je,ke), &
             flowDomsd(nn,1,sps)%sFaceJ(ie,0:je,ke), &
             flowDomsd(nn,1,sps)%sFaceK(ie,je,0:ke), stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)

        allocate(&
             flowDomsd(nn,1,sps)%w (0:ib,0:jb,0:kb,1:nw), &
             flowDomsd(nn,1,sps)%dw(0:ib,0:jb,0:kb,1:nw), &
             flowDomsd(nn,1,sps)%fw(0:ib,0:jb,0:kb,1:nw), &
             flowDomsd(nn,1,sps)%scratch(0:ib,0:jb,0:kb,5), &
             stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        allocate(flowDomsd(nn,1,sps)%p(0:ib,0:jb,0:kb), &
             flowDomsd(nn,1,sps)%gamma(0:ib,0:jb,0:kb),  &
             flowDomsd(nn,1,sps)%aa(0:ib,0:jb,0:kb), stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)

        allocate(&
        flowDomsd(nn,1,sps)%ux(il,jl,kl), &
        flowDomsd(nn,1,sps)%uy(il,jl,kl), &
        flowDomsd(nn,1,sps)%uz(il,jl,kl), &
        flowDomsd(nn,1,sps)%vx(il,jl,kl), &
        flowDomsd(nn,1,sps)%vy(il,jl,kl), &
        flowDomsd(nn,1,sps)%vz(il,jl,kl), &
        flowDomsd(nn,1,sps)%wx(il,jl,kl), &
        flowDomsd(nn,1,sps)%wy(il,jl,kl), &
        flowDomsd(nn,1,sps)%wz(il,jl,kl), &
        flowDomsd(nn,1,sps)%qx(il,jl,kl), &
        flowDomsd(nn,1,sps)%qy(il,jl,kl), &
        flowDomsd(nn,1,sps)%qz(il,jl,kl), stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)

        ! Nominally we would only need these if visous is True, but
        ! tapende INSISTS on puting rlvd = 0 and rev = 0, it will
        ! segfault if we don't allocate it.
        allocate(flowDomsd(nn,1,sps)%rlv(0:ib,0:jb,0:kb), stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        allocate(flowDomsd(nn,1,sps)%rev(0:ib,0:jb,0:kb),stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        allocate(&
             flowDomsd(nn,1,sps)%dtl (1:ie,1:je,1:ke), &
             flowDomsd(nn,1,sps)%radI(1:ie,1:je,1:ke),     &
             flowDomsd(nn,1,sps)%radJ(1:ie,1:je,1:ke),     &
             flowDomsd(nn,1,sps)%radK(1:ie,1:je,1:ke),stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        !allocate memory for boundary condition derivatives
        allocate(flowDomsd(nn,1,sps)%BCData(nBocos), stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        flowdomsd(nn,1,sps)%nBocos = flowdoms(nn,level,sps)%nbocos   
  
        bocoLoop: do mm=1,nBocos
           
           ! Store the cell range of the boundary subface
           ! a bit easier.
           
           iBeg = BCData(mm)%icbeg; iStop = BCData(mm)%icend
           jBeg = BCData(mm)%jcbeg; jStop = BCData(mm)%jcend

           allocate(flowDomsd(nn,1,sps)%BCData(mm)%norm(iBeg:iStop,jBeg:jStop,3), stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
           
           allocate(flowDomsd(nn,1,sps)%BCData(mm)%rface(iBeg:iStop,jBeg:jStop), stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
           
           allocate(flowDomsd(nn,1,sps)%BCData(mm)%Fp(&
                bcData(mm)%inBeg+1:bcData(mm)%inEnd, &
                bcData(mm)%jnBeg+1:bcData(mm)%jnEnd, 3), stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
 
           allocate(flowDomsd(nn,1,sps)%BCData(mm)%Fv(&
                bcData(mm)%inBeg+1:bcData(mm)%inEnd, &
                bcData(mm)%jnBeg+1:bcData(mm)%jnEnd, 3), stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)

           allocate(flowDomsd(nn,1,sps)%BCData(mm)%M(&
                bcData(mm)%inBeg+1:bcData(mm)%inEnd, &
                bcData(mm)%jnBeg+1:bcData(mm)%jnEnd, 3), stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)

           allocate(flowDomsd(nn,1,sps)%BCData(mm)%sepSensor(&
                bcData(mm)%inBeg+1:bcData(mm)%inEnd, &
                bcData(mm)%jnBeg+1:bcData(mm)%jnEnd), stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)

           allocate(flowDomsd(nn,1,sps)%BCData(mm)%Cavitation(&
                bcData(mm)%inBeg+1:bcData(mm)%inEnd, &
                bcData(mm)%jnBeg+1:bcData(mm)%jnEnd), stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
           
           allocate(flowDomsd(nn,1,sps)%BCData(mm)%oArea(&
                bcData(mm)%inbeg:bcData(mm)%inEnd, &
                bcData(mm)%jnbeg:bcData(mm)%jnEnd), stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)

           allocate(flowDomsd(nn,1,sps)%BCData(mm)%uSlip(iBeg:iStop,jBeg:jStop,3), stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)

           allocate(flowDomsd(nn,1,sps)%BCData(mm)%TNS_Wall(iBeg:iStop,jBeg:jStop), stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
        end do bocoLoop
        
        if (sps == 1) then
           allocate(flowDomsd(nn,1,sps)%bmti1(je,ke,nt1:nt2,nt1:nt2), &
                flowDomsd(nn,1,sps)%bmti2(je,ke,nt1:nt2,nt1:nt2), &
                flowDomsd(nn,1,sps)%bmtj1(ie,ke,nt1:nt2,nt1:nt2), &
                flowDomsd(nn,1,sps)%bmtj2(ie,ke,nt1:nt2,nt1:nt2), &
                flowDomsd(nn,1,sps)%bmtk1(ie,je,nt1:nt2,nt1:nt2), &
                flowDomsd(nn,1,sps)%bmtk2(ie,je,nt1:nt2,nt1:nt2), &
                flowDomsd(nn,1,sps)%bvti1(je,ke,nt1:nt2), &
                flowDomsd(nn,1,sps)%bvti2(je,ke,nt1:nt2), &
                flowDomsd(nn,1,sps)%bvtj1(ie,ke,nt1:nt2), &
                flowDomsd(nn,1,sps)%bvtj2(ie,ke,nt1:nt2), &
                flowDomsd(nn,1,sps)%bvtk1(ie,je,nt1:nt2), &
                flowDomsd(nn,1,sps)%bvtk2(ie,je,nt1:nt2), &
                stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
        end if

        allocate(flowDomsd(nn,1,sps)%d2Wall(2:il,2:jl,2:kl), &
             stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        allocate(flowDomsd(nn,1,sps)%viscSubface(nviscBocos), &
             stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        viscbocoLoop: do mm=1,nviscBocos
           
           iBeg = BCData(mm)%inBeg + 1
           iStop = BCData(mm)%inEnd
           
           jBeg = BCData(mm)%jnBeg + 1
           jStop = BCData(mm)%jnEnd
           
           allocate(flowDomsd(nn,1,sps)%viscSubface(mm)%tau(iBeg:iStop,jBeg:jStop,6), &
                stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
           
           allocate(flowDomsd(nn,1,sps)%viscSubface(mm)%q(iBeg:iStop,jBeg:jStop,6), &
                stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
           
        enddo viscbocoLoop

        ! Zero out all the derivative values we've just allocated
        call zeroADSeeds(nn, 1, sps)
     end do
  end do
  
  ! Also allocate a "color" array for the derivative calcs. Only do
  ! this on flowDomsd, only on the 1st timeInstance. This goes from
  ! 0:{i,j,k}b which is necessary for the double halo cells. The same
  ! array is used for the nodal colors, however, in that case only the
  ! 0:{i,j,k}e entries are needed. 
  do nn=1,nDom
     call setPointers(nn, level, 1)
     allocate(flowDomsd(nn,1,1)%color(0:ib,0:jb,0:kb),stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)
  end do
  ! Finally Allocate wtmp,dw_deriv and dwtmp in the flowDomsd structure
  ! Allocate Memory and copy out w and dw for reference

  do nn=1,nDom
     allocspectralLoop: do sps=1,nTimeIntervalsSpectral

        call setPointers(nn,level,sps)
        shockSensor => flowDoms(nn,level,sps)%shockSensor

        call block_res(nn, sps, .False., alpha, beta, liftIndex, force, moment, sepSensor, Cavitation)

        allocate(flowDomsd(nn,1,sps)%wtmp(0:ib,0:jb,0:kb,1:nw),stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)

        allocate(flowDomsd(nn,1,sps)%dwtmp(0:ib,0:jb,0:kb,1:nw),stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        allocate(flowDomsd(nn,1,sps)%dwtmp2(0:ib,0:jb,0:kb,1:nw),stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        allocate(flowDomsd(nn,1,sps)%xtmp(0:ie,0:je,0:ke,3),stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        allocate(flowDomsd(nn,1,sps)%dw_deriv(0:ib,0:jb,0:kb,1:nState,1:nState),stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        flowDomsd(nn,1,sps)%dw_deriv(:,:,:,:,:) = 0.0

        ! Set the values
        do l=1,nw
           do k=0,kb 
              do j=0,jb
                 do i=0,ib
                    flowdomsd(nn,1,sps)%wtmp(i,j,k,l)  = w(i,j,k,l)
                    flowdomsd(nn,1,sps)%dwtmp(i,j,k,l) = dw(i,j,k,l)
                 end do
              end do
           end do
        end do

        do l=1,3
           do k=0,ke
              do j=0,je
                 do i=0,ie
                    flowdomsd(nn,1,sps)%xtmp(i,j,k,l)  = x(i,j,k,l)
                 end do
              end do
           end do
        end do
             
        call initRes_block(1, nwf, nn, sps)
        
        ! Note: we have to divide by the volume for dwtmp2 since
        ! normally, dw would have been mulitpiled by 1/Vol in block_res 
        
        do l=1,nw
           do k=0,kb 
              do j=0,jb
                 do i=0,ib
                    flowdomsd(nn,1,sps)%dwtmp2(i,j,k,l) = dw(i,j,k,l)/vol(i,j,k)
                 end do
              end do
           end do
        end do
     end do allocspectralLoop
  end do

#ifndef USE_COMPLEX
  ! Finally allocate space for the BC pointers
  isizemax = 0
  jsizemax = 0
  do nn=1,nDom
     isizemax = max(isizemax, flowDoms(nn, 1, 1)%ie)
     isizemax = max(isizemax, flowDoms(nn, 1, 1)%je)

     jsizemax = max(jsizemax, flowDoms(nn, 1, 1)%je)
     jsizemax = max(jsizemax, flowDoms(nn, 1, 1)%ke)
  end do

  allocate(&
       ww0d(isizemax, jsizemax, nw), ww1d(isizemax, jsizemax, nw), &
       ww2d(isizemax, jsizemax, nw), ww3d(isizemax, jsizemax, nw), &
       pp0d(isizemax, jsizemax), pp1d(isizemax, jsizemax), &
       pp2d(isizemax, jsizemax), pp3d(isizemax, jsizemax), &
       rlv0d(isizemax, jsizemax), rlv1d(isizemax, jsizemax), &
       rlv2d(isizemax, jsizemax), rlv3d(isizemax, jsizemax), &
       rev0d(isizemax, jsizemax), rev1d(isizemax, jsizemax), &
       rev2d(isizemax, jsizemax), rev3d(isizemax, jsizemax), &
       ssid(isizemax, jsizemax,3), xxd(isizemax+1, jsizemax+1,3), stat=ierr)
  call EChk(ierr,__FILE__,__LINE__)
  allocate(&
       ww0(isizemax, jsizemax, nw), ww1(isizemax, jsizemax, nw), &
       ww2(isizemax, jsizemax, nw), ww3(isizemax, jsizemax, nw), &
       pp0(isizemax, jsizemax), pp1(isizemax, jsizemax), &
       pp2(isizemax, jsizemax), pp3(isizemax, jsizemax), &
       rlv0(isizemax, jsizemax), rlv1(isizemax, jsizemax), &
       rlv2(isizemax, jsizemax), rlv3(isizemax, jsizemax), &
       rev0(isizemax, jsizemax), rev1(isizemax, jsizemax), &
       rev2(isizemax, jsizemax), rev3(isizemax, jsizemax), &
       gamma0(isizemax, jsizemax), gamma1(isizemax, jsizemax), &
       gamma2(isizemax, jsizemax), gamma3(isizemax, jsizemax), &
       ssi(isizemax, jsizemax,3), xx(isizemax+1, jsizemax+1,3), stat=ierr)
  call EChk(ierr,__FILE__,__LINE__) 

  ! Now zero these
  ww0d = zero
  ww1d = zero
  ww2d = zero
  ww3d = zero

  pp0d = zero
  pp1d = zero
  pp2d = zero
  pp3d = zero

  rlv0d = zero
  rlv1d = zero
  rlv2d = zero
  rlv3d = zero

  rev0d = zero
  rev1d = zero
  rev2d = zero
  rev3d = zero
  ssid = zero
  xxd = zero
#endif 
end subroutine alloc_derivative_values
