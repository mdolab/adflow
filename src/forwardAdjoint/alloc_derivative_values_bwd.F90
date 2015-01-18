! This is a special function that is sued to alloc derivative values
! in blockpointers_d for use with the AD code.

subroutine alloc_derivative_values_bwd(level)

  use blockPointers_b ! This modules includes blockPointers
  use communication
  use inputtimespectral
  use flowvarrefstate
  use inputPhysics
  use BCTypes
  use cgnsGrid 
  use paramTurb
  use turbMod
  use inputADjoint
  use wallDistanceData
  use inputDIscretization
  implicit none

  ! Input parameters
  integer(kind=intType) :: level

  ! Local variables
  integer(kind=intType) :: sps,ierr,i,j,k,l, mm, nState, nn
  integer(kind=intType) :: iBeg, jBeg, iEnd, jEnd
  integer(kind=intType) :: massShape(2), max_face_size

  real(kind=realType) :: alpha, beta,  sepSensor, Cavitation
  real(kind=realType) :: force(3,nTimeIntervalsSpectral), moment(3, nTimeIntervalsSpectral)
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

  allocate(flowDomsb(nDom,1,nTimeIntervalsSpectral),stat=ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! winfd hasn't be allocated so we'll do it here
  allocate(winfb(10),stat=ierr) 
  call EChk(ierr,__FILE__,__LINE__)

  ! If we are not using RANS with wallDistance create a dummy xSurfVec
  ! since one doesn't exist yet
  if (.not. wallDistanceNeeded .or. .not. useApproxWallDistance) then 
     call VecCreateMPI(SUMB_COMM_WORLD, 1, PETSC_DETERMINE, xSurfVec(1), ierr)
  end if

  call VecDuplicate(xSurfVec(1), xSurfVecb, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointers(nn,level,sps)

        ! Allocate d2wall if not already done so
        if (.not. associated(flowDoms(nn, 1, sps)%d2wall)) then 
           allocate(flowDoms(nn, 1, sps)%d2wall(2:il, 2:jl, 2:kl))
           call EChk(ierr,__FILE__,__LINE__)
        end if
        
        allocate(flowDomsb(nn, 1, sps)%d2wall(2:il, 2:jl, 2:kl))

        allocate(flowDomsb(nn,1,sps)%x(0:ie,0:je,0:ke,3), stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        allocate(flowDomsb(nn,1,sps)%vol(0:ib,0:jb,0:kb),&
             flowDomsb(nn,1,sps)%si(0:ie,1:je,1:ke,3), &
             flowDomsb(nn,1,sps)%sj(1:ie,0:je,1:ke,3), &
             flowDomsb(nn,1,sps)%sk(1:ie,1:je,0:ke,3), stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        allocate(&
             flowDomsb(nn,1,sps)%rotMatrixI(il,2:jl,2:kl,3,3), &
             flowDomsb(nn,1,sps)%rotMatrixJ(2:il,jl,2:kl,3,3), &
             flowDomsb(nn,1,sps)%rotMatrixK(2:il,2:jl,kl,3,3),stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        allocate(flowDomsb(nn,1,sps)%s(ie,je,ke,3),      &
             flowDomsb(nn,1,sps)%sFaceI(0:ie,je,ke), &
             flowDomsb(nn,1,sps)%sFaceJ(ie,0:je,ke), &
             flowDomsb(nn,1,sps)%sFaceK(ie,je,0:ke), stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)

        allocate(&
             flowDomsb(nn,1,sps)%w (0:ib,0:jb,0:kb,1:nw), &
             flowDomsb(nn,1,sps)%dw(0:ib,0:jb,0:kb,1:nw), &
             flowDomsb(nn,1,sps)%fw(0:ib,0:jb,0:kb,1:nw), stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        allocate(flowDomsb(nn,1,sps)%p(0:ib,0:jb,0:kb), &
             flowDomsb(nn,1,sps)%gamma(0:ib,0:jb,0:kb), stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)

        allocate(&
        flowDomsb(nn,1,sps)%ux(il,jl,kl), &
        flowDomsb(nn,1,sps)%uy(il,jl,kl), &
        flowDomsb(nn,1,sps)%uz(il,jl,kl), &
        flowDomsb(nn,1,sps)%vx(il,jl,kl), &
        flowDomsb(nn,1,sps)%vy(il,jl,kl), &
        flowDomsb(nn,1,sps)%vz(il,jl,kl), &
        flowDomsb(nn,1,sps)%wx(il,jl,kl), &
        flowDomsb(nn,1,sps)%wy(il,jl,kl), &
        flowDomsb(nn,1,sps)%wz(il,jl,kl), &
        flowDomsb(nn,1,sps)%qx(il,jl,kl), &
        flowDomsb(nn,1,sps)%qy(il,jl,kl), &
        flowDomsb(nn,1,sps)%qz(il,jl,kl), stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)


        ! Nominally we would only need these if visous is True, but
        ! tapende INSISTS on puting rlvd = 0 and rev = 0, it will
        ! segfault if we don't allocate it.
        allocate(flowDomsb(nn,1,sps)%rlv(0:ib,0:jb,0:kb), stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        allocate(flowDomsb(nn,1,sps)%rev(0:ib,0:jb,0:kb),stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        allocate(&
             flowDomsb(nn,1,sps)%dtl (1:ie,1:je,1:ke), &
             flowDomsb(nn,1,sps)%radI(1:ie,1:je,1:ke),     &
             flowDomsb(nn,1,sps)%radJ(1:ie,1:je,1:ke),     &
             flowDomsb(nn,1,sps)%radK(1:ie,1:je,1:ke),stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        !allocate memory for boundary condition derivatives
        allocate(flowDomsb(nn,1,sps)%BCData(nBocos), stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        flowdomsb(nn,1,sps)%nBocos = flowdoms(nn,level,sps)%nbocos   
  
        bocoLoop: do mm=1,nBocos
           
           ! Store the cell range of the boundary subface
           ! a bit easier.
           
           iBeg = BCData(mm)%icbeg; iEnd = BCData(mm)%icend
           jBeg = BCData(mm)%jcbeg; jEnd = BCData(mm)%jcend

           allocate(flowDomsb(nn,1,sps)%BCData(mm)%norm(iBeg:iEnd,jBeg:jEnd,3), stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
           
           allocate(flowDomsb(nn,1,sps)%BCData(mm)%rface(iBeg:iEnd,jBeg:jEnd), stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
           
           allocate(flowDomsb(nn,1,sps)%BCData(mm)%Fp(&
                bcData(mm)%inBeg+1:bcData(mm)%inEnd, &
                bcData(mm)%jnBeg+1:bcData(mm)%jnEnd, 3), stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
 
           allocate(flowDomsb(nn,1,sps)%BCData(mm)%Fv(&
                bcData(mm)%inBeg+1:bcData(mm)%inEnd, &
                bcData(mm)%jnBeg+1:bcData(mm)%jnEnd, 3), stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)

           allocate(flowDomsb(nn,1,sps)%BCData(mm)%M(&
                bcData(mm)%inBeg+1:bcData(mm)%inEnd, &
                bcData(mm)%jnBeg+1:bcData(mm)%jnEnd, 3), stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
           
           allocate(flowDomsb(nn,1,sps)%BCData(mm)%sepSensor(&
                bcData(mm)%inBeg+1:bcData(mm)%inEnd, &
                bcData(mm)%jnBeg+1:bcData(mm)%jnEnd), stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
          
           allocate(flowDomsb(nn,1,sps)%BCData(mm)%Cavitation(&
                bcData(mm)%inBeg+1:bcData(mm)%inEnd, &
                bcData(mm)%jnBeg+1:bcData(mm)%jnEnd), stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)

           allocate(flowDomsb(nn,1,sps)%BCData(mm)%oArea(&
                bcData(mm)%inbeg:bcData(mm)%inEnd, &
                bcData(mm)%jnbeg:bcData(mm)%jnEnd), stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)

           allocate(flowDomsb(nn,1,sps)%BCData(mm)%uSlip(iBeg:iEnd,jBeg:jEnd,3), stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)

           allocate(flowDomsb(nn,1,sps)%BCData(mm)%TNS_Wall(iBeg:iEnd,jBeg:jEnd), stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
        end do bocoLoop
        
        if (sps == 1) then
           allocate(flowDomsb(nn,1,sps)%bmti1(je,ke,nt1:nt2,nt1:nt2), &
                flowDomsb(nn,1,sps)%bmti2(je,ke,nt1:nt2,nt1:nt2), &
                flowDomsb(nn,1,sps)%bmtj1(ie,ke,nt1:nt2,nt1:nt2), &
                flowDomsb(nn,1,sps)%bmtj2(ie,ke,nt1:nt2,nt1:nt2), &
                flowDomsb(nn,1,sps)%bmtk1(ie,je,nt1:nt2,nt1:nt2), &
                flowDomsb(nn,1,sps)%bmtk2(ie,je,nt1:nt2,nt1:nt2), &
                flowDomsb(nn,1,sps)%bvti1(je,ke,nt1:nt2), &
                flowDomsb(nn,1,sps)%bvti2(je,ke,nt1:nt2), &
                flowDomsb(nn,1,sps)%bvtj1(ie,ke,nt1:nt2), &
                flowDomsb(nn,1,sps)%bvtj2(ie,ke,nt1:nt2), &
                flowDomsb(nn,1,sps)%bvtk1(ie,je,nt1:nt2), &
                flowDomsb(nn,1,sps)%bvtk2(ie,je,nt1:nt2), &
                stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
       
           flowDomsb(nn,1,sps)%bmti1 = zero
           flowDomsb(nn,1,sps)%bmti2 = zero
           flowDomsb(nn,1,sps)%bmtj1 = zero
           flowDomsb(nn,1,sps)%bmtj2 = zero
           flowDomsb(nn,1,sps)%bmtk1 = zero
           flowDomsb(nn,1,sps)%bmtk2 = zero
           flowDomsb(nn,1,sps)%bvti1 = zero
           flowDomsb(nn,1,sps)%bvti2 = zero
           flowDomsb(nn,1,sps)%bvtj1 = zero
           flowDomsb(nn,1,sps)%bvtj2 = zero
           flowDomsb(nn,1,sps)%bvtk1 = zero
           flowDomsb(nn,1,sps)%bvtk2 = zero
           
        end if
        
        allocate(flowDomsb(nn,1,sps)%d2Wall(2:il,2:jl,2:kl), &
             stat=ierr)
        flowDomsb(nn,1,sps)%d2wall = zero
        call EChk(ierr,__FILE__,__LINE__)
        
        allocate(flowDomsb(nn,1,sps)%viscSubface(nviscBocos), &
             stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
           
        viscbocoLoop: do mm=1,nviscBocos
           
           iBeg = BCData(mm)%inBeg + 1
           iEnd = BCData(mm)%inEnd
           
           jBeg = BCData(mm)%jnBeg + 1
           jEnd = BCData(mm)%jnEnd
            
           allocate(flowDomsb(nn,1,sps)%viscSubface(mm)%tau(iBeg:iEnd,jBeg:jEnd,6), &
                stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
           
           allocate(flowDomsb(nn,1,sps)%viscSubface(mm)%q(iBeg:iEnd,jBeg:jEnd,6), &
                stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
           
        enddo viscbocoLoop
        
        ! Zero out all the derivative values we've just allocated
        call zeroADSeedsBwd(nn, 1, sps)
     end do
  end do
  
  ! Also allocate a "color" array for the derivative calcs. Only do
  ! this on flowDomsb, only on the 1st timeInstance
  do nn=1,nDom
     call setPointers(nn, level, 1)
     allocate(flowDomsb(nn,1,1)%color(0:ib,0:jb,0:kb),stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)
  end do
  ! Finally Allocate wtmp,dw_deriv and dwtmp in the flowDomsb structure
  ! Allocate Memory and copy out w and dw for reference

  do nn=1,nDom
     allocspectralLoop: do sps=1,nTimeIntervalsSpectral

        call setPointers(nn,level,sps)
        
        call block_res(nn, sps, .False., alpha, beta, liftIndex, force, moment, sepSensor, Cavitation)
     
        allocate(flowDomsb(nn,1,sps)%wtmp(0:ib,0:jb,0:kb,1:nw),stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)

        allocate(flowDomsb(nn,1,sps)%dwtmp(0:ib,0:jb,0:kb,1:nw),stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        allocate(flowDomsb(nn,1,sps)%dwtmp2(0:ib,0:jb,0:kb,1:nw),stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        allocate(flowDomsb(nn,1,sps)%xtmp(0:ie,0:je,0:ke,3),stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        allocate(flowDomsb(nn,1,sps)%dw_deriv(0:ib,0:jb,0:kb,1:nState,1:nState),stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        flowDomsb(nn,1,sps)%dw_deriv(:,:,:,:,:) = 0.0

        allocate(flowDomsb(nn,1,sps)%w_deriv(0:ib,0:jb,0:kb,1:nState,1:nState),stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        flowDomsb(nn,1,sps)%w_deriv(:,:,:,:,:) = 0.0

        ! Set the values
        do l=1,nw
           do k=0,kb 
              do j=0,jb
                 do i=0,ib
                    flowDomsb(nn,1,sps)%wtmp(i,j,k,l)  = w(i,j,k,l)
                    flowDomsb(nn,1,sps)%dwtmp(i,j,k,l) = dw(i,j,k,l)
                 end do
              end do
           end do
        end do

        do l=1,3
           do k=0,ke
              do j=0,je
                 do i=0,ie
                    flowDomsb(nn,1,sps)%xtmp(i,j,k,l)  = x(i,j,k,l)
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
                    flowDomsb(nn,1,sps)%dwtmp2(i,j,k,l) = dw(i,j,k,l)/vol(i,j,k)
                 end do
              end do
           end do
        end do
     end do allocspectralLoop
  end do
end subroutine alloc_derivative_values_bwd
