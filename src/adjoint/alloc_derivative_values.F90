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
  use adjointVars
  implicit none

  ! Input parameters
  integer(kind=intType) :: level

  ! Local variables
  integer(kind=intType) :: sps,ierr,i,j,k,l, mm, nn
  integer(kind=intType) :: iBeg, jBeg, iStop, jStop, isizemax, jsizemax
  integer(kind=intType) :: massShape(2), max_face_size

  real(kind=realType) :: alpha, beta, force(3, nTimeINtervalsSpectral), moment(3, nTimeIntervalsSpectral), sepSensor, Cavitation
  integer(kind=intType) :: liftIndex
 
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
     do sps=1, nTimeIntervalsSpectral
        call VecCreateMPI(SUMB_COMM_WORLD, 1, PETSC_DETERMINE, xSurfVec(1, sps), ierr)
     end do
  end if

  ! Duplicate the PETSc Xsurf Vec, but only on the first level:
  allocate(xSurfVecd(nTimeIntervalsSpectral))
  do sps=1, nTimeIntervalsSpectral
     call VecDuplicate(xSurfVec(1, sps), xSurfVecd(sps), ierr)
     call EChk(ierr,__FILE__,__LINE__)
  end do

  do nn=1, nDom
     do sps=1, nTimeIntervalsSpectral
        call setPointers(nn, level, sps)

        ! Allocate d2wall if not already done so
        if (.not. associated(flowDoms(nn, 1, sps)%d2wall)) then 
           allocate(flowDoms(nn, 1, sps)%d2wall(2:il, 2:jl, 2:kl))
           call EChk(ierr,__FILE__,__LINE__)
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
           
           allocate(flowDomsd(nn,1,sps)%BCData(mm)%Fp(iBeg:iStop, jBeg:jStop, 3),&
                flowDomsd(nn,1,sps)%BCData(mm)%Fv(iBeg:iStop, jBeg:jStop, 3),&
                flowDomsd(nn,1,sps)%BCData(mm)%area(iBeg:iStop, jBeg:jStop), stat=ierr)
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
     end do
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

#endif 
derivVarsAllocated = .True. 
end subroutine alloc_derivative_values
