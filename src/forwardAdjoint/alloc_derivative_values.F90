! This is a special function that is sued to alloc derivative values
! in blockpointers_d for use with the AD code.

subroutine alloc_derivative_values(nn, level)

  use blockPointers_d ! This modules includes blockPointers

  use inputtimespectral
  use flowvarrefstate
  use inputPhysics
  use BCTypes
  use cgnsGrid 
  implicit none

  ! Input parameters
  integer(kind=intType) :: nn, level

  ! Local variables
  integer(kind=intType) :: sps,ierr,i,j,k,l, mm
  integer(kind=intType) :: iBeg, jBeg, iEnd, jEnd
  integer(kind=intType) :: massShape(2), max_face_size

  real(kind=realType) :: alpha, beta, Lift, Drag, CL, CD
  real(kind=realType), dimension(3) :: Force, Moment, cForce, cMoment
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
  allocate(winfd(nw),stat=ierr)
  call EChk(ierr,__FILE__,__LINE__)

  do sps=1,nTimeIntervalsSpectral
     call setPointers(nn,level,sps)

     ! Allocate the tempHalo locations in BOTH the normal and AD calcs
     ! Number of nodes including halos = ie+1 = ib, etc
     ! max_face_size = 2*(ib*jb + ib*kb + jb*kb)
     ! allocate(flowDoms (nn,1,sps)%tempHalo(3,max_face_size))
     ! allocate(flowDomsd(nn,1,sps)%tempHalo(3,max_face_size))
     ! flowDoms (nn,1,sps)%tempHalo = zero
     ! flowDomsd(nn,1,sps)%tempHalo = zero
     
     allocate(flowDomsd(nn,1,sps)%x(0:ie,0:je,0:ke,3), stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)
     flowDomsd(nn,1,sps)%x = zero

     allocate(flowDomsd(nn,1,sps)%vol(0:ib,0:jb,0:kb),&
          flowDomsd(nn,1,sps)%si(0:ie,1:je,1:ke,3), &
          flowDomsd(nn,1,sps)%sj(1:ie,0:je,1:ke,3), &
          flowDomsd(nn,1,sps)%sk(1:ie,1:je,0:ke,3), stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)
     flowDomsd(nn,1,sps)%si = zero
     flowDomsd(nn,1,sps)%sj = zero
     flowDomsd(nn,1,sps)%sk = zero
     flowDomsd(nn,1,sps)%vol = zero

     allocate(flowDomsd(nn,1,sps)%rotMatrixI(il,2:jl,2:kl,3,3), &
          flowDomsd(nn,1,sps)%rotMatrixJ(2:il,jl,2:kl,3,3), &
          flowDomsd(nn,1,sps)%rotMatrixK(2:il,2:jl,kl,3,3),stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     allocate(flowDomsd(nn,1,sps)%s(ie,je,ke,3),      &
          flowDomsd(nn,1,sps)%sFaceI(0:ie,je,ke), &
          flowDomsd(nn,1,sps)%sFaceJ(ie,0:je,ke), &
          flowDomsd(nn,1,sps)%sFaceK(ie,je,0:ke), stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     flowDomsd(nn,1,sps)%s = zero
     flowDomsd(nn,1,sps)%sFaceI = zero
     flowDomsd(nn,1,sps)%sFaceJ = zero
     flowDomsd(nn,1,sps)%sFaceK = zero

     allocate(&
          flowDomsd(nn,1,sps)%w (0:ib,0:jb,0:kb,1:nw), &
          flowDomsd(nn,1,sps)%dw(0:ib,0:jb,0:kb,1:nw), &
          flowDomsd(nn,1,sps)%fw(0:ib,0:jb,0:kb,1:nw), stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     flowDomsd(nn,1,sps)%w = zero
     flowDomsd(nn,1,sps)%dw = zero
     flowDomsd(nn,1,sps)%fw = zero

     allocate(flowDomsd(nn,1,sps)%p(0:ib,0:jb,0:kb), &
          flowDomsd(nn,1,sps)%gamma(0:ib,0:jb,0:kb), stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     flowDomsd(nn,1,sps)%p = zero
     flowDomsd(nn,1,sps)%gamma = zero

     ! Nominally we would only need these if visous is True, but
     ! tapende INSISTS on puting rlvd = 0 and rev = 0, it will
     ! segfault if we don't allocate it.
     allocate(flowDomsd(nn,1,sps)%rlv(0:ib,0:jb,0:kb), stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)
     flowDomsd(nn,1,sps)%rlv = zero

     allocate(flowDomsd(nn,1,sps)%rev(0:ib,0:jb,0:kb),stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)
     flowDomsd(nn,1,sps)%rev = zero

     allocate(&
          flowDomsd(nn,1,sps)%dtl (1:ie,1:je,1:ke), &
          flowDomsd(nn,1,sps)%radI(1:ie,1:je,1:ke),     &
          flowDomsd(nn,1,sps)%radJ(1:ie,1:je,1:ke),     &
          flowDomsd(nn,1,sps)%radK(1:ie,1:je,1:ke),stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     flowDomsd(nn,1,sps)%dtl  = zero
     flowDomsd(nn,1,sps)%radI = zero
     flowDomsd(nn,1,sps)%radJ = zero
     flowDomsd(nn,1,sps)%radK = zero

     !allocate memory for boundary condition derivatives
     allocate(flowDomsd(nn,1,sps)%BCData(nBocos), stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     flowdomsd(nn,1,sps)%nBocos = flowdoms(nn,level,sps)%nbocos   
     bcdatad =>flowDomsd(nn,1,sps)%BCData

     bocoLoop: do mm=1,nBocos
        
        ! Store the cell range of the boundary subface
        ! a bit easier.

        iBeg = BCData(mm)%icbeg; iEnd = BCData(mm)%icend
        jBeg = BCData(mm)%jcbeg; jEnd = BCData(mm)%jcend

        allocate(BCDatad(mm)%norm(iBeg:iEnd,jBeg:jEnd,3), stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        bcDatad(mm)%norm = zero

        allocate(BCDatad(mm)%rface(iBeg:iEnd,jBeg:jEnd), stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        bcDatad(mm)%rface = zero

        allocate(BCDatad(mm)%F(iBeg:iEnd,jBeg:jEnd,3), stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        bcDatad(mm)%F = zero

        ! Determine the boundary condition we are having here
        ! and allocate the memory accordingly.

        select case (BCType(mm))

        case (NSWallAdiabatic)

           ! Adiabatic wall. Just allocate the memory for uSlip.

           allocate(BCDatad(mm)%uSlip(iBeg:iEnd,jBeg:jEnd,3), stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
           BCData(mm)%uSlip = zero

           !=======================================================

        case (NSWallIsothermal)

           ! Isothermal wall. Allocate the memory for uSlip
           ! and TNS_Wall.

           allocate(BCDatad(mm)%uSlip(iBeg:iEnd,jBeg:jEnd,3),  &
                BCDatad(mm)%TNS_Wall(iBeg:iEnd,jBeg:jEnd), &
                stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
           BCData(mm)%uSlip = zero
           !=======================================================
        end select
     end do bocoLoop

     ! It appears these values only require 1 sps instance
     !    sps1RansTest: if(sps == 1 .and. &
     !           equations == RANSEquations) then
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

     if (viscous) then
        allocate(flowDomsd(nn,1,sps)%d2Wall(2:il,2:jl,2:kl), &
             stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
     end if
  end do

  ! Also allocate a "color" array for the derivative calcs. Only do
  ! this on flowDomsd, only on the 1st timeInstance
  
  allocate(flowDomsd(nn,1,1)%color(0:ib,0:jb,0:kb),stat=ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  ! Finally Allocate wtmp,dw_deriv and dwtmp in the flowDomsd structure
  ! Allocate Memory and copy out w and dw for reference

  allocspectralLoop: do sps=1,nTimeIntervalsSpectral

     call setPointers(nn,level,sps)
     call block_res(nn, sps, .False., .False., &
          alpha, beta, liftIndex, Force, Moment, Lift, Drag, &
          cForce, cMoment, CL, CD)
     
     allocate(flowDomsd(nn,1,sps)%wtmp(0:ib,0:jb,0:kb,1:nw),stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     allocate(flowDomsd(nn,1,sps)%dwtmp(0:ib,0:jb,0:kb,1:nw),stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     allocate(flowDomsd(nn,1,sps)%dwtmp2(0:ib,0:jb,0:kb,1:nw),stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     allocate(flowDomsd(nn,1,sps)%xtmp(0:ie,0:je,0:ke,3),stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     allocate(flowDomsd(nn,1,sps)%dw_deriv(0:ib,0:jb,0:kb,1:nw,1:nw),stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     flowDomsd(nn,1,sps)%dw_deriv(:,:,:,:,:) = 0.0

     ! Set the values
     flowdomsd(nn,1,sps)%wtmp  = w
     flowdomsd(nn,1,sps)%dwtmp = dw
     flowdomsd(nn,1,sps)%xtmp  = x
   
     !call initRes_block(1,nwf,nn,sps)

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

end subroutine alloc_derivative_values
