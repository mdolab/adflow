! This is a special function that is sued to alloc derivative values
! in blockpointers_d for use with the AD code.

subroutine alloc_derivative_values(nn)
  
  use blockPointers_d ! This modules includes blockPointers

  use inputtimespectral
  use flowvarrefstate
  use inputPhysics
  use forwardAdjointVars

  implicit none

  integer(kind=intType) :: nn,sps,ierr,i,j,k,l

  ! First create a flowdoms-like structure that is of length
  ! ntimeintervalspectral for the current block

  allocate(flowDomsd(nTimeIntervalsSpectral),stat=ierr)
  call EChk(ierr,__FILE__,__LINE__)

  do sps=1,nTimeIntervalsSpectral
     call setPointers(nn,1,sps)

     allocate(flowDomsd(sps)%x(0:ie,0:je,0:ke,3), stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)
     flowDomsd(sps)%x = 0

     allocate(flowDomsd(sps)%si(0:ie,1:je,1:ke,3), &
          flowDomsd(sps)%sj(1:ie,0:je,1:ke,3), &
          flowDomsd(sps)%sk(1:ie,1:je,0:ke,3),stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)
     flowDomsd(sps)%si = 0
     flowDomsd(sps)%sj = 0
     flowDomsd(sps)%sk = 0

     allocate(flowDomsd(sps)%rotMatrixI(il,2:jl,2:kl,3,3), &
          flowDomsd(sps)%rotMatrixJ(2:il,jl,2:kl,3,3), &
          flowDomsd(sps)%rotMatrixK(2:il,2:jl,kl,3,3),stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     allocate(flowDomsd(sps)%s(ie,je,ke,3),      &
          flowDomsd(sps)%sFaceI(0:ie,je,ke), &
          flowDomsd(sps)%sFaceJ(ie,0:je,ke), &
          flowDomsd(sps)%sFaceK(ie,je,0:ke), stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     flowDomsd(sps)%s = 0
     flowDomsd(sps)%sFaceI = 0
     flowDomsd(sps)%sFaceJ = 0
     flowDomsd(sps)%sFaceK(ie,je,0:ke) = 0
     
     allocate(flowDomsd(sps)%w(0:ib,0:jb,0:kb,1:nw), &
          flowDomsd(sps)%dw(0:ib,0:jb,0:kb,1:nw), &
          flowDomsd(sps)%fw(0:ib,0:jb,0:kb,1:nw), &
          stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     flowDomsd(sps)%w = 0.0
     flowDomsd(sps)%dw = 0.0
     flowDomsd(sps)%fw = 0.0

     allocate(flowDomsd(sps)%p(0:ib,0:jb,0:kb), &
          flowDomsd(sps)%gamma(0:ib,0:jb,0:kb), stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     flowDomsd(sps)%p = 0.0
     flowDomsd(sps)%gamma = 0.0

     if( viscous ) then
        allocate(flowDomsd(sps)%rlv(0:ib,0:jb,0:kb), stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
     end if

     if( eddyModel ) then
        allocate(flowDomsd(sps)%rev(0:ib,0:jb,0:kb),stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
     end if

     allocate(flowDomsd(sps)%dtl(1:ie,1:je,1:ke), &
              flowDomsd(sps)%radI(1:ie,1:je,1:ke),     &
              flowDomsd(sps)%radJ(1:ie,1:je,1:ke),     &
              flowDomsd(sps)%radK(1:ie,1:je,1:ke),stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     flowDomsd(sps)%dtl = 0.0
     flowDomsd(sps)%radI = 0.0
     flowDomsd(sps)%radJ = 0.0
     flowDomsd(sps)%radK = 0.0

     ! It appears these values only require 1 sps instance
  !    sps1RansTest: if(sps == 1 .and. &
!           equations == RANSEquations) then
        if (sps == 1) then
        allocate(flowDomsd(sps)%bmti1(je,ke,nt1:nt2,nt1:nt2), &
                 flowDomsd(sps)%bmti2(je,ke,nt1:nt2,nt1:nt2), &
                 flowDomsd(sps)%bmtj1(ie,ke,nt1:nt2,nt1:nt2), &
                 flowDomsd(sps)%bmtj2(ie,ke,nt1:nt2,nt1:nt2), &
                 flowDomsd(sps)%bmtk1(ie,je,nt1:nt2,nt1:nt2), &
                 flowDomsd(sps)%bmtk2(ie,je,nt1:nt2,nt1:nt2), &
                 flowDomsd(sps)%bvti1(je,ke,nt1:nt2), &
                 flowDomsd(sps)%bvti2(je,ke,nt1:nt2), &
                 flowDomsd(sps)%bvtj1(ie,ke,nt1:nt2), &
                 flowDomsd(sps)%bvtj2(ie,ke,nt1:nt2), &
                 flowDomsd(sps)%bvtk1(ie,je,nt1:nt2), &
                 flowDomsd(sps)%bvtk2(ie,je,nt1:nt2), &
                 stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
     end if
     
     if (viscous) then
        allocate(flowDomsd(sps)%d2Wall(2:il,2:jl,2:kl), &
             stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
     end if
  end do
 
  ! Also allocate a "color" array for the derivative calcs. Only do
  ! this on flowDomsd, only on the 1st timeInstance
  
  allocate(flowDomsd(1)%color(0:ib,0:jb,0:kb),stat=ierr)
  call EChk(ierr,__FILE__,__LINE__)
 
  ! Finally Allocate wtmp,dw_deriv and dwtmp in the flowDomsd structure
  ! Allocate Memory and copy out w and dw for reference

  allocspectralLoop: do sps=1,nTimeIntervalsSpectral

     call setPointersAdj(nn,1,sps)
     call block_res(nn,sps)
     allocate(flowDomsd(sps)%wtmp(0:ib,0:jb,0:kb,1:nw),stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     allocate(flowDomsd(sps)%dw_deriv(0:ib,0:jb,0:kb,1:nw,1:nw),stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     allocate(flowDomsd(sps)%dwtmp(0:ib,0:jb,0:kb,1:nw),stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     allocate(flowDomsd(sps)%dwtmp2(0:ib,0:jb,0:kb,1:nw),stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)


     flowDomsd(sps)%dw_deriv(:,:,:,:,:) = 0.0
     ! Set the values
     do l=1,nw
        do k=0,kb 
           do j=0,jb
              do i=0,ib
                 flowdomsd(sps)%wtmp(i,j,k,l) = w(i,j,k,l)
              end do
           end do
        end do
     end do

     do l=1,nw
        do k=2,kl
           do j=2,jl
              do i=2,il
                 flowdomsd(sps)%dwtmp(i,j,k,l) = dw(i,j,k,l)
              end do
           end do
        end do
     end do
  
     call initRes_block(1,nwf,nn,sps)
   
     do l=1,nw
        do k=0,kb 
           do j=0,jb
              do i=0,ib
                 flowdomsd(sps)%dwtmp2(i,j,k,l) = dw(i,j,k,l)/vol(i,j,k)
              end do
           end do
        end do
     end do

  end do allocspectralLoop
  
end subroutine alloc_derivative_values
