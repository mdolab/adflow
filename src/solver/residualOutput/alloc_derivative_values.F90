! This is a special function that is sued to alloc derivative values
! in blockpointers_d for use with the AD code.

subroutine alloc_derivative_values(nn)
  
  use blockPointers_d ! This modules includes blockPointers
  use inputtimespectral
  use flowvarrefstate
  use inputPhysics
  implicit none

  integer(kind=intType) :: nn,sps,ierr

  ! First create a flowdoms-like structure that is of length
  ! ntimeintervalspectral for the current block

  allocate(flowDomsd(nTimeIntervalsSpectral),stat=ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Call setPointers to get the size info we need

  do sps=1,nTimeIntervalsSpectral
     call setPointers(nn,1,sps)

     allocate(flowDomsd(sps)%x(0:ie,0:je,0:ke,3), stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     allocate(flowDomsd(sps)%si(0:ie,1:je,1:ke,3), &
          flowDomsd(sps)%sj(1:ie,0:je,1:ke,3), &
          flowDomsd(sps)%sk(1:ie,1:je,0:ke,3),stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     allocate(flowDomsd(sps)%rotMatrixI(il,2:jl,2:kl,3,3), &
          flowDomsd(sps)%rotMatrixJ(2:il,jl,2:kl,3,3), &
          flowDomsd(sps)%rotMatrixK(2:il,2:jl,kl,3,3),stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     allocate(flowDomsd(sps)%s(ie,je,ke,3),      &
          flowDomsd(sps)%sFaceI(0:ie,je,ke), &
          flowDomsd(sps)%sFaceJ(ie,0:je,ke), &
          flowDomsd(sps)%sFaceK(ie,je,0:ke), stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     allocate(flowDomsd(sps)%w(0:ib,0:jb,0:kb,1:nw), &
          flowDomsd(sps)%dw(0:ib,0:jb,0:kb,1:nw), &
          flowDomsd(sps)%fw(0:ib,0:jb,0:kb,1:nw), &
          stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     allocate(flowDomsd(sps)%p(0:ib,0:jb,0:kb), &
          flowDomsd(sps)%gamma(0:ib,0:jb,0:kb), stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

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
     
     ! It appears these values only require 1 sps instance
     sps1RansTest: if(sps == 1 .and. &
          equations == RANSEquations) then
        
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
     end if sps1RansTest
     
     if (viscous) then
        allocate(flowDomsd(sps)%d2Wall(2:il,2:jl,2:kl), &
             stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
     end if
  end do
  
end subroutine alloc_derivative_values
