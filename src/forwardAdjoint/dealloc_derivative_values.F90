! This is a special function that is sued to dealloc derivative values
! in blockpointers_d for use with the AD code.

subroutine dealloc_derivative_values(nn)

  use blockPointers_d ! This modules includes blockPointers
  use inputtimespectral
  use flowvarrefstate
  use inputPhysics
  use forwardAdjointVars
  implicit none
  integer(kind=intType) :: nn,sps,ierr

  ! Reset w and dw -> Its like nothing happened...
  deallocatespectral: do sps=1,nTimeIntervalsSpectral
     call setPointersAdj(nn,1,sps)
     ! Reset w 
     flowDoms(nn,1,sps)%w = flowDomsd(sps)%wtmp

     ! Set dw
     flowDoms(nn,1,sps)%dw =flowDomsd(sps)%dwtmp
     
     ! Deallocate memtory
     deallocate(flowDomsd(sps)%dwtmp,stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)
     deallocate(flowDomsd(sps)%xtmp,stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)
     deallocate(flowDomsd(sps)%dwtmp2,stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)
     deallocate(flowDomsd(sps)%dw_deriv,stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)
     deallocate(flowDomsd(sps)%wtmp,stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
  end do deallocatespectral

  do sps=1,nTimeIntervalsSpectral
     call setPointers(nn,1,sps)

     deallocate(flowDomsd(sps)%x, stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     deallocate(flowDomsd(sps)%si, &
                flowDomsd(sps)%sj, &
                flowDomsd(sps)%sk,stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     deallocate(flowDomsd(sps)%rotMatrixI, &
              flowDomsd(sps)%rotMatrixJ, &
              flowDomsd(sps)%rotMatrixK,stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     deallocate(flowDomsd(sps)%s,      &
              flowDomsd(sps)%sFaceI, &
              flowDomsd(sps)%sFaceJ, &
              flowDomsd(sps)%sFaceK, stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     deallocate(flowDomsd(sps)%w, &
              flowDomsd(sps)%dw, &
              flowDomsd(sps)%fw, &
              stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     deallocate(flowDomsd(sps)%p, &
                flowDomsd(sps)%gamma, stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     if( viscous ) then
        deallocate(flowDomsd(sps)%rlv, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
     end if

     if( eddyModel ) then
        deallocate(flowDomsd(sps)%rev,stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
     end if

     deallocate(flowDomsd(sps)%dtl, &
                flowDomsd(sps)%radI,     &
                flowDomsd(sps)%radJ,     &
                flowDomsd(sps)%radK,stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     ! It appears these values only require 1 sps instance
!      sps1RansTest: if(sps == 1 .and. &
!           equations == RANSEquations) then
        if (sps==1) then
        deallocate(flowDomsd(sps)%bmti1,&
                 flowDomsd(sps)%bmti2,&
                 flowDomsd(sps)%bmtj1,&
                 flowDomsd(sps)%bmtj2,&
                 flowDomsd(sps)%bmtk1,&
                 flowDomsd(sps)%bmtk2,&
                 flowDomsd(sps)%bvti1,&
                 flowDomsd(sps)%bvti2,&
                 flowDomsd(sps)%bvtj1,&
                 flowDomsd(sps)%bvtj2,&
                 flowDomsd(sps)%bvtk1,&
                 flowDomsd(sps)%bvtk2,&
                 stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
     end if
     
     if (viscous) then
        deallocate(flowDomsd(sps)%d2Wall, &
             stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
     end if
  end do

  ! Deallocate the color array in flowDoms
 
  deallocate(flowDomsd(1)%color,stat=ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  ! Finally deallocate flowdomsd
  deallocate(flowdomsd,stat=ierr)
  call EChk(ierr,__FILE__,__LINE__)
 
end subroutine dealloc_derivative_values
