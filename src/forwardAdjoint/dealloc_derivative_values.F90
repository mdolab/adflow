! This is a special function that is sued to dealloc derivative values
! in blockpointers_d for use with the AD code.

subroutine dealloc_derivative_values(nn)

  use blockPointers_d ! This modules includes blockPointers
  use inputtimespectral
  use flowvarrefstate
  use inputPhysics
  use cgnsGrid

  implicit none
  integer(kind=intType) :: nn,sps,ierr,i

  ! Reset w and dw -> Its like nothing happened...
  deallocatespectral: do sps=1,nTimeIntervalsSpectral
     call setPointersAdj(nn,1,sps)
     ! Reset w 
     flowDoms(nn,1,sps)%w = flowDomsd(nn,1,sps)%wtmp

     ! Set dw
     flowDoms(nn,1,sps)%dw =flowDomsd(nn,1,sps)%dwtmp
     
     ! Deallocate memtory
     deallocate(&
          flowDomsd(nn,1,sps)%wtmp, &
          flowDomsd(nn,1,sps)%dwtmp, &
          flowDomsd(nn,1,sps)%dwtmp2, &
          flowDomsd(nn,1,sps)%xtmp, &
          flowDomsd(nn,1,sps)%dw_deriv, stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)
  end do deallocatespectral

  do sps=1,nTimeIntervalsSpectral
     call setPointers(nn,1,sps)

     ! Allocate the tempHalo locations in BOTH the normal and AD calcs
     deallocate(flowDoms (nn,1,sps)%tempHalo)
     deallocate(flowDomsd(nn,1,sps)%tempHalo)


     deallocate(flowDomsd(nn,1,sps)%x, stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     deallocate(flowDomsd(nn,1,sps)%si, &
                flowDomsd(nn,1,sps)%sj, &
                flowDomsd(nn,1,sps)%sk, &
                flowDomsd(nn,1,sps)%vol, stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     deallocate(flowDomsd(nn,1,sps)%rotMatrixI, &
              flowDomsd(nn,1,sps)%rotMatrixJ, &
              flowDomsd(nn,1,sps)%rotMatrixK,stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     deallocate(flowDomsd(nn,1,sps)%s,      &
              flowDomsd(nn,1,sps)%sFaceI, &
              flowDomsd(nn,1,sps)%sFaceJ, &
              flowDomsd(nn,1,sps)%sFaceK, stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     deallocate(flowDomsd(nn,1,sps)%w, &
              flowDomsd(nn,1,sps)%dw, &
              flowDomsd(nn,1,sps)%fw, &
              stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     deallocate(flowDomsd(nn,1,sps)%p, &
                flowDomsd(nn,1,sps)%gamma, stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     if( viscous ) then
        deallocate(flowDomsd(nn,1,sps)%rlv, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
     end if

     if( eddyModel ) then
        deallocate(flowDomsd(nn,1,sps)%rev,stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
     end if

     deallocate(flowDomsd(nn,1,sps)%dtl, &
                flowDomsd(nn,1,sps)%radI,     &
                flowDomsd(nn,1,sps)%radJ,     &
                flowDomsd(nn,1,sps)%radK,stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Set the pointer for BCData and deallocate the memory stored 

     BCDatad => flowDomsd(nn,1,sps)%BCData
     do i=1,flowDomsd(nn,1,sps)%nBocos
        
        if( associated(BCDatad(i)%norm) ) &
             deallocate(BCDatad(i)%norm, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        if( associated(BCDatad(i)%rface) ) &
             deallocate(BCDatad(i)%rface, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        if( associated(BCDatad(i)%uSlip) ) &
             deallocate(BCDatad(i)%uSlip, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        if( associated(BCDatad(i)%TNS_Wall) ) &
             deallocate(BCDatad(i)%TNS_Wall, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        if( associated(BCDatad(i)%ptInlet) ) &
             deallocate(BCDatad(i)%ptInlet, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        if( associated(BCDatad(i)%ttInlet) ) &
             deallocate(BCDatad(i)%ttInlet, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        if( associated(BCDatad(i)%htInlet) ) &
             deallocate(BCDatad(i)%htInlet, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        if( associated(BCDatad(i)%flowXdirInlet) ) &
             deallocate(BCDatad(i)%flowXdirInlet, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        if( associated(BCDatad(i)%flowYdirInlet) ) &
             deallocate(BCDatad(i)%flowYdirInlet, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        if( associated(BCDatad(i)%flowZdirInlet) ) &
             deallocate(BCDatad(i)%flowZdirInlet, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        if( associated(BCDatad(i)%rho) ) &
             deallocate(BCDatad(i)%rho, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        if( associated(BCDatad(i)%velx) ) &
             deallocate(BCDatad(i)%velx, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        if( associated(BCDatad(i)%vely) ) &
             deallocate(BCDatad(i)%vely, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        if( associated(BCDatad(i)%velz) ) &
             deallocate(BCDatad(i)%velz, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        if( associated(BCDatad(i)%ps) ) &
             deallocate(BCDatad(i)%ps, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        if( associated(BCDatad(i)%turbInlet) ) &
             deallocate(BCDatad(i)%turbInlet, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
  
        nullify(BCDatad(i)%norm)
        nullify(BCDatad(i)%rface)
        nullify(BCDatad(i)%uSlip)
        nullify(BCDatad(i)%TNS_Wall)
        nullify(BCDatad(i)%ptInlet)
        nullify(BCDatad(i)%ttInlet)
        nullify(BCDatad(i)%htInlet)
        nullify(BCDatad(i)%flowXdirInlet)
        nullify(BCDatad(i)%flowYdirInlet)
        nullify(BCDatad(i)%flowZdirInlet)
        nullify(BCDatad(i)%rho)
        nullify(BCDatad(i)%velx)
        nullify(BCDatad(i)%vely)
        nullify(BCDatad(i)%velz)
        nullify(BCDatad(i)%ps)
        nullify(BCDatad(i)%turbInlet)
        
     enddo
     if( associated(flowDomsd(nn,1,sps)%BCData) ) &
          deallocate(flowDomsd(nn,1,sps)%BCData, stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)
  
 
     if (sps==1) then
        deallocate(flowDomsd(nn,1,sps)%bmti1,&
             flowDomsd(nn,1,sps)%bmti2,&
             flowDomsd(nn,1,sps)%bmtj1,&
             flowDomsd(nn,1,sps)%bmtj2,&
             flowDomsd(nn,1,sps)%bmtk1,&
             flowDomsd(nn,1,sps)%bmtk2,&
             flowDomsd(nn,1,sps)%bvti1,&
             flowDomsd(nn,1,sps)%bvti2,&
             flowDomsd(nn,1,sps)%bvtj1,&
             flowDomsd(nn,1,sps)%bvtj2,&
             flowDomsd(nn,1,sps)%bvtk1,&
             flowDomsd(nn,1,sps)%bvtk2,&
             stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
     end if
     
     if (viscous) then
        deallocate(flowDomsd(nn,1,sps)%d2Wall, &
             stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
     end if
  end do

  ! Dealloc mass info
  deallocate(massFlowFamilyInvd,massFlowFamilyDissd,stat=ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Deallocate the color array in flowDoms
  deallocate(flowDomsd(nn,1,1)%color,stat=ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  ! Finally deallocate flowdomsd
  deallocate(flowdomsd,stat=ierr)
  call EChk(ierr,__FILE__,__LINE__)
 
end subroutine dealloc_derivative_values
