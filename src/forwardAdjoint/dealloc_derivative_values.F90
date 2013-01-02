! This is a special function that is sued to dealloc derivative values
! in blockpointers_d for use with the AD code.

subroutine dealloc_derivative_values(nn, level)

  use blockPointers_d ! This modules includes blockPointers
  use inputtimespectral
  use flowvarrefstate
  use inputPhysics
  use cgnsGrid
  use BCTypes
  implicit none

  ! Input Parameters
  integer(kind=intType) :: nn, level

  ! Local variables
  integer(kind=intType) :: sps,ierr,i,mm

  ! Reset w and dw -> Its like nothing happened...
  deallocatespectral: do sps=1,nTimeIntervalsSpectral
     call setPointers(nn,level,sps)
     ! Reset w 
     flowDoms(nn,level,sps)%w = flowDomsd(nn,1,sps)%wtmp

     ! Set dw
     flowDoms(nn,level,sps)%dw =flowDomsd(nn,1,sps)%dwtmp
     
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
     call setPointers(nn,level,sps)

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
     do mm=1,flowDomsd(nn,1,sps)%nBocos
        
        ! Norm is always allocated
        deallocate(BCDatad(mm)%norm,stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        select case (BCType(mm))

        case (NSWallAdiabatic)
           deallocate(BCDatad(mm)%uSlip, stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
        case (NSWallIsothermal)
           deallocate(BCDatad(mm)%uSlip,BCDatad(mm)%TNS_Wall,stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
        case (EulerWall,farField)
           deallocate(BCDatad(mm)%rface, stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
        case (SupersonicInflow, DomainInterfaceAll)
           deallocate(BCDatad(mm)%rho, BCDatad(mm)%velx, BCDatad(mm)%vely, &
                BCDatad(mm)%velz, BCDatad(mm)%ps, stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
           if(nt2 >= nt1) then
              deallocate(BCDatad(mm)%turbInlet, stat=ierr)
              call EChk(ierr,__FILE__,__LINE__)
           endif
        case (SubsonicInflow)
           deallocate(BCDatad(mm)%flowXdirInlet,BCDatad(mm)%flowYdirInlet, &
                BCDatad(mm)%flowZdirInlet, BCDatad(mm)%ptInlet, &
                BCDatad(mm)%ttInlet, BCDatad(mm)%htInlet,    &
                BCDatad(mm)%rho, BCDatad(mm)%velx, BCDatad(mm)%vely, &
                BCDatad(mm)%velz, stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
           if(nt2 >= nt1) then
              deallocate(BCDatad(mm)%turbInlet, stat=ierr)
              call EChk(ierr,__FILE__,__LINE__)
           endif
        case (SubsonicOutflow, MassBleedOutflow, &
             DomainInterfaceP)
           deallocate(BCDatad(mm)%ps,stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
        case (DomainInterfaceRhoUVW)
           deallocate(BCDatad(mm)%rho, BCDatad(mm)%velx, BCDatad(mm)%vely, &
                BCDatad(mm)%velz,stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
           if(nt2 >= nt1) then
              deallocate(BCDatad(mm)%turbInlet, stat=ierr)
              call EChk(ierr,__FILE__,__LINE__)
           endif
        case (DomainInterfaceTotal)
           deallocate(BCDatad(mm)%flowXdirInlet, BCDatad(mm)%flowYdirInlet, &
                BCDatad(mm)%flowZdirInlet, BCDatad(mm)%ptInlet, &
                BCDatad(mm)%ttInlet, BCDatad(mm)%htInlet, stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
           if(nt2 >= nt1) then
              deallocate(BCDatad(mm)%turbInlet, stat=ierr)
              call EChk(ierr,__FILE__,__LINE__)
           endif
        case (domainInterfaceRho)
           deallocate(BCDatad(mm)%rho, stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
        end select
     enddo

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
