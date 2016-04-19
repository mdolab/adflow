! This is a special function that is sued to dealloc derivative values
! in blockpointers_d for use with the AD code.

subroutine dealloc_derivative_values(level)

  use blockPointers
  use inputtimespectral
  use flowvarrefstate
  use inputPhysics
  use cgnsGrid
  use BCTypes
  use communication
  use wallDistanceData
#ifndef USE_COMPLEX
  use bcroutines_b
#endif
  use adjointVars
  implicit none

  ! Input Parameters
  integer(kind=intType) :: level

  ! Local variables
  integer(kind=intType) :: sps, ierr, i, j, k, l, mm, nn

  do nn=1,nDom
     ! Deallocate the color array in flowDoms
     deallocate(flowDomsd(nn,1,1)%color,stat=ierr)
     call EChk(ierr,__FILE__,__LINE__)
  end do

  do nn=1,nDom
     ! Reset w and dw -> Its like nothing happened...
     deallocatespectral: do sps=1,nTimeIntervalsSpectral
        call setPointers(nn,level,sps)

        ! Deallocate memtory
        deallocate(&
             flowDomsd(nn,1,sps)%wtmp, &
             flowDomsd(nn,1,sps)%dwtmp, &
             flowDomsd(nn,1,sps)%dwtmp2, &
             flowDomsd(nn,1,sps)%dw_deriv, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
     end do deallocatespectral

     do sps=1,nTimeIntervalsSpectral
        call setPointers(nn,level,sps)

        ! deallocate shockSensor in flowDoms *NOT* flowDomsd
        deallocate(flowDoms(nn,1,sps)%shockSensor, stat=ierr)

        call EChk(ierr,__FILE__,__LINE__)

        deallocate(flowDomsd(nn,1,sps)%x, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        deallocate(flowDomsd(nn,1,sps)%si, &
             flowDomsd(nn,1,sps)%sj, &
             flowDomsd(nn,1,sps)%sk, &
             flowDomsd(nn,1,sps)%vol, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        deallocate(&
             flowDomsd(nn,1,sps)%rotMatrixI, &
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
             flowDomsd(nn,1,sps)%scratch, &
             stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        deallocate(flowDomsd(nn,1,sps)%p, &
             flowDomsd(nn,1,sps)%gamma, &
             flowDomsd(nn,1,sps)%aa, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        deallocate(flowDomsd(nn,1,sps)%rlv, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        deallocate(flowDomsd(nn,1,sps)%rev,stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        deallocate(&
             flowDomsd(nn,1,sps)%radI,     &
             flowDomsd(nn,1,sps)%radJ,     &
             flowDomsd(nn,1,sps)%radK,stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)

        deallocate(&
        flowDomsd(nn,1,sps)%ux, &
        flowDomsd(nn,1,sps)%uy, &
        flowDomsd(nn,1,sps)%uz, &
        flowDomsd(nn,1,sps)%vx, &
        flowDomsd(nn,1,sps)%vy, &
        flowDomsd(nn,1,sps)%vz, &
        flowDomsd(nn,1,sps)%wx, &
        flowDomsd(nn,1,sps)%wy, &
        flowDomsd(nn,1,sps)%wz, &
        flowDomsd(nn,1,sps)%qx, &
        flowDomsd(nn,1,sps)%qy, &
        flowDomsd(nn,1,sps)%qz, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)

        ! Deallocate allocated boundayr data
        do mm=1,nBocos
           deallocate(flowDomsd(nn,1,sps)%BCData(mm)%norm,stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
           
           deallocate(flowDomsd(nn,1,sps)%BCData(mm)%rface, stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)

           deallocate(flowDomsd(nn,1,sps)%BCData(mm)%F, stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)

           deallocate(flowDomsd(nn,1,sps)%BCData(mm)%dualArea, stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)

           deallocate(flowDomsd(nn,1,sps)%BCData(mm)%uSlip, stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)

           deallocate(flowDomsd(nn,1,sps)%BCData(mm)%TNS_Wall,stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
        enddo
        
        deallocate(flowDomsd(nn,1,sps)%BCData, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        if (sps==1) then
           deallocate(&
                flowDomsd(nn,1,sps)%bmti1,&
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
        
        deallocate(flowDomsd(nn,1,sps)%d2Wall, &
             stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        viscbocoLoop: do mm=1,nViscBocos
           deallocate(flowDomsd(nn,1,sps)%viscSubface(mm)%tau, stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
           
           deallocate(flowDomsd(nn,1,sps)%viscSubface(mm)%q, stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
           
        end do viscbocoLoop
        
        deallocate(flowDomsd(nn,1,sps)%viscSubFace, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
           
     end do
  end do

  ! Also dealloc winfd
  deallocate(winfd, stat=ierr)
  call EChk(ierr,__FILE__,__LINE__)
 
  ! Finally deallocate flowdomsd
  deallocate(flowdomsd,stat=ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! And the petsc vector(s)
  if (.not. wallDistanceNeeded) then 
     do sps=1, nTimeIntervalsSpectral
        call VecDestroy(xSurfVec(1, sps), ierr)
     end do
  end if

  do sps=1, nTimeIntervalsSpectral
     call VecDestroy(xSurfVecd(sps), ierr)
     call EChk(ierr,__FILE__,__LINE__)
  end do
  deallocate(xSurfVecd)

#ifndef USE_COMPLEX
  ! Deallocate reverse mode space for bcpointers
  deallocate(ww0, ww1, ww2, ww3, pp0, pp1, pp2, pp3, rlv0, rlv1, rlv2, rlv3, &
       rev0, rev1, rev2, rev3, gamma0, gamma1, gamma2, gamma3, ssi, xx)
  deallocate(ww0d, ww1d, ww2d, ww3d, pp0d, pp1d, pp2d, pp3d, rlv0d, rlv1d, rlv2d, rlv3d, &
       rev0d, rev1d, rev2d, rev3d, ssid, xxd)
#endif
derivVarsAllocated = .False.
end subroutine dealloc_derivative_values
