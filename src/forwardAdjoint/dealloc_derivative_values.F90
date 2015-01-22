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

        ! Reset w and dw                            
        do l=1,nw
           do k=0,kb 
              do j=0,jb
                 do i=0,ib
                    w(i,j,k,l) = flowdomsd(nn,1,sps)%wtmp(i,j,k,l)
                    dw(i,j,k,l) = flowdomsd(nn,1,sps)%dwtmp(i,j,k,l)
                 end do
              end do
           end do
        end do
          
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
        
        deallocate(flowDomsd(nn,1,sps)%dtl, &
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
           
           deallocate(flowDomsd(nn,1,sps)%BCData(mm)%Fp, stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
           
           deallocate(flowDomsd(nn,1,sps)%BCData(mm)%Fv, stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
           
           deallocate(flowDomsd(nn,1,sps)%BCData(mm)%M, stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)

           deallocate(flowDomsd(nn,1,sps)%BCData(mm)%sepSensor, stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
 
           deallocate(flowDomsd(nn,1,sps)%BCData(mm)%Cavitation, stat=ierr)
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
     call VecDestroy(xSurfVec(1), ierr)
  end if

  call VecDestroy(xSurfVecd, ierr)
  call EChk(ierr,__FILE__,__LINE__)
end subroutine dealloc_derivative_values
