! This is a special function that is sued to dealloc derivative values
! in blockpointers_d for use with the AD code.

subroutine dealloc_derivative_values_bwd(level)

  use blockPointers_b ! This modules includes blockPointers
  use inputtimespectral
  use flowvarrefstate
  use inputPhysics
  use cgnsGrid
  use BCTypes
  use wallDistanceData
  implicit none

  ! Input Parameters
  integer(kind=intType) :: level

  ! Local variables
  integer(kind=intType) :: sps, ierr, i, j, k, l, mm, nn

  do nn=1,nDom
     ! Deallocate the color array in flowDoms
     deallocate(flowDomsb(nn,1,1)%color,stat=ierr)
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
                    w(i,j,k,l) = flowdomsb(nn,1,sps)%wtmp(i,j,k,l)
                    dw(i,j,k,l) = flowdomsb(nn,1,sps)%dwtmp(i,j,k,l)
                 end do
              end do
           end do
        end do
          
        ! Deallocate memtory
        deallocate(&
             flowDomsb(nn,1,sps)%wtmp, &
             flowDomsb(nn,1,sps)%dwtmp, &
             flowDomsb(nn,1,sps)%dwtmp2, &
             flowDomsb(nn,1,sps)%xtmp, &
             flowDomsb(nn,1,sps)%w_deriv, &
             flowDomsb(nn,1,sps)%dw_deriv, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
     end do deallocatespectral

     do sps=1,nTimeIntervalsSpectral
        call setPointers(nn,level,sps)

        deallocate(flowDomsb(nn,1,sps)%x, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        deallocate(flowDomsb(nn,1,sps)%si, &
             flowDomsb(nn,1,sps)%sj, &
             flowDomsb(nn,1,sps)%sk, &
             flowDomsb(nn,1,sps)%vol, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        deallocate(&
             flowDomsb(nn,1,sps)%rotMatrixI, &
             flowDomsb(nn,1,sps)%rotMatrixJ, &
             flowDomsb(nn,1,sps)%rotMatrixK,stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        deallocate(flowDomsb(nn,1,sps)%s,      &
             flowDomsb(nn,1,sps)%sFaceI, &
             flowDomsb(nn,1,sps)%sFaceJ, &
             flowDomsb(nn,1,sps)%sFaceK, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)

        deallocate(flowDomsb(nn,1,sps)%w, &
             flowDomsb(nn,1,sps)%dw, &
             flowDomsb(nn,1,sps)%fw, &
             stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        deallocate(flowDomsb(nn,1,sps)%p, &
             flowDomsb(nn,1,sps)%gamma, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        deallocate(flowDomsb(nn,1,sps)%rlv, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        deallocate(flowDomsb(nn,1,sps)%rev,stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        deallocate(flowDomsb(nn,1,sps)%dtl, &
             flowDomsb(nn,1,sps)%radI,     &
             flowDomsb(nn,1,sps)%radJ,     &
             flowDomsb(nn,1,sps)%radK,stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        ! Deallocate allocated boundayr data
        do mm=1,nBocos
           deallocate(flowDomsb(nn,1,sps)%BCData(mm)%norm,stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
           
           deallocate(flowDomsb(nn,1,sps)%BCData(mm)%rface, stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
           
           deallocate(flowDomsb(nn,1,sps)%BCData(mm)%Fp, stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
           
           deallocate(flowDomsb(nn,1,sps)%BCData(mm)%Fv, stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
           
           deallocate(flowDomsb(nn,1,sps)%BCData(mm)%M, stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)

           deallocate(flowDomsb(nn,1,sps)%BCData(mm)%sepSensor, stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)

           deallocate(flowDomsb(nn,1,sps)%BCData(mm)%Cavitation, stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)

           deallocate(flowDomsb(nn,1,sps)%BCData(mm)%uSlip, stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)

           deallocate(flowDomsb(nn,1,sps)%BCData(mm)%TNS_Wall,stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
        enddo
        
        deallocate(flowDomsb(nn,1,sps)%BCData, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        if (sps==1) then
           deallocate(&
                flowDomsb(nn,1,sps)%bmti1,&
                flowDomsb(nn,1,sps)%bmti2,&
                flowDomsb(nn,1,sps)%bmtj1,&
                flowDomsb(nn,1,sps)%bmtj2,&
                flowDomsb(nn,1,sps)%bmtk1,&
                flowDomsb(nn,1,sps)%bmtk2,&
                flowDomsb(nn,1,sps)%bvti1,&
                flowDomsb(nn,1,sps)%bvti2,&
                flowDomsb(nn,1,sps)%bvtj1,&
                flowDomsb(nn,1,sps)%bvtj2,&
                flowDomsb(nn,1,sps)%bvtk1,&
                flowDomsb(nn,1,sps)%bvtk2,&
                stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
        end if
        
        deallocate(flowDomsb(nn,1,sps)%d2Wall, &
             stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        viscbocoLoop: do mm=1,nviscBocos
           deallocate(flowDomsb(nn,1,sps)%viscSubface(mm)%tau, stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
           
           deallocate(flowDomsb(nn,1,sps)%viscSubface(mm)%q, stat=ierr)
           call EChk(ierr,__FILE__,__LINE__)
           
        end do viscbocoLoop
        
        deallocate(flowDomsb(nn,1,sps)%viscSubFace, stat=ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
     end do
  end do


  ! Also dealloc winfd
  deallocate(winfb, stat=ierr)
  call EChk(ierr,__FILE__,__LINE__)
 
  ! Finally deallocate flowdomsb
  deallocate(flowdomsb,stat=ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! And the petsc vector
  call VecDestroy(xSurfVecb, ierr)
  call EChk(ierr,__FILE__,__LINE__)

end subroutine dealloc_derivative_values_bwd
