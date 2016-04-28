subroutine forcesAndMomentsZipper(cFp, cFv, cMp, cMv)

  use communication
  use blockPointers
  use BCTypes
  use flowVarRefState
  use inputPhysics
  use costFunctions
  use overset
  use inputTimeSpectral

  implicit none

  real(kind=realType), dimension(3), intent(out) :: cFp, cFv
  real(kind=realType), dimension(3), intent(out) :: cMp, cMv
  integer(kind=intType) :: i, j, k, l, ii, ierr, nn, sps
  real(kind=realType), dimension(:), pointer :: xx, pp, vv, xVolume
  real(kind=realType), dimension(3) :: x1, x2, x3, pp1, pp2, pp3, pp4
  real(kind=realType), dimension(3) :: vv1, vv2, vv3, vv4, ss
  real(kind=realType), dimension(3) :: Fp, Fv, Mp, Mv
  real(kind=realType) :: fact, scaleDim

  ! Fill up xVolumeVec 
  call VecGetArrayF90(globalNodes, xVolume, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  ii = 0
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointers(nn, 1, 1)
        do k=1, kl
           do j=1, jl
              do i=1, il
                 do l= 1,3
                    ii = ii + 1
                    xVolume(ii) = X(i, j, k, l)
                 end do
              end do
           end do
        end do
     end do
  end do
  call vecRestoreArrayF90(globalNodes, xVolume, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  ! Perform the scatter from the global x vector to zipperNodes
  call VecScatterBegin(nodeZipperScatter, globalNodes, &
       zipperNodes, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterEnd(nodeZipperScatter, globalNodes, &
       zipperNodes, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  
  ! Perform the scatter from the global surface tractions to the
  ! zipper tractions

  call VecScatterBegin(tractionZipperScatter, globalPressureTractions, &
       zipperPressureTractions, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterEnd(tractionZipperScatter, globalPressureTractions, &
       zipperPressureTractions, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)


  call VecScatterBegin(tractionZipperScatter, globalViscousTractions, &
       zipperViscousTractions, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterEnd(tractionZipperScatter, globalViscousTractions, &
       zipperViscousTractions, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)


  ! Puck out pointers for the nodes and tractions

  if (myid == 0) then 
     Fp = zero
     Fv = zero
     Mp = zero
     Mv = zero
     

     call VecGetArrayF90(zipperPressureTractions, pp, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecGetArrayF90(zipperViscousTractions, vv, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecGetArrayF90(zipperNodes, xx, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     ! Number of triangles is the length of zipperNodes /9
     do i=1, size(xx)/12

        x1 = xx((i-1)*9+1:(i-1)*9+3)
        x2 = xx((i-1)*9+4:(i-1)*9+6)
        x3 = xx((i-1)*9+6:(i-1)*9+9)

        vv1 = vv((i-1)*12+1:(i-1)*12+3)
        vv2 = vv((i-1)*12+4:(i-1)*12+6)
        vv3 = vv((i-1)*12+7:(i-1)*12+9)
        vv4 = vv((i-1)*12+10:(i-1)*12+12)

        ! Using the uv...stored somewhere 
        !v = zero ! (1-uv(1))*v1 etc

        pp1 = pp((i-1)*12+1:(i-1)*12+3)
        pp2 = pp((i-1)*12+4:(i-1)*12+6)
        pp3 = pp((i-1)*12+7:(i-1)*12+9)
        pp4 = pp((i-1)*12+10:(i-1)*12+12)

        ! Using the uv...stored somewhere 
        !p = zero ! (1-uv(1))*v1 etc

        
        ! Compute area
        

        ! Add to Fp and Fv
        !Fp = Fp + 
        !Fv = Fv + 

        !Mp = Mp + 
        !Mv = Mv + 

     end do

     call VecRestoreArrayF90(zipperPressureTractions, pp, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     

     call VecRestoreArrayF90(zipperViscousTractions, vv, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     call VecRestoreArrayF90(zipperNodes, xx, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     

     fact = two/(gammaInf*pInf*MachCoef*MachCoef &
          *surfaceRef*LRef*LRef*scaleDim)

     cFp(:) =  fact*Fp
     cFv(:) =  fact*Fv

     fact = fact/(lengthRef*LRef)
     cMp(:) = Mp(:)*fact
     cMv(:) = Mv(:)*fact

  end if

end subroutine forcesAndMomentsZipper
