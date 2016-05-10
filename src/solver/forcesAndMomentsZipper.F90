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
  integer(kind=intType) :: i, j, k, l, ii, ierr, nn, sps, iTri, nl, ng
  real(kind=realType), dimension(:), pointer :: xx, pp, vv, xVolume
  real(kind=realType), dimension(3) :: x1, x2, x3, xc
  real(kind=realType), dimension(3) :: pp1, pp2, pp3, pp4, pres
  real(kind=realType), dimension(3) :: vv1, vv2, vv3, vv4, ss, vis, norm
  real(kind=realType), dimension(3) :: Fp, Fv, Mp, Mv, Ftmp, Mtmp
  real(kind=realType) :: fact, scaleDim, u, v
  real(kind=realType), dimension(3) :: refPoint
  real(kind=realType) :: areaSum
  character(60) :: fout
  ! ------------------------------------------------------------------

  ! Set the actual scaling factor such that ACTUAL forces are computed
  scaleDim = pRef/pInf

  ! Determine the reference point for the moment computation in
  ! meters.

  refPoint(1) = LRef*pointRef(1)
  refPoint(2) = LRef*pointRef(2)
  refPoint(3) = LRef*pointRef(3)

  ! Initialize the force and moment coefficients to 0 

  cFp(1:3) = zero
  cFv(1:3) = zero
  cMp(1:3) = zero
  cMv(1:3) = zero

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

  ! Perform assembly on global force vectors
  call VecAssemblyBegin(globalPressureTractions, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecAssemblyEnd(globalPressureTractions, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecAssemblyBegin(globalViscousTractions, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecAssemblyEnd(globalViscousTractions, ierr)
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

  !! Test output
  !! ---------------------
  !write(fout,'(i7)')myid
  !fout = 'vec.'//trim(adjustl(fout))
  !fout = trim(adjustl(fout))
  !call PetscViewerASCIIOpen(sumb_comm_world,fout,viewer,ierr)
  !call VecView(zipperViscousTractions,viewer,ierr)
  !call PetscViewerDestroy(viewer,ierr)

  !call VecGetLocalSize(zipperViscousTractions,nl,ierr)
  !print*,'myid, nl ', myid, nl
  !! ---------------------

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
     ! or length of zipperPressuretraction / 3

     areaSum = zero

     do i=1, size(xx)/9
     !do i=1, size(pp)/3

        x1 = xx((i-1)*9+1:(i-1)*9+3)
        x2 = xx((i-1)*9+4:(i-1)*9+6)
        x3 = xx((i-1)*9+7:(i-1)*9+9)

        ! Just use cell center values
        pres(:) = pp((i-1)*3+1:i*3)
        vis(:) = vv((i-1)*3+1:i*3)
        
        ! Compute area
        call cross_prod(x2-x1, x3-x1, norm)
        ss = half * norm

        areaSum = areaSum + norm2(ss)

        ! Compute cell center
        xc(:) = third*(x1 + x2 + x3) - refPoint(:)
        
        ! Add to Fp and Fv
        Fp = Fp + pres * norm2(ss)
        Fv = Fv + vis  * norm2(ss)

        ! Add to Mp and Mv
        Ftmp = pres * norm2(ss)
        call cross_prod(xc, Ftmp, Mtmp)
        Mp = Mp + Mtmp

        Ftmp = vis * norm2(ss)
        call cross_prod(xc, Ftmp, Mtmp)
        Mv = Mv + Mtmp

     end do
     !write(7000,*)areaSum

     !! Debug triangles
     !! ----------------------------------------------------------
     !open(unit=101, file='tempTri.dat', form='formatted')
     !rewind(101)
     !write(101,*) 'TITLE = "Triangles"'
     !write(101,*) 'Variables = "X", "Y", "Z"'
     !write(101,*) "Zone T=Triangles_petsc"
     !write (101,*) "Nodes = ", size(xx)/3, " Elements= ", size(xx)/9, " ZONETYPE=FETRIANGLE"
     !write (101,*) "DATAPACKING=POINT"
     !! Node data
     !do i=1, size(xx)/9
     !   x1 = xx((i-1)*9+1:(i-1)*9+3)
     !   x2 = xx((i-1)*9+4:(i-1)*9+6)
     !   x3 = xx((i-1)*9+7:(i-1)*9+9)
     !   write(101,'(3(E20.12,x))') x1(1:3)
     !   write(101,'(3(E20.12,x))') x2(1:3)
     !   write(101,'(3(E20.12,x))') x3(1:3)
     !end do
     !! Cell data (conn)
     !do i=1, size(xx)/9
     !   write(101, '(3(i5,x))') (i-1)*3+1, (i-1)*3+2, (i-1)*3+3
     !end do
     !close(101)
     !! ----------------------------------------------------------

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
