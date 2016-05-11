subroutine updateXSurf(level)

  use wallDistanceData
  use blockPointers
  use inputTimeSpectral
  implicit none

  ! Input Parameters
  integer(kind=intType), intent(in) :: level

  ! Working Parameters
  integer(kind=intType) :: ii, i,j,k,l, nn, sps, ierr

  ! Fill up xVolumeVec 
  call VecGetArrayF90(xVolumeVec(level), xVolume, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  ii = 0
  do nn=1,nDom
     do sps=1,nTimeIntervalsSpectral
        call setPointers(nn, level, sps)
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
  call vecRestoreArrayF90(xVolumeVec(level), xVolume, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  ! Perform the scatter from the global x vector to xSurf. SPS loop since the xSurfVec is done by SPS instance.
  do sps=1, nTimeIntervalsSpectral
     call VecScatterBegin(wallScatter(level, sps), xVolumeVec(level), &
          xSurfVec(level, sps), INSERT_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     call VecScatterEnd(wallScatter(level, sps), xVolumeVec(level), &
          xSurfVec(level, sps), INSERT_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr,__FILE__,__LINE__)
  end do

end subroutine updateXSurf
