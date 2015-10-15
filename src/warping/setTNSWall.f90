subroutine setTNSWall(tnsw, npts, sps_in)

  use BCTypes
  use blockPointers
  use flowVarRefState
  use inputTimeSpectral
  use communication
  use inputPhysics
  implicit none
  !
  !      Arguments.
  !
  integer(kind=intType), intent(in) :: npts, sps_in
  real(kind=realType), intent(in) :: tnsw(npts)
  !
  !      Local variables.
  !
  integer(kind=intType) :: mm, nn, i, j, ii, sps
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd

  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************

  sps = sps_in

  ii = 0 
  domains: do nn=1,nDom
     call setPointers(nn,1_intType,sps)

     ! Loop over the number of viscous boundary subfaces of this block.
     ! According to preprocessing/viscSubfaceInfo, visc bocos are numbered
     ! before other bocos. Therefore, mm_nViscBocos == mm_nBocos
     bocos: do mm=1,nBocos
        mask: if(bcData(mm)%mask == 1) then 

           if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
                BCType(mm) == NSWallIsothermal) then
              do j=BCData(mm)%jnBeg,BCData(mm)%jnEnd
                 do i=BCData(mm)%inBeg,BCData(mm)%inEnd
                    ii = ii + 1
                    BCData(mm)%TNS_Wall(i,j) = tnsw(ii)/TRef
                 end do
              end do
           end if

        end if mask
     end do bocos
  end do domains

  ! Interpolate TNS_Wall to coarse grids
  ! call setBCDataCoarseGrid

end subroutine setTNSWall
