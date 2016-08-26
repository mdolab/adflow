subroutine getHeatFlux(hflux, npts, sps_in)
  use constants
  use blockPointers
  use flowVarRefState
  use inputTimeSpectral
  use communication
  use inputPhysics
  implicit none
  !
  !      Local variables.
  !
  integer(kind=intType), intent(in) :: npts, sps_in
  real(kind=realType), intent(out) :: hflux(npts)

  integer(kind=intType) :: mm, nn, i, j, ii, sps
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd

  !       Begin execution                                                

  sps = sps_in

  ii = 0 
  domains: do nn=1,nDom
     call setPointers(nn,1_intType,sps)
     call heatFluxes
     
     ! Loop over the number of viscous boundary subfaces of this block.
     ! According to preprocessing/viscSubfaceInfo, visc bocos are numbered
     ! before other bocos. Therefore, mm_nViscBocos == mm_nBocos
     bocos: do mm=1,nBocos

        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then
           do j=BCData(mm)%jnBeg,BCData(mm)%jnEnd
              do i=BCData(mm)%inBeg,BCData(mm)%inEnd
                 ii = ii + 1
                 hflux(ii) = BCData(mm)%sHeatFlux(i, j)
              end do
           end do
        end if

     end do bocos
  end do domains
end subroutine getHeatFlux

subroutine heatFluxes
  use constants
  use blockPointers
  use flowVarRefState
  use inputPhysics
  use bcroutines
  use costFunctions
  implicit none
  !
  !      Local variables.
  !
  integer(kind=intType) :: nn, i, j, ii
  real(kind=realType) :: fact, scaleDim
  real(kind=realType) :: qw, qA
  logical :: heatedSubface

  !
  !       Begin execution                                                
  !
  ! ! Set the actual scaling factor such that ACTUAL heat flux is computed
  ! ! The factor is determined from stanton number
  ! scaleDim = pRef*sqrt(pRef/rhoRef)

  ! ! Loop over the boundary subfaces of this block.
  ! bocos: do nn=1,nBocos
  !    ! Only include this patch if necessary
  !    mask: if(bcData(nn)%mask == 1) then 
  !       allwalls: if(BCType(nn) == EulerWall .or. &
  !            BCType(nn) == NSWallAdiabatic .or. &
  !            BCType(nn) == NSWallIsothermal) then

  !          ! Check if the face is a wall with nonzero heat flux, i.e. isothermal.
  !          heatedSubface = .true.
  !          if(BCType(nn) == EulerWall .or. &
  !               BCType(nn) == NSWallAdiabatic) heatedSubface = .false.

  !          ! Set a bunch of pointers depending on the face id to make
  !          ! a generic treatment possible. The routine setBcPointers
  !          ! is not used, because quite a few other ones are needed.
  !          call setBCPointers(nn, .True.)

  !          select case (BCFaceID(nn))
  !          case (iMin)
  !             fact = -one
  !          case (iMax)
  !             fact = one
  !          case (jMin)
  !             fact = -one
  !          case (jMax)
  !             fact = one
  !          case (kMin)
  !             fact = -one
  !          case (kMax)
  !             fact = one
  !          end select

  !          ! Loop over the quadrilateral faces of the subface. Note that
  !          ! the nodal range of BCData must be used and not the cell
  !          ! range, because the latter may include the halo's in i and
  !          ! j-direction. The offset +1 is there, because inBeg and jnBeg
  !          ! refer to nodal ranges and not to cell ranges. The loop
  !          ! (without the AD stuff) would look like:
  !          !
  !          ! do j=(BCData(nn)%jnBeg+1),BCData(nn)%jnEnd
  !          !    do i=(BCData(nn)%inBeg+1),BCData(nn)%inEnd

  !          bcData(nn)%dualArea  = zero
  !          bcData(nn)%sHeatFlux = zero
  !          do ii=0,(BCData(nn)%jnEnd - bcData(nn)%jnBeg)*(bcData(nn)%inEnd - bcData(nn)%inBeg) -1 
  !             i = mod(ii, (bcData(nn)%inEnd-bcData(nn)%inBeg)) + bcData(nn)%inBeg + 1
  !             j = ii/(bcData(nn)%inEnd-bcData(nn)%inBeg) + bcData(nn)%jnBeg + 1

  !             ! Scatter a quarter of the area to each node:
  !             qA = fourth*sqrt(ssi(i, j, 1)**2 + ssi(i, j, 2)**2 + ssi(i, j, 3)**2)
  !             BCData(nn)%dualArea(i-1, j-1) = BCData(nn)%dualArea(i-1, j-1) + qA
  !             BCData(nn)%dualArea(i  , j-1) = BCData(nn)%dualArea(i  , j-1) + qA
  !             BCData(nn)%dualArea(i-1, j  ) = BCData(nn)%dualArea(i-1, j  ) + qA
  !             BCData(nn)%dualArea(i  , j  ) = BCData(nn)%dualArea(i  , j  ) + qA
  !          enddo
  !          !
  !          !          **************************************************************
  !          !          *                                                            *
  !          !          * Integration of the viscous forces.                         *
  !          !          * Only for viscous boundaries.                               *
  !          !          *                                                            *
  !          !          **************************************************************
  !          !
  !          heatedwall: if( heatedSubface ) then
  !             ! Loop over the quadrilateral faces of the subface and
  !             ! compute nodal total heat flux (area included)
  !             do ii=0,(BCData(nn)%jnEnd - bcData(nn)%jnBeg)*(bcData(nn)%inEnd - bcData(nn)%inBeg) -1 
  !                i = mod(ii, (bcData(nn)%inEnd-bcData(nn)%inBeg)) + bcData(nn)%inBeg + 1
  !                j = ii/(bcData(nn)%inEnd-bcData(nn)%inBeg) + bcData(nn)%jnBeg + 1

  !                ! Compute the normal heat flux on the face. Inward positive.
  !                qw = -fact*scaleDim* &
  !                     sqrt(ssi(i,j,1)**2 + ssi(i,j,2)**2 + ssi(i,j,3)**2) * &
  !                     ( viscSubface(nn)%q(i,j,1)*BCData(nn)%norm(i,j,1) &
  !                     + viscSubface(nn)%q(i,j,2)*BCData(nn)%norm(i,j,2) &
  !                     + viscSubface(nn)%q(i,j,3)*BCData(nn)%norm(i,j,3))

  !                ! Scatter a quarter of the heat flux to each node:
  !                qw = fourth*qw
  !                BCData(nn)%sHeatFlux(i-1,j-1) = BCData(nn)%sHeatFlux(i-1,j-1) + qw
  !                BCData(nn)%sHeatFlux(i  ,j-1) = BCData(nn)%sHeatFlux(i  ,j-1) + qw
  !                BCData(nn)%sHeatFlux(i-1,j  ) = BCData(nn)%sHeatFlux(i-1,j  ) + qw
  !                BCData(nn)%sHeatFlux(i  ,j  ) = BCData(nn)%sHeatFlux(i  ,j  ) + qw
  !             enddo
  !          else
  !             bcData(nn)%sHeatFlux = zero
  !          end if heatedwall

  !          ! Get the nodal heat flux (per unit area)
  !          do j= BCData(nn)%jnBeg, BCData(nn)%jnEnd
  !             do i=BCData(nn)%inBeg, BCData(nn)%inEnd
  !                bcData(nn)%sHeatFlux(i, j) = &
  !                     bcData(nn)%sHeatFlux(i, j) / bcData(nn)%dualArea(i, j)
  !             end do
  !          end do

  !          call resetBCPointers(nn, .True.)
  !       end if allwalls
  !    else
  !       bcData(nn)%sHeatFlux = zero
  !    end if mask
  ! enddo bocos
end subroutine heatFluxes
