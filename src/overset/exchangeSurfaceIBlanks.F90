subroutine exchangeSurfaceIblanks(level, sps, commPattern, internal)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * ExchangeIblank exchanges the 1 to 1 internal halo's for the    *
  !      * given level and sps instance.                                  *
  !      *                                                                *
  !      ******************************************************************
  !
  use constants
  use blockPointers
  use communication
  use utils, only : setPointers
  use haloExchange, only : whalo1to1intgeneric
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: level, sps

  type(commType),          dimension(*), intent(in) :: commPattern
  type(internalCommType), dimension(*), intent(in) :: internal
  
  ! Local
  integer(kind=intType) :: i, j, k, ii, nn, mm
  integer(kind=intType), dimension(:), allocatable :: iBlankSave
  integer(kind=intType), dimension(:, :), pointer :: ibp

  ! Just cheat by saving iBlank iblank array, resusing itand 
  ii = 0
  do nn=1, nDom
     call setPointers(nn, level, sps)
     ii = ii + (ib+1)*(jb+1)*(kb+1)
  end do

  allocate(iBlankSave(ii))
  ii = 0
  do nn=1, nDom
     call setPointers(nn, level, sps)
     do k=0,kb
        do j=0,jb
           do i=0,ib
              ii =ii + 1
              iBlankSave(ii) = iblank(i,j,k)
              iblank(i,j,k) = 1
           end do
        end do
     end do
     
     ! Push the surface iblank back to the volume:
     bocoLoop: do mm=1, nBocos
        wallType: if (BCType(mm) == EulerWall .or. &
             BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then
           
           select case (BCFaceID(mm))
           case (iMin)
              ibp => iblank(2, :, :)
           case (iMax)
              ibp => iblank(il, :, :)
           case (jMin)
              ibp => iblank(:, 2, :)
           case (jMax)
              ibp => iblank(:, jl, :)
           case (kMin)
              ibp => iblank(:, :, 2)
           case (kMax)
              ibp => iblank(:, :, kl)
           end select

           ! NO HALOS!
           do j=BCData(mm)%jnBeg+1, BCData(mm)%jnEnd
              do i=BCData(mm)%inBeg+1, BCData(mm)%inEnd

                 ! Remember to account for the pointer offset since
                 ! the iblank starts at zero
                 ibp(i+1, j+1) = BCData(mm)%iBlank(i,j)
              end do
           end do
        end if wallType
     end do bocoLoop
  end do
     
  ! Exchange iblanks
  domainLoop:do nn=1, nDom
     flowDoms(nn, level, sps)%intCommVars(1)%var => &
          flowDoms(nn, level, sps)%iblank(:, :, :)
  end do domainLoop
  
  ! Run the generic integer exchange
  call wHalo1to1IntGeneric(1, level, sps, commPattern, internal)

  ii = 0
  do nn=1, nDom
     call setPointers(nn, level, sps)

     ! Extract the surface iblank from the volume. 
     bocoLoop2: do mm=1, nBocos
        wallType2: if (BCType(mm) == EulerWall .or. &
             BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then
           
           select case (BCFaceID(mm))
           case (iMin)
              ibp => iblank(2, :, :)
           case (iMax)
              ibp => iblank(il, :, :)
           case (jMin)
              ibp => iblank(:, 2, :)
           case (jMax)
              ibp => iblank(:, jl, :)
           case (kMin)
              ibp => iblank(:, :, 2)
           case (kMax)
              ibp => iblank(:, :, kl)
           end select

           ! INCLUDE THE HALOS!
           do j=BCData(mm)%jnBeg, BCData(mm)%jnEnd+1
              do i=BCData(mm)%inBeg, BCData(mm)%inEnd+1
                 ! Remember to account for the pointer offset since
                 ! the iblank starts at zero
                 BCData(mm)%iBlank(i,j) = ibp(i+1, j+1)
               end do
            end do
         end if wallType2
      end do bocoLoop2

     ! Restore the saved array
     do k=0,kb
        do j=0,jb
           do i=0,ib
              ii =ii + 1
              iBlank(i,j,k) = iBlankSave(ii)
           end do
        end do
     end do
  end do
  deallocate(iblankSave)
end subroutine exchangeSurfaceIblanks
