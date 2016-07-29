subroutine exchangeSurfaceDelta(level, sps, commPattern, internal)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * ExchangeSurfaceDelta exchanges surface delta to fill up halo   *
  !      * surface cells from adjacent blocks.                            *
  !      *                                                                *
  !      ******************************************************************
  !
  use blockPointers
  use bctypes
  use communication
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: level, sps

  type(commType),          dimension(*), intent(in) :: commPattern
  type(internalCommType), dimension(*), intent(in) :: internal
  
  ! Local
  integer(kind=intType) :: i, j, k, ii, nn, mm
  real(kind=realType), dimension(:), allocatable :: pSave
  real(kind=realType), dimension(:, :), pointer :: deltaPtr

  ! Just cheat by exchangint pressure. saving Pressure, dumping deltaPtr into the pressure,
  ! exchanging that and then restoring the pressure

  do nn=1, nDom
     call setPointers(nn, level, sps)

     ! Allocate pointer space for the integer flag communication
     allocate(flowDoms(nn, level, sps)%realCommVars(1)%var(1:ib+1, 1:jb+1, 1:kb+1))

     ! Push the surface iblank back to the generic volume variable rVar1
     bocoLoop: do mm=1, nBocos
        wallType: if (BCType(mm) == EulerWall .or. &
             BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then
           
           select case (BCFaceID(mm))
           case (iMin)
              deltaPtr => flowDoms(nn, level, sps)%realCommVars(1)%var(2+1, :, :)
           case (iMax)
              deltaPtr => flowDoms(nn, level, sps)%realCommVars(1)%var(il+1, :, :)
           case (jMin)
              deltaPtr => flowDoms(nn, level, sps)%realCommVars(1)%var(:, 2+1,  :)
           case (jMax)
              deltaPtr => flowDoms(nn, level, sps)%realCommVars(1)%var(:, jl+1, :)
           case (kMin)
              deltaPtr => flowDoms(nn, level, sps)%realCommVars(1)%var(:, :, 2+1 )
           case (kMax)
              deltaPtr => flowDoms(nn, level, sps)%realCommVars(1)%var(:, :, jl+1)
           end select
           
           ! NO HALOS!
           do j=BCData(mm)%jnBeg+1, BCData(mm)%jnEnd
              do i=BCData(mm)%inBeg+1, BCData(mm)%inEnd
                 
                 ! Remember to account for the pointer offset since
                 ! the iblank starts at zero
                 deltaPtr(i+1, j+1) = BCData(mm)%delta(i, j)
              end do
           end do
        end if wallType
     end do bocoLoop
  end do
  
  ! Exchange the variane
  call whalo1to1RealGeneric(1, level, sps, commPatternCell_1st, internalCell_1st)

  ! Copy back out
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
              deltaPtr => flowDoms(nn, level, sps)%realCommVars(1)%var(2+1, :, :)
           case (iMax)
              deltaPtr => flowDoms(nn, level, sps)%realCommVars(1)%var(il+1, :, :)
           case (jMin)
              deltaPtr => flowDoms(nn, level, sps)%realCommVars(1)%var(:, 2+1,  :)
           case (jMax)
              deltaPtr => flowDoms(nn, level, sps)%realCommVars(1)%var(:, jl+1, :)
           case (kMin)
              deltaPtr => flowDoms(nn, level, sps)%realCommVars(1)%var(:, :, 2+1 )
           case (kMax)
              deltaPtr => flowDoms(nn, level, sps)%realCommVars(1)%var(:, :, jl+1)
           end select
           
           ! INCLUDE THE HALOS!
           do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
              do i=BCData(mm)%icBeg, BCData(mm)%icEnd

                 ! Remember to account for the pointer offset since
                 ! the iblank starts at zero
                 BCData(mm)%delta(i,j) = deltaPtr(i+1, j+1)
               end do
            end do
         end if wallType2
      end do bocoLoop2

      ! Now deallocate this pointer
      deallocate(flowDoms(nn, level, sps)%realCommVars(1)%var)
   end do
 end subroutine exchangeSurfaceDelta
