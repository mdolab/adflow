subroutine getArea(totalArea)

  use BCTypes
  use blockPointers
  use cgnsGrid
  use communication
  use flowVarRefState
  use mdDataLocal
  use block
  implicit none
  !
  !      Local variables.
  !
  integer :: ierr
  real(kind=realType)   :: totalArea
  integer(kind=intType) :: mm, nn, i, j
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd

  real(kind=realType), dimension(:,:,:), pointer :: ss
  real(kind=realType) ::  a
  integer(kind=intType) :: iset(3),size,highInd,lowInd,ndof

  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************

  a = 0.0
  domains: do nn=1,nDom
     call setPointers(nn,1_intType,1_intType)

     ! Loop over the number of boundary subfaces of this block.
     bocos: do mm=1,nBocos

        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or.&
             BCType(mm) == NSWallIsothermal) then
           select case (BCFaceID(mm))

           case (iMin)
              ss => si( 1,:,:,:)
           case (iMax)
              ss => si(il,:,:,:)
           case (jMin)
              ss => sj(:, 1,:,:)
           case (jMax)
              ss => sj(:,jl,:,:)
           case (kMin)
              ss => sk(:,:, 1,:)
           case (kMax)
              ss => sk(:,:,kl,:)
           end select

           jBeg = BCData(mm)%jnBeg + 1; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg + 1; iEnd = BCData(mm)%inEnd

           do j=jBeg, jEnd ! This is a face loop
              do i=iBeg, iEnd ! This is a face loop 
                 a = a + sqrt(ss(i,j,1)*ss(i,j,1) + ss(i,j,2)*ss(i,j,2) + ss(i,j,3)*ss(i,j,3))
              end do
           end do
        end if
     end do bocos
  end do domains

  call MPI_Allreduce (a,totalArea,1,sumb_real,MPI_SUM,SUMB_COMM_WORLD,ierr)
 
end subroutine getArea

