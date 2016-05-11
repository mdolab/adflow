subroutine setIblankArray(level, sps)

  use blockPointers
  use communication
  implicit none

  ! Input/Output
  integer(kind=intType), intent(in) :: level, sps
  
  ! Working
  integer(kind=intType) :: i, j, k, nn
  logical :: isHole, isFlooded, isFloodSeed
  integer(kind=intType) :: nCompute, nFringe, nBlank, nFloodSeed, nFlooded
  integer(kind=intType) :: counts(5), ierr
  nCompute = 0
  nFringe = 0
  nBlank = 0
  nFloodSeed = 0
  nFlooded = 0

  do nn=1, nDom
     call setPointers(nn, level, sps)
     iBlank = 1
     do k=2, kl
        do j=2, jl
           do i=2, il

              if (fringes(i, j, k)%donorProc /= -1) then
                 iblank(i, j, k) = -1
                 nFringe = nFringe + 1

              else if (isFloodSeed(fringes(i, j, k)%status)) then
                 iBlank(i, j, k) = -3
                 nFloodSeed = nFloodSeed + 1

              else if (isFlooded(fringes(i, j, k)%status)) then
                 iBlank(i, j, k) = -2
                 nFlooded = nFlooded + 1
           
              else if (isHole(fringes(i, j, k)%status)) then 
                 iBlank(i, j, k) = 0
                 nBlank = nBlank + 1

              else
                 ! Compute cell
                 nCompute = nCompute + 1
              end if

           end do
        end do
     end do
  end do

  ! Update the iblank info. 
  domainLoop:do nn=1, nDom
     flowDoms(nn, level, sps)%intCommVars(1)%var => &
          flowDoms(nn, level, sps)%iblank(:, :, :)
  end do domainLoop
  
  ! Run the generic integer exchange
  call wHalo1to1IntGeneric(1, level, sps, commPatternCell_2nd, internalCell_2nd)

  call mpi_reduce((/nCompute, nFringe, nBlank, nFlooded, nFloodSeed/), &
       counts, 5, sumb_integer, MPI_SUM, 0, sumb_comm_world, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  if (myid == 0) then 
     print *, '+--------------------------------+'
     print *, '| Compute   Cells:', counts(1)
     print *, '| Fringe    Cells:', counts(2)
     print *, '| Blanked   Cells:', counts(3)
     print *, '| Flooded   Cells:', counts(4)
     print *, '| FloodSeed Cells:', counts(5)
     print *, '+--------------------------------+' 
  end if
end subroutine setIblankArray

subroutine dumpIblank(level, sps)

  use blockPointers
  use communication
  implicit none

  ! Input/Output
  integer(kind=intType), intent(in) :: level, sps
  
  ! Working
  integer(kind=intType) :: i, j, k, nn
  real(kind=realType) :: xp(3)
  character(80) :: fileName
  
  write (fileName,"(a,I2.2,a)") "proc_", myid, ".dat"
  open(unit=19,file=trim(fileName),form='formatted')

  do nn=1, nDom
     call setPointers(nn, level, sps)
     do k=2, kl
        do j=2, jl
           do i=2, il
              
              ! Compute the cell center:
              xp = eighth*(&
                   x(i-1, j-1, k-1, :) + &
                   x(i  , j-1, k-1, :) + &
                   x(i-1, j  , k-1, :) + &
                   x(i  , j  , k-1, :) + &
                   x(i-1, j-1, k  , :) + &
                   x(i  , j-1, k  , :) + &
                   x(i-1, j  , k  , :) + &
                   x(i  , j  , k  , :))

              write(19, "(E18.10, E18.10, E18.10, I3)"), xp(1), xp(2), xp(3), iblank(i, j, k)
           end do
        end do
     end do
  end do

  close(19)

end subroutine dumpIblank
