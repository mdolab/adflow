! Debugging routine for writing the dual grids along with the volume
! in CGNS format to help debug. 

subroutine writeDualMesh(fileName)


  ! This is a debugging routine for writing out meshes *as they are
  ! partioned*. This can be useful for debugging overset issues.
  use constants
  use communication, only : adflow_comm_world, myid, nProc
  use blockPointers, only : ie, je, ke, il, jl, kl, x, globalCell, vol, &
       nDom, iblank
  implicit none

  include 'cgnslib_f.h'


  character(len=*), intent(in) :: fileName
  integer(kind=intType) :: nDomTotal, iProc, nn, i, j, k, iDim, iDom, ierr, ii
  integer(kind=intType) :: iii,jjj,kkk
  integer(kind=intType) :: bufSize, maxSize, ibufSize, imaxSize
  integer(kind=intType), dimension(3, nDom) :: localDim
  integer(kind=intType), dimension(:), allocatable :: nDomProc, cumDomProc
  integer(kind=intType), dimension(:, :), allocatable :: dims
  real(kind=realType), dimension(:), allocatable :: buffer
  real(kind=realType), dimension(:, :, :, :), allocatable :: xtmp
  integer(kind=intType) :: ier, zoneCOunter, sizes(9), base, zoneID, coordID, cg, zone
  integer(kind=intType) :: ifield, iSol
  character*40 :: tmpStr, zoneName
  character*32 :: coorNames(3)
  integer status(MPI_STATUS_SIZE) 

  coorNames(1) = "CoordinateX"
  coorNames(2) = "CoordinateY"
  coorNames(3) = "CoordinateZ"

  ! Gather the dimensions of all blocks to everyone
  call mpi_allreduce(nDom, nDomTotal, 1, adflow_integer, MPI_SUM, &
       adflow_comm_world, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  ! Store the sizes of the local blocks
  do nn=1,nDom
     call setPointers(nn, 1, 1)
     localDim(1, nn) = ie
     localDim(2, nn) = je
     localDim(3, nn) = ke
  end do

  ! Allocate the space we need for the numbers and cumulative form
  allocate(nDomProc(0:nProc-1), cumDomProc(0:nProc), dims(3, nDomTotal))

  ! Receive the number of domains from each proc using an allgather.
  call mpi_allgather(nDom, 1, adflow_integer, nDomProc, 1, adflow_integer, &
       adflow_comm_world, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  ! Compute the cumulative format:
  cumDomProc(0) = 0
  do iProc=1, nProc
     cumDomProc(iProc) = cumDomProc(iProc-1) + nDomProc(iProc-1)
  end do

  ! We will also allgather all of the block sizes which will make
  ! things a little easier since everyone will know the proper sizes
  ! for the sends
  call mpi_allgatherV(localDim, nDom*3, adflow_integer, dims, 3*nDomProc, &
       3*cumDomProc, adflow_integer, adflow_comm_world, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  maxSize = 0
  do i=1,nDomTotal
     maxSize = max(maxSize, dims(1, i)*dims(2,i)*dims(3,i)*5)
  end do

  allocate(buffer(maxSize))

  if (myid == 0) then 
     call cg_open_f(fileName, mode_write, cg, ier)
     base = 1
     call cg_base_write_f(cg, "Base#1", 3, 3, base, ier)

     zoneCounter = 0
     ! Write my own blocks first
     do nn=1,nDom
        call setPointers(nn, 1, 1)
        
        sizes(1) = ie
        sizes(2) = je
        sizes(3) = ke
        sizes(4) = il
        sizes(5) = jl
        sizes(6) = kl
        sizes(7) = 0
        sizes(8) = 0
        sizes(9) = 0

999     FORMAT('domain.', I5.5)
        zoneCounter = zoneCounter + 1
        write(zonename, 999) zoneCounter 

        call cg_zone_write_f(cg, base, zonename, sizes, Structured, zoneID, ier)
        
        allocate(xtmp(sizes(1), sizes(2), sizes(3), 5))

        do k=1, ke
           do j=1, je
              do i=1, ie
                 xtmp(i,j,k,1:3) = eighth*(&
                      x(i-1, j-1, k-1, :) + &
                      x(i  , j-1, k-1, :) + &
                      x(i-1, j  , k-1, :) + &
                      x(i  , j  , k-1, :) + &
                      x(i-1, j-1, k  , :) + &
                      x(i  , j-1, k  , :) + &
                      x(i-1, j  , k  , :) + &
                      x(i  , j  , k  , :))
                 xtmp(i,j,k,4) = vol(i,j,k)
                 if (globalCell(i,j,k) >=0)  then 
                    xtmp(i,j,k,5) = dble(iblank(i,j,k))
                 else
                    xtmp(i,j,k,5) = zero
                 end if
              end do
           end do
        end do

        do idim=1, 3
           call cg_coord_write_f(cg, base, zoneID, realDouble, coorNames(idim), &
                xtmp(:, :, :, idim), coordID, ier)
        end do

        call cg_sol_write_f(cg, base, zoneID, "flowSolution", Vertex, iSol, ier)

        call cg_field_write_f(cg, base, zoneID, iSol, realDouble, "volume", &
        xtmp(:, :, :, 4), iField, ier)

        call cg_field_write_f(cg, base, zoneID, iSol, realDouble, "iBlank", &
             xtmp(:, :, :, 5), iField, ier)

        deallocate(xtmp)
     end do
     
     ! Now loop over the remaining blocks...receiving each and writing:

     do iProc=1, nProc-1
        do nn=1, nDomProc(iProc)
           iDom = cumDomProc(iProc) + nn
           bufSize = dims(1, iDom)*dims(2, iDom)*dims(3,iDom)*5

           call MPI_Recv(buffer, bufSize, adflow_real, iProc, iProc, &
                adflow_comm_world, status, ierr)

           zoneCounter = zoneCounter + 1
           write(zonename, 999) zoneCounter 
              sizes(1) = dims(1, iDom)
              sizes(2) = dims(2, iDom)
              sizes(3) = dims(3, iDom)
              sizes(4) = dims(1, iDom)-1
              sizes(5) = dims(2, iDom)-1
              sizes(6) = dims(3, iDom)-1
              sizes(7) = 0
              sizes(8) = 0
              sizes(9) = 0
           call cg_zone_write_f(cg, base, zonename, sizes, Structured, zoneID, ier)
           ii = 0
           allocate(xtmp(sizes(1), sizes(2), sizes(3), 5))
           do k=1, sizes(3)
              do j=1, sizes(2)
                 do i=1, sizes(1)
                    xtmp(i,j,k,1) = buffer(ii+1)
                    xtmp(i,j,k,2) = buffer(ii+2)
                    xtmp(i,j,k,3) = buffer(ii+3)
                    xtmp(i,j,k,4) = buffer(ii+4)
                    xtmp(i,j,k,5) = buffer(ii+5)
                    ii = ii + 5
                 end do
              end do
           end do

           do idim=1, 3
              call cg_coord_write_f(cg, base, zoneID, realDouble, coorNames(idim), &
                   xtmp(:, :, :, idim), coordID, ier)
           end do
           
           call cg_sol_write_f(cg, base, zoneID, "flowSolution", Vertex, iSol, ier)
           call cg_field_write_f(cg, base, zoneID, iSol, realDouble, "volume", &
                xtmp(:, :, :, 4), iField, ier)
           call cg_field_write_f(cg, base, zoneID, iSol, realDouble, "iBlank", &
                xtmp(:, :, :, 5), iField, ier)

           deallocate(xtmp)
        end do
     end do
  
  else 

     ! Pack and send my stuff:
     do nn=1, nDom
        call setPointers(nn, 1, 1)
        ii = 0
        do k=1, ke
           do j=1, je
              do i=1, ie
                 do iDim=1,3
                    buffer(ii+idim) =  eighth*(&
                         x(i-1, j-1, k-1, idim) + &
                         x(i  , j-1, k-1, idim) + &
                         x(i-1, j  , k-1, idim) + &
                         x(i  , j  , k-1, idim) + &
                         x(i-1, j-1, k  , idim) + &
                         x(i  , j-1, k  , idim) + &
                         x(i-1, j  , k  , idim) + &
                         x(i  , j  , k  , idim))
                 end do
                 buffer(ii+4) = vol(i,j,k)
                 if (globalCell(i,j,k) >0) then 
                    buffer(ii+5) = dble(iblank(i,j,k))
                 else
                    buffer(ii+5) = zero
                 end if
                 
                 ii = ii + 5
              end do
           end do
        end do

        call mpi_send(buffer, ii, adflow_real, 0, myid, &
             adflow_comm_world, ierr)
     end do
  end if

  deallocate(buffer, nDomProc, cumDomProc, dims)


end subroutine writeDualMesh
