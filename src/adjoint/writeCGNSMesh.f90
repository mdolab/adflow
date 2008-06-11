!________________________________________________________________________
!
!     File:          writeCGNSMesh.f90
!     Author:        Joaquim R. R. A. Martins
!     Starting date: 04-06-2006
!     Last modified: 04-06-2006
!________________________________________________________________________

      subroutine writeCGNSMesh(fileName)
!________________________________________________________________________
!
!     Writes a single block mesh in CGNS format
!________________________________________________________________________

      use precision
      use blockPointers

      implicit none

      include 'cgnslib_f.h'

!
!     IO variables
!
      character fileName*32
      
!
!     Local variables
!
      integer(kind=intType) :: i, j, k
      integer(kind=intType), dimension(3,3) :: iSize
      integer(kind=intType) :: icelldim, iphysdim, ier
      integer(kind=intType) :: index_file, index_base
      integer(kind=intType) :: index_zone, index_coord
      integer(kind=intType) :: index_field, index_flow
      character basename*32, zonename*32, solname*32
      real(kind=realType), dimension(:), allocatable :: coor       
      
      integer(kind=intType) :: ll, icoor, ix, ierr
      character coorName*32

!File Parameters
      integer :: unit = 8,unit1 = 30,ierror
      character(len = 12)::outfile
      
      outfile = "cgnstest.txt"
      
      open (UNIT=unit,File=outfile,status='replace',action='write',iostat=ierror)
      if(ierror /= 0)                        &
           call terminate("verifyResiduals", &
           "Something wrong when &
           &calling open")


      write(*,*)'Calling writeCGNSMesh...'
!________________________________________________________________________
!
!     Begin execution
!________________________________________________________________________

      ! Open CGNS file
      call cg_open_f(fileName, MODE_WRITE, index_file, ier)
      ! Create base
      basename = 'Base'
      icelldim = 3
      iphysdim = 3
      call cg_base_write_f(index_file, basename, icelldim, iphysdim, &
      index_base, ier)
      ! Define zone name
      zonename = 'Zone 1'
      ! Vertex size
      iSize(1,1) = il
      iSize(2,1) = jl
      iSize(3,1) = kl
!**************
      write(*,*)'isize', iSize(:,1)
      ! Cell size
      iSize(1,2) = iSize(1,1) - 1 
      iSize(2,2) = iSize(2,1) - 1
      iSize(3,2) = iSize(3,1) - 1
      ! Boundary vertex size (always zero for structured grids)
      iSize(1,3) = 0
      iSize(2,3) = 0
      iSize(3,3) = 0
      
      ! Create zone
      call cg_zone_write_f(index_file, index_base, zonename, iSize, &
      Structured, index_zone, ier)

      ! Allocate temporary coordinates
      ll = il * jl * kl
      print *, "ll = ", ll
      allocate(coor(ll), stat=ierr)
      if (ierr /= 0) then
         print *, "ierr = 1"
         stop
      end if
      
      do ix = 1, 3
         icoor = 1
         do k = 1, kl
            do j = 1, jl
               do i = 1, il
                  !coor(icoor) = x(i+1, j+1, k+1, ix)
                  coor(icoor) = x(i, j, k, ix)
                  icoor = icoor + 1
!********************************
                  write(unit,20) coor(icoor-1),x(i+1, j+1, k+1, ix),i,j,k,ix
20                format(1x,'x ',2f18.10,4I4)
               end do
            end do
         end do
         !*****************
         write(*,*) 'Icoord',icoor
         
         select case(ix)
         case(1)
            coorName = 'CoordinateX'
         case(2)
            coorName = 'CoordinateY'
         case(3)
            coorName = 'CoordinateZ'
         end select
         
         ! Write grid coordinates
         call cg_coord_write_f(index_file, index_base, index_zone, RealDouble, &
              coorName, coor, index_coord, ierr)
      end do
      deallocate(coor)
      
      ! Close CGNS file
      call cg_close_f(index_file, ier)
      
      end subroutine writeCGNSMesh
