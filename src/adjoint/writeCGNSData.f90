!________________________________________________________________________
!
!     File:          writeCGNSData.f90
!     Author:        Joaquim R. R. A. Martins
!     Starting date: 04-06-2006
!     Last modified: 04-07-2006
!________________________________________________________________________

      subroutine writeCGNSData(fileName, dataName, data)
!________________________________________________________________________
!
!     Writes a single block mesh in CGNS format
!________________________________________________________________________

      use precision
      use blockPointers
      use flowvarrefstate

      implicit none

      include 'cgnslib_f.h'

!
!     IO variables
!
      real(kind=realType), dimension(nx,ny,nz,nw) :: data
      character fileName*32, dataName*32
      
!
!     Local variables
!
      integer(kind=intType) :: iw, ierr
      integer(kind=intType) :: indexFile, indexField, indexFlow
      integer(kind=intType) :: indexBase, indexZone
      character solName*32, dataNameTmp*34, nwc*2

!________________________________________________________________________
!
!     Begin execution
!________________________________________________________________________


      ! Open CGNS file
      call cg_open_f(fileName, MODE_MODIFY, indexFile, ierr)

      solName = "Solution"
      indexBase = 1
      indexZone = 1

      ! Create flow solution node
      call cg_sol_write_f(indexFile, indexBase, indexZone, solName, &
                          CellCenter, indexFlow, ierr)
     ! call cg_sol_write_f(indexFile, indexBase, indexZone, solName, &
     !                     Vertex, indexFlow, ierr)

      do iw = 1,nw
         write (nwc, '(I1)') iw
         dataNameTmp = trim(dataName)//nwc
         call cg_field_write_f(indexFile, indexBase, indexZone, & 
                               indexFlow, RealDouble, &
                               dataNameTmp, data(:,:,:,iw), &
                               indexField, ierr)
      end do

      ! Close CGNS file
      call cg_close_f(indexFile, ierr)
            
      end subroutine writeCGNSData
