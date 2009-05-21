!
! ***********************************
! *  File: integratedWarpDerivParallel.f90
! *  Author: C.A.(Sandy) Mader
! *  Started: 21-05-2009
! *  Modified: 21-05-2009
! ***********************************

subroutine integratedWarpDerivParallel!(ncoords,xyzface,indices_new)

  ! Take the current mesh and complete a warp based on the current
  ! Perturbed Coordinates

  use blockpointers
  use communication
  use mdDataLocal
  use warpingPETSc
!  use ADjointPETSc
  implicit none
  !Subroutine Arguments
!  integer(kind=intType)::ncoords
!  real(kind=realType), dimension(3,ncoords)::xyzface
!  integer(kind=intType),dimension(4,ncoords)::indices_new
  
  ! Local Arguments
  
  integer(kind=intType)::level=1,idxvol,idxsurf
  integer(kind=intType)::nn,sps=1,imax,jmax,kmax,i,j,k,n,mm,ll,nnn
  real(kind=realType), dimension(:,:,:,:),allocatable::xyznew,xyz0,xyznewd
  integer(kind=intType),dimension(6)::IFACEPTB
  integer(kind=intType),dimension(12)::IEDGEPTB
!  real(kind=realType), dimension(:,:,:,:,:,:,:),allocatable::xyzderiv

  real(kind=realType), dimension(2) :: time
  real(kind=realType)               :: timeAdjLocal, timeAdj

!!$  integer :: unitWarp = 8,ierror
!!$  character(len = 20)::outfile,testfile
!!$  
!!$  write(testfile,100) myid!12
!!$100 format (i5)  
!!$  testfile=adjustl(testfile)
!!$  write(outfile,101) trim(testfile)!testfile
!!$101 format("ADParaWarpfile",a,".out")
!!$  !outfile = "CSMachfile.txt"
!!$  unitWarp = 8+myID
!!$  !outfile = "CSMachfile.txt"
!!$  
!!$  open (UNIT=unitWarp,File=outfile,status='replace',action='write',iostat=ierror)
!!$  if(ierror /= 0)                        &
!!$       call terminate("integradtedWarpDeriv", &
!!$       "Something wrong when &
!!$       &calling open")

#ifndef USE_NO_PETSC  
  
  ! Send some feedback to screen.

  if( PETScRank==0 ) &
       write(*,10) "Assembling dXv/dXs matrix..."
  
  ! Get the initial time.
  
  call cpu_time(time(1))
  
  !zero the matrix for dXvdXs ADD call
  call MatZeroEntries(dXvdXs,PETScIerr)
  
  if( PETScIerr/=0 ) &
       call terminate("integratedWarpDerivParallel", "Error in MatZeroEntries dXvdXs")
  
  !loop over the Global surface points on this process to calculate the derivatives 
  !loop over domains
  do nn = 1,nDom
     call setPointers(nn,1,sps)
   

 !    print *,'derivative sizes',IMAX+1,JMAX+1,KMAX+1,3,ndom,ncoords,3
 !    allocate(xyzderiv(0:IMAX+1,0:JMAX+1,0:KMAX+1,3,ndom,ncoords,3))
     
     !loop over new coordinates array
     print *,'number of local surface nodes',myID,mdNGlobalSurfNodesLocal(myID+1)
     do mm = 1,mdNGlobalSurfNodesLocal(myID+1)
        !Check to see that coordinate is in this block. if so, update
        if( mdSurfGlobalIndLocal(4,mm)==nn)then
           print *,'getting sensitivites for surface point:',mm,'on Process',myID
           
           !only local block needs to be perturbed. Index sychronization will take care of the rest
           
           IMAX = IL
           JMAX = JL
           KMAX = KL
           ! SAVE NEW AND INITIAL XYZ VALUES TO BE PASSED TO WARPBLK
           
           print *,'allocate xyz0'
           ALLOCATE(XYZ0(3,0:IMAX+1,0:JMAX+1,0:KMAX+1),XYZNEW(3,0:IMAX+1,0:JMAX+1,0:KMAX+1),XYZNEWd(3,0:IMAX+1,0:JMAX+1,0:KMAX+1))
           
           
           
           do ll = 1,3
              xyz0 = 0
              xyznew = 0
              
              XYZ0(1,1:IMAX,1:JMAX,1:KMAX) = XInit(1:IMAX,1:JMAX,1:KMAX,1)
              XYZ0(2,1:IMAX,1:JMAX,1:KMAX) = XInit(1:IMAX,1:JMAX,1:KMAX,2)
              XYZ0(3,1:IMAX,1:JMAX,1:KMAX) = XInit(1:IMAX,1:JMAX,1:KMAX,3)
              XYZNEW(1,1:IMAX,1:JMAX,1:KMAX) = X(1:IMAX,1:JMAX,1:KMAX,1)
              XYZNEW(2,1:IMAX,1:JMAX,1:KMAX) = X(1:IMAX,1:JMAX,1:KMAX,2)
              XYZNEW(3,1:IMAX,1:JMAX,1:KMAX) = X(1:IMAX,1:JMAX,1:KMAX,3)
              print *,'call warp local'
              
              xyznewd(:,:,:,:)  = 0.0
              xyznewd(ll, mdSurfGlobalIndLocal(1,mm), mdSurfGlobalIndLocal(2,mm), mdSurfGlobalIndLocal(3,mm)) = 1.0
              print *,'perturbed node?', xyznewd(ll, mdSurfGlobalIndLocal(1,mm), mdSurfGlobalIndLocal(2,mm), mdSurfGlobalIndLocal(3,mm)), ll, mdSurfGlobalIndLocal(1,mm), mdSurfGlobalIndLocal(2,mm), mdSurfGlobalIndLocal(3,mm)
              !determine the explicitly and implicitly perturbed faces and edges
              !call flagImplicitEdgesAndFaces(ifaceptb,iedgeptb)
              call flagImplicitEdgesAndFacesDeriv(xyznewd,ifaceptb,iedgeptb)
              
              print*,'warpingblock',nn,ll!,indices_new(1,mm),indices_new(2,mm),indices_new(3,mm)
              call WARP_LOCAL_D(xyznew, xyznewd, xyz0, ifaceptb, iedgeptb, imax&
                   &  , jmax, kmax)
              
              print *,'assign derivatives'
              ! ASSIGN THESE derivative values
              DO I=1,IMAX
                 DO J=1,JMAX
                    DO K=1,KMAX
                       do n = 1,3
                          idxvol = globalNode(i,j,k)*3+n
                          idxsurf= mdSurfGlobalIndLocal(5,mm)*3+11
                          if (xyznewd(n,I,J,K).ne.0.0)then
                             call MatSetValues(dXvdXs, 1, idxvol-1, 1, idxsurf-1,   &
                                               xyznewd(n,I,J,K), ADD_VALUES, PETScIerr)
                             if( PETScIerr/=0 ) &
                                  print *,'matrix setting error'!call errAssemb("MatSetValues", "verifydrdw")
                          endif
                       enddo
                    END DO
                 END DO
              END DO
              
           end do
           deALLOCATE(XYZ0,XYZNEW,xyznewd)
        endif
     end do
     !deallocate(xyzderiv)
  end do
!
!     ******************************************************************
!     *                                                                *
!     * Complete the PETSc matrix assembly process.                    *
!     *                                                                *
!     ******************************************************************
!
      ! MatAssemblyBegin - Begins assembling the matrix. This routine
      !  should be called after completing all calls to MatSetValues().
      !
      ! Synopsis
      !
      ! #include "petscmat.h" 
      ! PetscErrorCode PETSCMAT_DLLEXPORT MatAssemblyBegin(Mat mat, &
      !                                            MatAssemblyType type)
      !
      ! Collective on Mat
      !
      ! Input Parameters
      !   mat  - the matrix
      !   type - type of assembly, either MAT_FLUSH_ASSEMBLY or
      !          MAT_FINAL_ASSEMBLY
      ! Notes
      ! MatSetValues() generally caches the values. The matrix is ready
      !  to use only after MatAssemblyBegin() and MatAssemblyEnd() have
      !  been called. Use MAT_FLUSH_ASSEMBLY when switching between
      !  ADD_VALUES and INSERT_VALUES in MatSetValues(); use
      !  MAT_FINAL_ASSEMBLY for the final assembly before using the
      !  matrix.
      !
      ! see .../petsc/docs/manualpages/Mat/MatAssemblyBegin.html

      call MatAssemblyBegin(dXvdXs,MAT_FINAL_ASSEMBLY,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("integratedWarpDerivParallel", &
                       "Error in MatAssemblyBegin dXvdXs")

      ! MatAssemblyEnd - Completes assembling the matrix. This routine
      !                  should be called after MatAssemblyBegin().
      !
      ! Synopsis
      !
      ! #include "petscmat.h" 
      ! PetscErrorCode PETSCMAT_DLLEXPORT MatAssemblyEnd(Mat mat,&
      !                                            MatAssemblyType type)
      !
      ! Collective on Mat
      !
      ! Input Parameters
      !   mat  - the matrix
      !   type - type of assembly, either MAT_FLUSH_ASSEMBLY or
      !          MAT_FINAL_ASSEMBLY
      !
      ! see .../petsc/docs/manualpages/Mat/MatAssemblyEnd.html

      call MatAssemblyEnd  (dXvdXs,MAT_FINAL_ASSEMBLY,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("integratedWarpDerivParallel", &
                       "Error in MatAssemblyEnd dXvdXs")

      ! Let PETSc know that the dXvdXs matrix retains the same nonzero 
      ! pattern, in case the matrix is assembled again, as for a new
      ! point in the design space.

      ! MatSetOption - Sets a parameter option for a matrix.
      !   Some options may be specific to certain storage formats.
      !   Some options determine how values will be inserted (or added).
      !   Sorted,row-oriented input will generally assemble the fastest.
      !   The default is row-oriented, nonsorted input.
      !
      ! Synopsis
      !
      ! #include "petscmat.h" 
      ! call MatSetOption(Mat mat,MatOption op,PetscErrorCode ierr)
      !
      ! Collective on Mat
      !
      ! Input Parameters
      !   mat    - the matrix
      !   option - the option, one of those listed below (and possibly
      !     others), e.g., MAT_ROWS_SORTED, MAT_NEW_NONZERO_LOCATION_ERR
      !
      ! see .../petsc/docs/manualpages/Mat/MatSetOption.html
      ! or PETSc users manual, pp.52

      call MatSetOption(dXvdXs,MAT_NO_NEW_NONZERO_LOCATIONS,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("integratedWarpDerivParallel", &
                       "Error in MatSetOption dXvdXs")

      ! Get new time and compute the elapsed time.

      call cpu_time(time(2))
      timeAdjLocal = time(2)-time(1)

      ! Determine the maximum time using MPI reduce
      ! with operation mpi_max.

      call mpi_reduce(timeAdjLocal, timeAdj, 1, sumb_real, &
                      mpi_max, 0, PETSC_COMM_WORLD, PETScIerr)

      if( PETScRank==0 ) &
        write(*,20) "Assembling dR/dx matrix time (s) =", timeAdj
!
!     ******************************************************************
!     *                                                                *
!     * Visualize the assembled matrix.                                *
!     *                                                                *
!     ******************************************************************
!
      ! MatView - Visualizes a matrix object.
      !
      ! Synopsis
      !
      ! #include "petscmat.h" 
      ! PetscErrorCode PETSCMAT_DLLEXPORT MatView(Mat mat, &
      !                                              PetscViewer viewer)
      !
      ! Collective on Mat
      !
      ! Input Parameters
      !   mat    - the matrix
      !   viewer - visualization context
      !
      ! Notes
      ! The available visualization contexts include
      !  PETSC_VIEWER_STDOUT_SELF  - standard output (default)
      !  PETSC_VIEWER_STDOUT_WORLD - synchronized standard output where
      !                         only the first processor opens the file.
      !                         All other processors send their data to
      !                         the first processor to print.
      !  PETSC_VIEWER_DRAW_WORLD- graphical display of nonzero structure
      !
      ! see .../petsc/docs/manualpages/Mat/MatView.html
      ! or PETSc users manual, pp.57,148

      if( debug ) then
        !call MatView(dXvdXs,PETSC_VIEWER_DRAW_WORLD,PETScIerr)
        call MatView(dXvdXs,PETSC_VIEWER_STDOUT_WORLD,PETScIerr)
        if( PETScIerr/=0 ) &
          call terminate("integratedWarpDerivParallel", "Error in MatView")
        !pause
      endif

      ! Flush the output buffer and synchronize the processors.

      call f77flush()
      call mpi_barrier(PETSC_COMM_WORLD, PETScIerr)


!      !now extract and write to a file
!                          write(unitWarp,13) XYZNEWd(n,I,J,K),i,j,k,n,nnn,nn,mm,ll
!13                        format(1x,'WarpSurf',f18.10,8I8)
!


      ! Output formats.

   10 format(a,1x,i3,1x,a,1x,i3)
   20 format(a,1x,f8.2)

  print *,'warp derivativesparallel finished'
#endif
end subroutine integratedWarpDerivParallel
   
