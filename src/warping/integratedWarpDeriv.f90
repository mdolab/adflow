!
! ***********************************
! *  File: integratedWarpDeriv.f90
! *  Author: C.A.(Sandy) Mader
! *  Started: 19-05-2009
! *  Modified: 19-05-2009
! ***********************************

subroutine integratedWarpDeriv(ncoords,xyzface,indices_new)

  ! Take the current mesh and complete a warp based on the current
  ! Perturbed Coordinates

  use blockpointers
  use communication
  use mdDataLocal
  use mdData, only: mdNSurfNodesCompact
  use warpingPETSc
  implicit none
  !Subroutine Arguments
  integer(kind=intType)::ncoords
  real(kind=realType), dimension(3,ncoords)::xyzface
  integer(kind=intType),dimension(5,ncoords)::indices_new
  
  ! Local Arguments
  
  integer(kind=intType)::level=1,idxsurf,idxvol
  integer(kind=intType)::nn,sps=1,imax,jmax,kmax,i,j,k,n,mm,ll,nnn,ii,length
  real(kind=realType), dimension(:,:,:,:),allocatable::xyznew,xyz0,xyznewd
  integer(kind=intType),dimension(6)::IFACEPTB
  integer(kind=intType),dimension(12)::IEDGEPTB
  real(kind=realType),dimension(:),allocatable::xyzref

  real(kind=realType), dimension(2) :: time
  real(kind=realType)               :: timeAdjLocal, timeAdj
!  real(kind=realType), dimension(:,:,:,:,:,:,:),allocatable::xyzderiv
  real(kind=realType)::xref
  integer :: unitWarp = 8,ierror
  character(len = 16)::outfile,testfile
  
  write(testfile,100) myid!12
100 format (i5)  
  testfile=adjustl(testfile)
  write(outfile,101) trim(testfile)!testfile
101 format("ADWarpfile",a,".out")
  !outfile = "CSMachfile.txt"
  unitWarp = 8+myID
  !outfile = "CSMachfile.txt"
  
  open (UNIT=unitWarp,File=outfile,status='replace',action='write',iostat=ierror)
  if(ierror /= 0)                        &
       call terminate("integradtedWarpDeriv", &
       "Something wrong when &
       &calling open")
    

  print *,'in integrated warp AD'
  !Loop over each global surface node in turn
  do ii = 1,mdNSurfNodesCompact
     
     !loop over all three directions
     do ll = 1,3

        !set the face back to the current surface
        call updateFaces(ncoords,xyzFace,indices_new)

        !determine the length of the local compacted surface index array
        length = size(mdSurfGlobalIndLocal(1,:))

        !allocate temporary storage for the current values
        allocate(xyzref(length))
        
        do nnn=1,nDom
           !print*,'warpingblock',nnn
           call setPointersadj(nnn,level,sps)     
           !print *,'flag implicites'
           
           IMAX = IL
           JMAX = JL
           KMAX = KL
           ! LOOP THROUGH ALL local BLOCKS AND CALL WARPBLK WHERE APPROPRIATE
           ! SAVE NEW AND INITIAL XYZ VALUES TO BE PASSED TO WARPBLK
           
           !print *,'allocate xyz0'
           ALLOCATE(XYZ0(3,0:IMAX+1,0:JMAX+1,0:KMAX+1),XYZNEW(3,0:IMAX+1,0:JMAX+1,0:KMAX+1),XYZNEWd(3,0:IMAX+1,0:JMAX+1,0:KMAX+1))
           
           xyz0 = 0
           xyznew = 0
                 
           XYZ0(1,1:IMAX,1:JMAX,1:KMAX) = XInit(1:IMAX,1:JMAX,1:KMAX,1)
           XYZ0(2,1:IMAX,1:JMAX,1:KMAX) = XInit(1:IMAX,1:JMAX,1:KMAX,2)
           XYZ0(3,1:IMAX,1:JMAX,1:KMAX) = XInit(1:IMAX,1:JMAX,1:KMAX,3)
           XYZNEW(1,1:IMAX,1:JMAX,1:KMAX) = X(1:IMAX,1:JMAX,1:KMAX,1)
           XYZNEW(2,1:IMAX,1:JMAX,1:KMAX) = X(1:IMAX,1:JMAX,1:KMAX,2)
           XYZNEW(3,1:IMAX,1:JMAX,1:KMAX) = X(1:IMAX,1:JMAX,1:KMAX,3)
           !print *,'call warp local'
           
           !Zero the derivative seed
           xyznewd(:,:,:,:)  = 0.0
           
           !allocate temporary storage for the current values
           do j = 1,length
              if(mdSurfGlobalIndLocal(5,j)== ii-1)then
                 if(mdSurfGlobalIndLocal(4,j)== nnn)then
                 !update those values
                 xyznewd(ll,mdSurfGlobalIndLocal(1,j),mdSurfGlobalIndLocal(2,j),mdSurfGlobalIndLocal(3,j)) = 1.0
                 xyzref(j) = xyznew(ll,mdSurfGlobalIndLocal(1,j),mdSurfGlobalIndLocal(2,j),mdSurfGlobalIndLocal(3,j))
                    xyznew(ll,mdSurfGlobalIndLocal(1,j),mdSurfGlobalIndLocal(2,j),mdSurfGlobalIndLocal(3,j)) = xyzref(j)+ 1.0e-12
                 endif
              endif
           enddo
           
           !write(unitWarp,11)ifaceptb,iedgeptb
11         format(1x,'Faceb ',6I2,' edgeb ',12I2)
           call flagImplicitEdgesAndFacesDeriv(xyznewd,ifaceptb,iedgeptb)
           !write(unitWarp,12)ifaceptb,iedgeptb
           
           call WARP_LOCAL_D(xyznew, xyznewd, xyz0, ifaceptb, iedgeptb, imax&
                &  , jmax, kmax)

           !allocate temporary storage for the current values
           do j = 1,length
              if(mdSurfGlobalIndLocal(5,j)== ii-1)then
                 if(mdSurfGlobalIndLocal(4,j)== nnn)then
                    !replace perturbed values
                    xyznew(ll,mdSurfGlobalIndLocal(1,j),mdSurfGlobalIndLocal(2,j),mdSurfGlobalIndLocal(3,j)) = xyzref(j)
                 endif
              endif
           enddo
           
           !print *,'assign derivatives'
           call setpointersadj(nnn,level,sps) 
           
           IMAX = IL
           JMAX = JL
           KMAX = KL
           ! ASSIGN THESE derivative values
           DO I=1,IMAX
              DO J=1,JMAX
                 DO K=1,KMAX
                    do n = 1,3
                       idxvol = globalNode(i,j,k)*3+n
                       !idxsurf= (mdSurfIndLocal(5,mm)-1)*3+ll
                       idxsurf= (ii-1)*3+ll
                       !if(XYZNEWd(n,I,J,K).ne. zero)then
                       if(XYZNEWd(n,I,J,K)> 1e-10)then
                          ! write(unitWarp,12)ifaceptb,iedgeptb!'face',ifaceptb,'edge',iedgeptb
12                        format(1x,'Face ',6I2,' edge ',12I2)
                          write(unitWarp,13) idxsurf,idxvol,i,j,k,n,nnn,ii,ll,XYZNEWd(n,I,J,K)
13                        format(1x,'WarpSurf',9I8,f18.10)
                          !write(unitWarp,13) idxsurf,idxvol,mdSurfGlobalIndLocal(5,mm)*3+ll,i,j,k,n,nnn,nn,mm,ll,XYZNEWd(n,I,J,K)
!13                        format(1x,'WarpSurf',11I8,f18.10)
                       endif
                       !print *,'assigning',i,j,k
                       !xyzderiv(i,j,k,n,nn,mm,ll) =  XYZNEWd(n,I,J,K)
                       !if( abs(xyzderiv(i,j,k,n,nn,mm,ll))>1.0e-12)then
                       !   print *,'Derivatives',i,j,k,n,nn,mm,ll,xyzderiv(i,j,k,n,nn,mm,ll)
                       !endif
                       if (XYZNEWd(n,I,J,K).ne.0.0)then
                          call MatSetValues(dXvdXs, 1, idxvol-1, 1, idxsurf-1,   &
                               XYZNEWd(n,I,J,K), INSERT_VALUES, PETScIerr)
                          if( PETScIerr/=0 ) &
                               print *,'matrix setting error'!call errAssemb("MatSetValues", "verifydrdw")
                       endif
!!$                       X(I,J,K,1) = XYZNEW(1,I,J,K)
!!$                       X(I,J,K,2) = XYZNEW(2,I,J,K)
!!$                       X(I,J,K,3) = XYZNEW(3,I,J,K)
                    enddo
                 END DO
              END DO
           END DO
           
           
           deALLOCATE(XYZ0,XYZNEW,xyznewd)
        end do
        deallocate(xyzref)
     end do
     
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

      call MatAssemblyEnd(dXvdXs,MAT_FINAL_ASSEMBLY,PETScIerr)

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
        write(*,20) "Assembling dXv/dXs matrix time (s) =", timeAdj
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
         call MatView(dXvdXs,PETSC_VIEWER_DRAW_WORLD,PETScIerr)
         !call MatView(dXvdXsFD,PETSC_VIEWER_STDOUT_WORLD,PETScIerr)
         if( PETScIerr/=0 ) &
              call terminate("integratedWarpDerivParallel", "Error in MatView")
         !pause
      endif

      ! Flush the output buffer and synchronize the processors.

      call f77flush()
      call mpi_barrier(PETSC_COMM_WORLD, PETScIerr)
      ! Output formats.

   10 format(a,1x,i3,1x,a,1x,i3)
   20 format(a,1x,f8.2)
  print *,'warp derivatives finished'
end subroutine integratedWarpDeriv
   
