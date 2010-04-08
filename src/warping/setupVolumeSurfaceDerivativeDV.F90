!
! ***********************************
! *  File: setupVolumeSurfaceDerivativesDV.F90
! *  Author: C.A.(Sandy) Mader
! *  Started: 15-07-2009
! *  Modified: 25-01-2010
! ***********************************

subroutine setupVolumeSurfaceDerivativesDV

  ! Find the derivatives of the mesh wrt the surface at the current point

  use blockpointers
  use communication
  use mdDataLocal
  use mdData, only: mdNSurfNodesCompact,mdGlobalSurfxx
  use warpingPETSc
  use inputTimeSpectral !nTimeIntervalsSpectral
  !use ADjointPETSc, only: PETScOne, value
  use section
  use monitor
  implicit none
  !Subroutine Arguments

  
  ! Local Arguments
  
  integer(kind=intType)::level=1,idxvol,idxsurf
  integer(kind=intType)::nn,sps=1,imax,jmax,kmax,i,j,k,n,mm,ll,nnn
  real(kind=realType), dimension(:,:,:,:),allocatable::xyznew,xyz0,xyznewd
  integer(kind=intType),dimension(6)::IFACEPTB
  integer(kind=intType),dimension(12)::IEDGEPTB

  real(kind=realType), dimension(2) :: time
  real(kind=realType)               :: timeAdjLocal, timeAdj
  real(kind=realType)::xref 

  real(kind=realType), dimension(nSections) :: dt, t
  real(kind=realType) :: tNew, tOld
  real(kind=realType) :: xp, yp, zp
  
  real(kind=realType), dimension(3)   :: rotationPoint
  real(kind=realType), dimension(3,3) :: rotationMatrix

  real(kind=realType), allocatable,dimension(:,:)::surfaceNodes
 
#ifndef USE_NO_PETSC  
  
  ! Send some feedback to screen.

  if( PETScRank==0 ) &
       write(*,10) "Assembling dXv/dXsDV Parallel matrix..."

!!$  ! ensure that the current derivatives are up to date
  allocate(surfaceNodes(3,mdNSurfNodesCompact))
  surfaceNodes = mdGlobalSurfxx(:,:,1)
  !print *,'surface Nodes',surfaceNodes
  call updateFacesGlobal(mdNSurfNodesCompact,surfaceNodes)
  
  ! Get the initial time.
  
  call cpu_time(time(1))
  
  !!zero the matrix for dXvdXsPara ADD call
  !call MatZeroEntries(dXvdXsPara,PETScIerr)
  !zero the matrix for dXvdXs ADD call
  call MatZeroEntries(dXvdXsDV,PETScIerr)
  
  if( PETScIerr/=0 ) &
       call terminate("setupVolumeSurfaceDerivativesDV", "Error in MatZeroEntries dXvdXsDV")
  
  !loop over the Global surface points on this process to calculate the derivatives 
  !if (myID==0) print *,'global surface',mdGlobalSurfxx(:,:,1)
  !Set some time variables for the rotations
  timeUnsteady = zero

  !loop over time instances
  do sps = 1,nTimeIntervalsSpectral
     
     !compute the rotation matrix for this time instance

     ! Compute the corresponding times for this spectral solution
     ! and call updateCoorFineMesh to determine the coordinates.
     
     do nn=1,nSections
        t(nn) = (sps-1)*dt(nn)
     enddo
     
     ! Determine the time values of the old and new time level.
     ! It is assumed that the rigid body rotation of the mesh is only
     ! used when only 1 section is present.
     
     tNew = timeUnsteady + timeUnsteadyRestart
     tOld = tNew - t(1)
     
     ! Compute the rotation matrix of the rigid body rotation as
     ! wellas the rotation point; the latter may vary in time due
     ! to rigid body translation.
     
     call rotMatrixRigidBody(tNew, tOld, rotationMatrix, rotationPoint)

     !loop over domains
     do nn = 1,nDom
        call setPointersAdj(nn,1,sps)
      
        !loop over new coordinates array
        !print *,'number of local surface nodes',myID,mdNGlobalSurfNodesLocal(myID+1)
        do mm = 1,mdNGlobalSurfNodesLocal(myID+1)
           !Check to see that coordinate is in this block. if so, update
           !print *,'mdSurfGlobalIndLoca',mdSurfGlobalIndLocal(4,mm),nn
           if( mdSurfGlobalIndLocal(4,mm)==nn)then
              
              !only local block needs to be perturbed. Index sychronization will take care of the rest
              
              IMAX = IL
              JMAX = JL
              KMAX = KL
              
              ! SAVE NEW AND INITIAL XYZ VALUES TO BE PASSED TO WARPBLK
              
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
                 
                 !zero the AD perturbation
                 xyznewd(:,:,:,:)  = 0.0
                 
                 !set the AD perturbation. A small real perturbation is set as 
                 !well to keep out of the funny region for deltax<eps in the
                 !warping algorithm
                 xyznewd(ll, mdSurfGlobalIndLocal(1,mm), mdSurfGlobalIndLocal(2,mm), mdSurfGlobalIndLocal(3,mm)) = 1.0
                 xref = xyznew(ll, mdSurfGlobalIndLocal(1,mm), mdSurfGlobalIndLocal(2,mm), mdSurfGlobalIndLocal(3,mm))
                 xyznew(ll, mdSurfGlobalIndLocal(1,mm), mdSurfGlobalIndLocal(2,mm), mdSurfGlobalIndLocal(3,mm)) = xref+ 1.0e-12
                 
                 !determine the explicitly and implicitly perturbed
                 !faces and edges
                 call flagImplicitEdgesAndFacesDeriv(xyznewd,ifaceptb,iedgeptb)
                 !print *,'warpl_local_d',sps,nn,ifaceptb,iedgeptb
                 !Warp the block
                 call WARP_LOCAL_D(xyznew, xyznewd, xyz0, ifaceptb, iedgeptb, imax&
                      &  , jmax, kmax)
                 
                 !reset the small real perturbation
                 xyznew(ll,mdSurfGlobalIndLocal(1,mm),mdSurfGlobalIndLocal(2,mm),mdSurfGlobalIndLocal(3,mm)) = xref
                 
                 ! ASSIGN THESE derivative values
                 DO I=1,IMAX
                    DO J=1,JMAX
                       DO K=1,KMAX
                          !add the derivative rotation in
                          xp = xyznewd(1,i,j,k)
                          yp = xyznewd(2,i,j,k)
                          zp = xyznewd(3,i,j,k) 
                          
                          xyznewd(1,i,j,k)= rotationMatrix(1,1)*xp &
                               + rotationMatrix(1,2)*yp &
                               + rotationMatrix(1,3)*zp 
                          xyznewd(2,i,j,k)= rotationMatrix(2,1)*xp &
                               + rotationMatrix(2,2)*yp &
                               + rotationMatrix(2,3)*zp 
                          xyznewd(3,i,j,k)= rotationMatrix(3,1)*xp &
                               + rotationMatrix(3,2)*yp &
                               + rotationMatrix(3,3)*zp

                          do n = 1,3
                             !no need to alter sps here since in this case
                             !all time instances are independent, the overall
                             !sps loop takes care of this...
                             idxvol = globalNode(i,j,k)*3+n
                             !However, we need to add sps variation here
                             !idxsurf= mdNSurfNodesCompact*(sps-1)+&
                             !     mdSurfGlobalIndLocal(5,mm)*3+ll!1
                             !Sps variation now included above
                             !idxsurf= (mdSurfGlobalIndLocal(5,mm)-1)*3+ll!1
                             idxsurf= (mdSurfGlobalIndLocal(5,mm))*3+ll!1
                             
                             if (xyznewd(n,I,J,K).ne.0.0)then
                                
!!$                             call MatSetValues(dXvdXsPara, 1, idxvol-1, 1, idxsurf-1,   &
!!$                                               xyznewd(n,I,J,K), ADD_VALUES, PETScIerr)
                                call MatSetValues(dXvdXsDV, 1, idxvol-1, 1, idxsurf-1,   &
                                     xyznewd(n,I,J,K), ADD_VALUES, PETScIerr)
                              !  call MatSetValues(dXvdXsDV, 1, 0, 1, idxsurf-1,   &
                              !       mdGlobalSurfxx(ll,mdSurfGlobalIndLocal(5,mm)+1,sps), INSERT_VALUES, PETScIerr)
                                
                                if( PETScIerr/=0 ) &
                                     print *,'matrix setting error'

                             endif

                          enddo
                       END DO
                    END DO
                 END DO
                 !print *,'idxsurf',idxsurf,mdSurfGlobalIndLocal(5,mm),ll
              end do
              deALLOCATE(XYZ0,XYZNEW,xyznewd)
           endif
        end do
     end do
  enddo
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

      !call MatAssemblyBegin(dXvdXsPara,MAT_FINAL_ASSEMBLY,PETScIerr)
      call MatAssemblyBegin(dXvdXsDV,MAT_FINAL_ASSEMBLY,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("setupVolumeSurfaceDerivativesDV", &
                       "Error in MatAssemblyBegin dXvdXsDV")

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

      !call MatAssemblyEnd  (dXvdXsPara,MAT_FINAL_ASSEMBLY,PETScIerr)
      call MatAssemblyEnd  (dXvdXsDV,MAT_FINAL_ASSEMBLY,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("setupVolumeSurfaceDerivativesDV", &
                       "Error in MatAssemblyEnd dXvdXsDV")

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

      !call MatSetOption(dXvdXsPara,MAT_NO_NEW_NONZERO_LOCATIONS,PETScIerr)
      !call MatSetOption(dXvdXsDV,MAT_NO_NEW_NONZERO_LOCATIONS,PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("setupVolumeSurfaceDerivativesDV", &
                       "Error in MatSetOption dXvdXsDV")

      ! Get new time and compute the elapsed time.

      call cpu_time(time(2))
      timeAdjLocal = time(2)-time(1)

      ! Determine the maximum time using MPI reduce
      ! with operation mpi_max.

      call mpi_reduce(timeAdjLocal, timeAdj, 1, sumb_real, &
                      mpi_max, 0, PETSC_COMM_WORLD, PETScIerr)

      if( PETScRank==0 ) &
        write(*,20) "Assembling dXv/dXsDV Parallel matrix time (s) =", timeAdj
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
         !call MatView(dXvdXsPara,PETSC_VIEWER_DRAW_WORLD,PETScIerr)
         !call MatView(dXvdXsPara,PETSC_VIEWER_STDOUT_WORLD,PETScIerr)
         !all MatView(dXvdXsDV,PETSC_VIEWER_DRAW_WORLD,PETScIerr)
         call MatView(dXvdXsDV,PETSC_VIEWER_STDOUT_WORLD,PETScIerr)
         !call MatView(dXvdXsDV,PETSC_VIEWER_STDOUT_SELF,PETScIerr)
         if( PETScIerr/=0 ) &
          call terminate("setupVolumeSurfaceDerivativesDV", "Error in MatViewDV")
        pause
      endif

      ! Flush the output buffer and synchronize the processors.
      deallocate(surfaceNodes)
      call f77flush()
      call mpi_barrier(PETSC_COMM_WORLD, PETScIerr)

      ! Output formats.

   10 format(a,1x,i3,1x,a,1x,i3)
   20 format(a,1x,f8.2)

#endif
    end subroutine setupVolumeSurfaceDerivativesDV

