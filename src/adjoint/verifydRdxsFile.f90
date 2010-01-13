!
!     ******************************************************************
!     *                                                                *
!     * File:          verifydRdxsfile.f90                              *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 07-22-2009                                      *
!     * Last modified: 07-22-2009                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine verifydRdxsfile
!
!     ******************************************************************
!     *                                                                *
!     *  This subroutine computes the values for dRdw and compares     *
!     *  the Tapenade result to the finite-difference results.         *
!     *                                                                *
!     ******************************************************************


      use blockPointers ! block (nDoms,flowDoms), globalCell
      use flowvarrefstate
      use communication
      use iteration     ! groundLevel
      use inputTimeSpectral ! spaceDiscr,nTimeIntervalsSpectral
      use inputIO
      use cgnsGrid        ! cgnsNFamilies

      !from old verify routine
      use ADjointPETSc, only: drdx,petscone,insert_values,petscierr,mat_final_assembly,petsc_viewer_draw_world,petsc_viewer_stdout_world,add_values,mat_initial_matrix,PETSC_DEFAULT_DOUBLE_PRECISION
      use warpingPETSc, only: dxvdxs,dRdXs!,drdxs
      !use FDPETSc, only: DRDWFD
      use precision
      !use blockPointers
      !use flowvarrefstate
      !use iteration
      use inputIteration
      use mdData
      !implicit none
     

      implicit none

!
!     Subroutine arguments
      integer(kind=intType) :: level = 1


!
!     Local variables 
!
      integer(kind=intType) :: i, j, k, n,nn,nnn,sps = 1
      integer(kind=intType) ::m,idxres,idxnode,famid,sps2
      real(kind=realType), dimension(10) :: time
      real(kind=realType) :: timeRes

  
      character fileName*32, dataName*32
      real(kind=realType) :: timeAdj, timeFD, timeResAdj,value
      integer :: ierr, testnode
 
      integer :: unitdRdxs = 8,ierror
      character(len = 20)::outfile,testfile
      write(testfile,100) myid!12
100   format (i5)  
      testfile=adjustl(testfile)
      write(outfile,101) trim(testfile)!testfile
101   format("ADdRdxsfile",a,".out")
      unitdrdxs = 8+myID

      
      open (UNIT=unitdRdxs,File=outfile,status='replace',action='write',iostat=ierror)
      if(ierror /= 0)                        &
           call terminate("verifydRdxFile", &
           "Something wrong when &
           &calling open")

      !setup drdx
      call setupGradientMatrixSpatial(level)
      
      !****************************
      !Now set up dxvdxs
      !*********************
       if(cgnsNfamilies > 0) then
         famID = 1
      else
         famID = 0
      endif
      call initializeWarping(famID)
      
      call setupVolumeSurfaceDerivatives

      !multiply and store in drdxs
      call MatMatMult(dRdx,dXvdXs,MAT_INITIAL_MATRIX,PETSC_DEFAULT_DOUBLE_PRECISION,dRdXs, PETScIerr) 


      !now extract and write to a file
      do sps = 1,nTimeIntervalsSpectral
         do nn = 1,mdNSurfNodesCompact
            do m = 1, 3
               idxnode   = (nn-1)*3+m+mdNSurfNodesCompact*(sps-1) 
               do sps2 = 1,nTimeIntervalsSpectral
                  do nnn = 1,ndom
                     call setPointersAdj(nnn,1,sps2)
                     DO I=2,Il
                        DO J=2,Jl
                           DO K=2,Kl
                              do n = 1,nw
                                 idxres = globalCell(i,j,k)*nw+n
                                 if ((idxres-1)>=0 .and. (idxnode-1)>=0)then
                                    call MatGetValues(drdxs,1,idxres-1,1,idxnode-1,value,PETScIerr)
                                    !if(value.ne.0)then
                                    if(abs(value)>1e-10)then
                                       !write(unitWarp,12)ifaceptb,iedgeptb !'face',ifaceptb,'edge',iedgeptb
                                       !12                                     format(1x,'Face',6I2,'edge',12I2)
                                       write(unitdrdxs,13) idxnode,idxres,m,nn,sps,n,k,j,i,nnn,sps2,value
                                       !write(unitWarp,13) xderiv,i,j,k,n,nnn,nn,mm,ll
13                                     format(1x,'drdxs',9I8,f18.10)
                                    endif
                                 end if
                              enddo
                           END DO
                        END DO
                     END DO
                     call setPointersAdj(nnn,1,sps)
                  end do
               enddo
            end do
         enddo
      end do
      call mpi_barrier(SUmb_comm_world, ierr)
      
      close(unitdrdxs)
      
      call mpi_barrier(SUmb_comm_world, ierr)
      
      
111   format(4I4, 3ES22.14)
      
      
      
    end subroutine verifydRdxsfile
