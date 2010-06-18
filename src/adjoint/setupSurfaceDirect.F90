!
!     ******************************************************************
!     *                                                                *
!     * File:          setupSurfaceDirect.F90                          *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 02-03-2010                                      *
!     * Last modified: 02-03-2010                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine setupSurfaceDirect
!
!     ******************************************************************
!     *                                                                *
!     * Setup up the direct method for a surface node                  *
!     *                                                                *
!     ******************************************************************
!
#ifndef USE_NO_PETSC
      use ADjointPETSc
      use ADjointVars
      !use warpingPETSc
      use WarpingPETSc, only: dXvdXsDV,drdxsDV,didxs2,djdxs2
      use mdData, only: mdNSurfNodesCompact
      use blockpointers !globalnode
      use inputTimeSpectral !nTimeIntervalsSpectral
      use flowvarrefstate !nw
      use communication   !sumb_comm_world
      implicit none

!
!     Subroutine arguments.
!
!
!     Local variables.
!
      real(kind=realType), dimension(2) :: time
      real(kind=realType)               :: timeAdjLocal, timeAdj

!      integer(kind=intType) :: idx!tmp for fd

      integer(kind=inttype)::i,j,k,n,idxres,idxnode,sps,sps2,nnn,m,nn
      character(len=2*maxStringLen) :: errorMessage

      character fileName*32, dataName*32
      !real(kind=realType) :: timeAdj, timeFD, timeResAdj,
      real(kind=realType)::values
      integer :: ierr, testnode
 
      integer :: unitdRdxs = 8,ierror,unitpvr
      character(len = 25)::outfile,testfile
      write(testfile,100) myid!12
100   format (i5)  
      testfile=adjustl(testfile)
      write(outfile,101) trim(testfile)!testfile
101   format("ADdRdxsfiledir",a,".out")
      unitdrdxs = 8+myID

      
      open (UNIT=unitdRdxs,File=outfile,status='replace',action='write',iostat=ierror)
      if(ierror /= 0)                        &
           call terminate("verifydRdxFile", &
           "Something wrong when &
           &calling open")

      write(testfile,102) myid!12
102   format (i5)  
      testfile=adjustl(testfile)
      write(outfile,103) trim(testfile)!testfile
103   format("ADpvrfiledir",a,".out")
      unitpvr = 18+myID

      
      open (UNIT=unitpvr,File=outfile,status='replace',action='write',iostat=ierror)
      if(ierror /= 0)                        &
           call terminate("setupsurfacedirect", &
           "Something wrong when &
           &calling open")

 
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!

!
      ! Send some feedback to screen.

      if( PETScRank==0 ) &
        write(*,10) "Setting up direct surface solution..."

      !create the PETSc Vector aswell
      ! Create the vector. Depending on either this is a sequential or 
      ! parallel run,  PETSc automatically generates the apropriate
      ! vector type over all processes in SUMB_PETSC_COMM_WORLD.

      call VecCreate(SUMB_PETSC_COMM_WORLD, selector, PETScIerr)


      if( PETScIerr/=0 ) &
        call terminate("setupSurfaceDirect", "Error in VecCreate selector")
      
      ! Set the local size and let PETSc determine its global size

      call VecSetSizes(selector,PETSC_DECIDE,3*mdNSurfNodesCompact,PETScIerr)
     
      if( PETScIerr/=0 ) then
        write(errorMessage,99) &
              "Error in VecSetSizes selctor for global size", mdNSurfNodesCompact
        call terminate("setupSurfaceDirect", errorMessage)
      endif

      ! Set the vector from options.

      call VecSetFromOptions(selector, PETScIerr)

      if( PETScIerr/=0 ) &
        call terminate("setupSurfaceDirect", &
                       "Error in VecSetFromOptions selector")

      ! Get the initial time.

      call cpu_time(time(1))

      !set column of interest
      call VecSetValue(selector, 0, 1.0 ,INSERT_VALUES, PETScIerr)

      call VecAssemblyBegin(selector,PETScIerr)
       
      if( PETScIerr/=0 ) &
           call terminate("setupASjointRHS", "Error in VecAssemblyBegin")  
      
      call VecAssemblyEnd  (selector,PETScIerr)
       
      if( PETScIerr/=0 ) &
           call terminate("setupADjointRHS", "Error in VecAssemblyEnd")
       
      if( debug ) then
         call VecView(selector,PETSC_VIEWER_DRAW_WORLD,PETScIerr)
         !call VecView(pvr,PETSC_VIEWER_STDOUT_WORLD,PETScIerr)
         if( PETScIerr/=0 ) &
              call terminate("setupADjointRHS", "Error in VecView")
         pause
      endif

      !Alternate solution      
      !multiply and store in drdxs
      call MatMatMult(dRdx,dXvdXsDV,MAT_INITIAL_MATRIX,PETSC_DEFAULT_DOUBLE_PRECISION,dRdXsDV, PETScIerr) 
      sps = 1
      print *,'sps',sps
      do nn = 1,mdNSurfNodesCompact
         print *,'nn',nn
         do m = 1, 3
            idxnode   = (nn-1)*3+m 
            do sps2 = 1,nTimeIntervalsSpectral
               do nnn = 1,ndom
                  call setPointersAdj(nnn,1,sps2)
                  DO I=2,Il
                     DO J=2,Jl
                        DO K=2,Kl
                           do n = 1,nw
                              idxres = globalCell(i,j,k)*nw+n
                              if ((idxres-1)>=0 .and. (idxnode-1)>=0)then
                                 call MatGetValues(drdxsdv,1,idxres-1,1,idxnode-1,values,PETScIerr)
                                 !print *,'value',values
                                 !if(value.ne.0)then
                                 if(abs(values)>1e-10)then
                                    !write(unitWarp,12)ifaceptb,iedgeptb !'face',ifaceptb,'edge',iedgeptb
                                    !12                                     format(1x,'Face',6I2,'edge',12I2)
                                    write(unitdrdxs,13) idxnode,idxres,m,nn,n,k,j,i,nnn,sps2,values
                                    !write(unitWarp,13) xderiv,i,j,k,n,nnn,nn,mm,ll
13                                  format(1x,'drdxs',10I8,f18.10)
                                 endif
                              end if
                           enddo
                        END DO
                     END DO
                  END DO
                  !print *,'spsloop',sps
                  call setPointersAdj(nnn,1,sps)
               end do
            enddo
         end do
      enddo
      
      call mpi_barrier(SUmb_comm_world, ierr)
      
      close(unitdrdxs)
      
      call mpi_barrier(SUmb_comm_world, ierr)

      call MatMult(dRdXsDV,selector,pvr,PETScIerr)
  
      if( debug ) then
         !call VecView(selector,PETSC_VIEWER_DRAW_WORLD,PETScIerr)
         call VecView(pvr,PETSC_VIEWER_STDOUT_WORLD,PETScIerr)
         if( PETScIerr/=0 ) &
              call terminate("setupADjointRHS", "Error in VecView")
         !pause
      endif
      
      print *,'sps',sps
      do nn = 1,1!mdNSurfNodesCompact
         print *,'nn',nn
         do m = 1,1! 3
            idxnode   = (nn-1)*3+m 
            do sps2 = 1,nTimeIntervalsSpectral
               do nnn = 1,ndom
                  call setPointersAdj(nnn,1,sps2)
                  DO I=2,Il
                     DO J=2,Jl
                        DO K=2,Kl
                           do n = 1,nw
                              idxres = globalCell(i,j,k)*nw+n
                              if ((idxres-1)>=0 .and. (idxnode-1)>=0)then
                                 call VecGetValues(pvr,1,idxres-1,values,PETScIerr)
                                 !print *,'value',values
                                 !if(value.ne.0)then
                                 if(abs(values)>1e-10)then
                                    !write(unitWarp,12)ifaceptb,iedgeptb !'face',ifaceptb,'edge',iedgeptb
                                    !12                                     format(1x,'Face',6I2,'edge',12I2)
                                    write(unitpvr,14) idxnode,idxres,m,nn,n,k,j,i,nnn,sps2,values
                                    !write(unitWarp,13) xderiv,i,j,k,n,nnn,nn,mm,ll
14                                  format(1x,'pvr',10I8,f18.10)
                                 endif
                              end if
                           enddo
                        END DO
                     END DO
                  END DO
                  !print *,'spsloop',sps
                  call setPointersAdj(nnn,1,sps)
               end do
            enddo
         end do
      enddo

      call mpi_barrier(SUmb_comm_world, ierr)
      close(unitpvr)

      call f77flush()
      call mpi_barrier(SUMB_PETSC_COMM_WORLD, PETScIerr)

      ! Output format.

   10 format(a)
   20 format(a,1x,f8.2)
   99 format(a,1x,i6)
!
#endif

    end subroutine setupSurfaceDirect
