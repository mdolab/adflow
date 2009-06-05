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
  implicit none
  !Subroutine Arguments
  integer(kind=intType)::ncoords
  real(kind=realType), dimension(3,ncoords)::xyzface
  integer(kind=intType),dimension(5,ncoords)::indices_new
  
  ! Local Arguments
  
  integer(kind=intType)::level=1,idxsurf,idxvol
  integer(kind=intType)::nn,sps=1,imax,jmax,kmax,i,j,k,n,mm,ll,nnn
  real(kind=realType), dimension(:,:,:,:),allocatable::xyznew,xyz0,xyznewd
  integer(kind=intType),dimension(6)::IFACEPTB
  integer(kind=intType),dimension(12)::IEDGEPTB
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
    
  !loop over the surface points on this process to calculate the derivatives 
  !loop over domains
  do nn = 1,nDom
     call setPointersadj(nn,1,sps)
   

 !    print *,'derivative sizes',IMAX+1,JMAX+1,KMAX+1,3,ndom,ncoords,3
 !    allocate(xyzderiv(0:IMAX+1,0:JMAX+1,0:KMAX+1,3,ndom,ncoords,3))
     
     !loop over new coordinates array
     do mm = 1,2!ncoords
        !Check to see that coordinate is in this block. if so, update
        if(indices_new(4,mm)==nn)then
           print *,'getting sensitivites for surface point:',mm
           do ll = 1,3
!!$              !reset the mesh coordinates to initial values
!!$              do nnn=1,nDom
!!$                 call setPointersadj(nnn,level,sps)
!!$                 DO I=1,il!IMAX
!!$                    DO J=1,jl!JMAX
!!$                       DO K=1,kl!KMAX
!!$                          X(I,J,K,1) = Xinit(I,J,K,1)
!!$                          X(I,J,K,2) = Xinit(I,J,K,2)
!!$                          X(I,J,K,3) = Xinit(I,J,K,3)
!!$                       END DO
!!$                    END DO
!!$                 END DO
!!$              enddo
!!$              call setPointersadj(nn,level,sps)
              
              !update the faces based on the new surface
              print *,'update faces'
              call updateFaces(ncoords,xyzFace,indices_new)
              
              print *,'exchangecoordinates'
              !Syncronize the faces to propogate the pertbation to adjacent blocks
              !  call exchangeCoor(level)
              call synchronizeBlockFaces(level,1)
              
              do nnn=1,nDom
                 print*,'warpingblock',nnn
                 call setPointersadj(nnn,level,sps)     
                 print *,'flag implicites'
                 
                 IMAX = IL
                 JMAX = JL
                 KMAX = KL
                 ! LOOP THROUGH ALL local BLOCKS AND CALL WARPBLK WHERE APPROPRIATE
                 ! SAVE NEW AND INITIAL XYZ VALUES TO BE PASSED TO WARPBLK
                 
                 print *,'allocate xyz0'
                 ALLOCATE(XYZ0(3,0:IMAX+1,0:JMAX+1,0:KMAX+1),XYZNEW(3,0:IMAX+1,0:JMAX+1,0:KMAX+1),XYZNEWd(3,0:IMAX+1,0:JMAX+1,0:KMAX+1))
                 
                 
                 
                 
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
                 !if(indices_new(4,mm)==nnn)then
                 xyznewd(ll,indices_new(1,mm),indices_new(2,mm),indices_new(3,mm)) = 1.0
                 xref = xyznew(ll,indices_new(1,mm),indices_new(2,mm),indices_new(3,mm))
                 xyznew(ll,indices_new(1,mm),indices_new(2,mm),indices_new(3,mm)) = xref+ 1.0e-12
                 !endif
                 !determine the explicitly and implicitly perturbed faces and edges
                 !call flagImplicitEdgesAndFaces(ifaceptb,iedgeptb)
                 !write(unitWarp,11)ifaceptb,iedgeptb
11               format(1x,'Faceb ',6I2,' edgeb ',12I2)
                 call flagImplicitEdgesAndFacesDeriv(xyznewd,ifaceptb,iedgeptb)
                 write(unitWarp,12)ifaceptb,iedgeptb
                 
                 print*,'warpingblock',nn,ll,indices_new(1,mm),indices_new(2,mm),indices_new(3,mm)
                 call WARP_LOCAL_D(xyznew, xyznewd, xyz0, ifaceptb, iedgeptb, imax&
                      &  , jmax, kmax)
                 xyznew(ll,indices_new(1,mm),indices_new(2,mm),indices_new(3,mm)) = xref
                 print *,'assign derivatives'
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
                             idxsurf= (mdSurfIndLocal(5,mm)-1)*3+ll
                             !if(XYZNEWd(n,I,J,K).ne. zero)then
                             if(XYZNEWd(n,I,J,K)> 1e-10)then
                               ! write(unitWarp,12)ifaceptb,iedgeptb!'face',ifaceptb,'edge',iedgeptb
12                              format(1x,'Face ',6I2,' edge ',12I2)
                                write(unitWarp,13) idxsurf,idxvol,mdSurfGlobalIndLocal(5,mm)*3+ll,i,j,k,n,nnn,nn,mm,ll,XYZNEWd(n,I,J,K)
13                              format(1x,'WarpSurf',11I8,f18.10)
                             endif
                             !print *,'assigning',i,j,k
                            !xyzderiv(i,j,k,n,nn,mm,ll) =  XYZNEWd(n,I,J,K)
                             !if( abs(xyzderiv(i,j,k,n,nn,mm,ll))>1.0e-12)then
                             !   print *,'Derivatives',i,j,k,n,nn,mm,ll,xyzderiv(i,j,k,n,nn,mm,ll)
                             !endif
!!$                       X(I,J,K,1) = XYZNEW(1,I,J,K)
!!$                       X(I,J,K,2) = XYZNEW(2,I,J,K)
!!$                       X(I,J,K,3) = XYZNEW(3,I,J,K)
                          enddo
                       END DO
                    END DO
                 END DO
                 
                 
                 deALLOCATE(XYZ0,XYZNEW,xyznewd)
              end do
           end do
        endif
     end do
     !deallocate(xyzderiv)
  end do
  print *,'warp derivatives finished'
end subroutine integratedWarpDeriv
   
