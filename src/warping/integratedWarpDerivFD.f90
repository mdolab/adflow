!
! ***********************************
! *  File: integratedWarpderivFD.f90
! *  Author: C.A.(Sandy) Mader
! *  Started: 19-05-2009
! *  Modified: 19-05-2009
! ***********************************

subroutine integratedWarpDerivFD(ncoords,xyzface,indices_new)

  ! Take the current mesh and complete a warp based on the current
  ! Perturbed Coordinates

  use blockpointers
  use communication
  implicit none
  !Subroutine Arguments
  integer(kind=intType)::ncoords
  real(kind=realType), dimension(3,ncoords)::xyzface
  integer(kind=intType),dimension(4,ncoords)::indices_new
  
  ! Local Arguments
  
  integer(kind=intType)::level=1
  integer(kind=intType)::nn,sps=1,imax,jmax,kmax,i,j,k,n,mm,ll,nnn
  real(kind=realType), dimension(:,:,:,:),allocatable::xyznew,xyz0,xyznewd
  integer(kind=intType),dimension(6)::IFACEPTB
  integer(kind=intType),dimension(12)::IEDGEPTB
!  real(kind=realType), dimension(:,:,:,:,:,:,:),allocatable::xyzderiv
 ! real(kind=realType), dimension(:,:,:,:,:,),allocatable::xyzderiv
  real(kind=realType)::xyzref,deltax = 1e-6,xderiv

  integer :: unitWarp = 8,ierr
  character(len = 16)::outfile,testfile
  
  write(testfile,100) myid!12
100 format (i5)  
  testfile=adjustl(testfile)
  write(outfile,101) trim(testfile)!testfile
101 format("FDWarpfile",a,".out")
  !outfile = "CSMachfile.txt"
  unitWarp = 8+myID
  !outfile = "CSMachfile.txt"
  
  open (UNIT=unitWarp,File=outfile,status='replace',action='write',iostat=ierr)
  if(ierr /= 0)                        &
       call terminate("integradtedWarpDeriv", &
       "Something wrong when &
       &calling open")
  !Allocate memory for the FD Blocks
  do nn = 1,nDom
     allocate(flowDoms(nn,level,sps)%xplus(3,0:Il+1,0:Jl+1,0:Kl+1), stat=ierr)
     if(ierr /= 0)                         &
          call terminate("integratedWarpDerivFD", &
          "Memory allocation failure for flowDoms%xPlus")
     allocate(flowDoms(nn,level,sps)%xminus(3,0:Il+1,0:Jl+1,0:Kl+1), stat=ierr)
     if(ierr /= 0)                         &
          call terminate("integratedWarpDerivFD", &
          "Memory allocation failure for flowDoms%xminus")
     !call setPointers(nn,level,sps)
     !allocate(xplus(3,0:Il+1,0:Jl+1,0:Kl+1))
     !allocate(xminus(3,0:Il+1,0:Jl+1,0:Kl+1))
     !print *,'xplus',xplus(1,1,1,1)
  enddo
  !loop over the surface points on this process to calculate the derivatives 
  !loop over domains
  do nn = 1,nDom
     call setPointers(nn,level,sps)
     print *,'flowdoms',flowdoms(nn,level,sps)%xplus(1,1,1,1)
     print *,'xplus1',xplus(1,1,1,1)

 !    print *,'derivative sizes',IMAX+1,JMAX+1,KMAX+1,3,ndom,ncoords,3
 !    allocate(xyzderiv(0:IMAX+1,0:JMAX+1,0:KMAX+1,3,ndom,ncoords,3))
     
     !loop over new coordinates array
     do mm = 1,ncoords
        !Check to see that coordinate is in this block. if so, update
        if(indices_new(4,mm)==nn)then
           print *,'getting sensitivites for surface point:',mm
           do ll = 1,3
              xyzref = xyzface(ll,mm)
              
              xyzface(ll,mm) = xyzref+deltax
              !reset the mesh coordinates to initial values
              do nnn=1,nDom
                 call setPointers(nnn,level,sps)
                 print *,'xplus2',xplus(1,1,1,1)
                 DO I=1,il!IMAX
                    DO J=1,jl!JMAX
                       DO K=1,kl!KMAX
                          X(I,J,K,1) = Xinit(I,J,K,1)
                          X(I,J,K,2) = Xinit(I,J,K,2)
                          X(I,J,K,3) = Xinit(I,J,K,3)
                       END DO
                    END DO
                 END DO
              enddo
              call setPointers(nn,level,sps)
              print *,'xplus3',xplus(1,1,1,1)
              !update the faces based on the new surface
              print *,'update faces'
              call updateFaces(ncoords,xyzFace,indices_new)
              
              print *,'exchangecoordinates'
              !Syncronize the faces to propogate the pertbation to adjacent blocks
              !  call exchangeCoor(level)
              call synchronizeBlockFaces(level,1)
              
              do nnn=1,nDom
                 print*,'warpingblock',nnn
                 call setPointers(nnn,level,sps)    
                 print *,'xplus4',xplus(1,1,1,1)
                 print *,'flag implicites'
                 !determine the explicitly and implicitly perturbed faces and edges
                 call flagImplicitEdgesAndFaces(ifaceptb,iedgeptb)
                 
                 IMAX = IL
                 JMAX = JL
                 KMAX = KL
                 ! LOOP THROUGH ALL local BLOCKS AND CALL WARPBLK WHERE APPROPRIATE
                 ! SAVE NEW AND INITIAL XYZ VALUES TO BE PASSED TO WARPBLK
                 
                 print *,'allocate xyz0'
                 ALLOCATE(XYZ0(3,0:IMAX+1,0:JMAX+1,0:KMAX+1),XYZNEW(3,0:IMAX+1,0:JMAX+1,0:KMAX+1))
                 
                 xyz0 = 0
                 xyznew = 0
                 
                 XYZ0(1,1:IMAX,1:JMAX,1:KMAX) = XInit(1:IMAX,1:JMAX,1:KMAX,1)
                 XYZ0(2,1:IMAX,1:JMAX,1:KMAX) = XInit(1:IMAX,1:JMAX,1:KMAX,2)
                 XYZ0(3,1:IMAX,1:JMAX,1:KMAX) = XInit(1:IMAX,1:JMAX,1:KMAX,3)
                 XYZNEW(1,1:IMAX,1:JMAX,1:KMAX) = X(1:IMAX,1:JMAX,1:KMAX,1)
                 XYZNEW(2,1:IMAX,1:JMAX,1:KMAX) = X(1:IMAX,1:JMAX,1:KMAX,2)
                 XYZNEW(3,1:IMAX,1:JMAX,1:KMAX) = X(1:IMAX,1:JMAX,1:KMAX,3)
                 print *,'call warp local'
                 
                 print*,'warpingblock',nn,ll,indices_new(1,mm),indices_new(2,mm),indices_new(3,mm)
                 call WARP_LOCAl(xyznew, xyz0, ifaceptb, iedgeptb, imax&
                      &  , jmax, kmax)
                 
                 call setPointers(nnn,level,sps) 
                 print *,'xplus5',xplus(1,1,1,1)
                 xplus(:,:,:,:) = xyznew
                 
                           
                 deALLOCATE(XYZ0,XYZNEW)
              end do
              
              xyzface(ll,mm)=xyzref-deltax

              do nnn=1,nDom
                 call setPointers(nnn,level,sps)
                
                 DO I=1,il!IMAX
                    DO J=1,jl!JMAX
                       DO K=1,kl!KMAX
                          X(I,J,K,1) = Xinit(I,J,K,1)
                          X(I,J,K,2) = Xinit(I,J,K,2)
                          X(I,J,K,3) = Xinit(I,J,K,3)
                       END DO
                    END DO
                 END DO
              enddo
              call setPointers(nn,level,sps)
              
              !update the faces based on the new surface
              print *,'update faces'
              call updateFaces(ncoords,xyzFace,indices_new)
              
              print *,'exchangecoordinates'
              !Syncronize the faces to propogate the pertbation to adjacent blocks
              !  call exchangeCoor(level)
              call synchronizeBlockFaces(level,1)
              
              do nnn=1,nDom
                 print*,'warpingblock',nnn
                 call setPointers(nnn,level,sps)     
                 print *,'flag implicites'
                 !determine the explicitly and implicitly perturbed faces and edges
                 call flagImplicitEdgesAndFaces(ifaceptb,iedgeptb)
                 
                 IMAX = IL
                 JMAX = JL
                 KMAX = KL
                 ! LOOP THROUGH ALL local BLOCKS AND CALL WARPBLK WHERE APPROPRIATE
                 ! SAVE NEW AND INITIAL XYZ VALUES TO BE PASSED TO WARPBLK
                 
                 print *,'allocate xyz0'
                 ALLOCATE(XYZ0(3,0:IMAX+1,0:JMAX+1,0:KMAX+1),XYZNEW(3,0:IMAX+1,0:JMAX+1,0:KMAX+1))
                 
                 xyz0 = 0
                 xyznew = 0
                 
                 XYZ0(1,1:IMAX,1:JMAX,1:KMAX) = XInit(1:IMAX,1:JMAX,1:KMAX,1)
                 XYZ0(2,1:IMAX,1:JMAX,1:KMAX) = XInit(1:IMAX,1:JMAX,1:KMAX,2)
                 XYZ0(3,1:IMAX,1:JMAX,1:KMAX) = XInit(1:IMAX,1:JMAX,1:KMAX,3)
                 XYZNEW(1,1:IMAX,1:JMAX,1:KMAX) = X(1:IMAX,1:JMAX,1:KMAX,1)
                 XYZNEW(2,1:IMAX,1:JMAX,1:KMAX) = X(1:IMAX,1:JMAX,1:KMAX,2)
                 XYZNEW(3,1:IMAX,1:JMAX,1:KMAX) = X(1:IMAX,1:JMAX,1:KMAX,3)
                 print *,'call warp local'
                 
                 print*,'warpingblock',nn,ll,indices_new(1,mm),indices_new(2,mm),indices_new(3,mm)
                 call WARP_LOCAl(xyznew, xyz0, ifaceptb, iedgeptb, imax&
                      &  , jmax, kmax)

                 xminus = xyznew
                 
                 deALLOCATE(XYZ0,XYZNEW)

              end do
              xyzface(ll,mm)=xyzref
         
              print *,'assign derivatives'
              ! ASSIGN THESE derivative values
              do nnn=1,nDom
                 print*,'warpingblock',nnn,level,sps
                 call setPointers(nnn,level,sps)     
                 IMAX = IL
                 JMAX = JL
                 KMAX = KL
                 DO I=1,IMAX
                    DO J=1,JMAX
                       DO K=1,KMAX
                          do n = 1,3
                             !print *,'indices',i,j,k,n,nnn
                             xderiv = (xplus(n,i,j,k)-xminus(n,i,j,k))/(2*deltax)
                             write(unitWarp,13) xderiv,i,j,k,n,nnn,nn,mm,ll
13                           format(1x,'WarpSurf',f18.10,8I8)
                             !print *,'assigning',i,j,k
                             
                          enddo
                       END DO
                    END DO
                 END DO
              end do
           end do
        endif
     end do
     !deallocate(xyzderiv)
  end do
   !Allocate memory for the FD Blocks
  do nn = 1,nDom
     call setPointers(nn,level,sps)
     deallocate(xplus)
     deallocate(xminus)
  enddo
  print *,'warp derivatives finished'
end subroutine integratedWarpDerivFD
   
