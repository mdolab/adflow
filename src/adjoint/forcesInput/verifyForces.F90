subroutine verifyForces(pts,npts,nTS)

  ! This routine does three thing:
  ! 1. Check that getForces and computeForcesAdj give the same results
  ! 2. Check that dFdw matrix is correct
  ! 3. Check that dFdx matrix is correct

  use adjointpetsc
  use BCTypes
  use blockPointers
  use flowvarrefstate ! nw
  use communication
  use inputPhysics
  use inputTimeSpectral
  ! Subroutine Arguments
  implicit none

  real(kind=realType),intent(inout)  :: pts(3,npts,nts) ! Used for fd
                                                    ! check so need
                                                    ! out
  integer(kind=intType), intent(in) :: npts,nTS
  
  ! Local Variables
  real(kind=realType)        :: forces0(3,npts,nts)
  real(kind=realType)        :: forces(3,npts,nts)
  real(kind=realType)        :: forcesAdj(3,npts,nts)
  real(kind=realType)        :: deriv(3,npts,nts)
  real(kind=realType)        :: grid_pts(3,3,3)
  real(kind=realType)        :: force(3),moment(3),refPoint(3)
  real(kind=realType)        :: wAdj(2,2,2,nw),wRef
  integer(kind=intType) :: nn,mm,rowStart,rowEnd,ierr,icol,irow,colStart,colEnd
  integer(kind=intType) :: i,j,k,ii,jj,kk,iii,jjj,kkk,l,sps
  integer(kind=intType) :: iBeg,jBeg,iEnd,jEnd,iStride,jStride
  integer(kind=intType) :: ipt,idim,jpt,jdim
  integer(kind=intType) :: lower_left,lower_right,upper_left,upper_right
  real(kind=realType)   :: fact,diff,tol,norm,h,vec_value,rel_err
  real(kind=realType)   :: vval
  real(kind=realType)   :: max_rel_err,max_err,ad_val,fd_val
  integer(kind=intType) :: i_err,j_err,k_err,l_err,err_count

  real(kind=realType) :: cFp(3),cFv(3),cMp(3),cMv(3),yplusmax
  real(kind=realType) :: origForceSum(3,nts),&
       getForceSum(3,nts), &
       ADForceSum(3,nts)
  real(kind=realType) :: origForceTotal(3,nts),&
       getForceTotal(3,nts), &
       ADForceTotal(3,nts)
  !logical :: rightHanded
  ! Get the reference set of forces 

  call getForces(forces0,pts,npts,nTS)

  ! Now compute the forces using the computeForcesAdj routine
  do sps = 1,nTimeIntervalsSpectral
     ii = 0
     refPoint(:) = 0.0
     domains: do nn=1,nDom
        call setPointers(nn,1_intType,sps)
        rightHanded = flowDoms(nn,1_intType,sps)%rightHanded

        bocos: do mm=1,nBocos
           if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or.&
                BCType(mm) == NSWallIsothermal) then

              jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
              iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd
              iStride = iEnd-iBeg+1
              jStride = jEnd-jBeg+1

              do j=jBeg,jEnd
                 do i=iBeg,iEnd

                    grid_pts(:,:,:) = 0.0
                    wAdj(:,:,:,:)   = 0.0

                    do iii =1,2
                       do jjj=1,2
                          lower_left  = ii + iii + (jjj-1)*iStride-istride-1
                          lower_right = ii + iii + (jjj-1)*iStride-istride
                          upper_left  = ii + iii + (jjj  )*iStride-istride-1
                          upper_right = ii + iii + (jjj  )*iStride-istride

                          if (lower_left > 0) then
                             grid_pts(:,iii  ,jjj  ) = pts(:,lower_left,sps)
                          end if

                          if (lower_right > 0) then
                             grid_pts(:,iii+1,jjj  ) = pts(:,lower_right,sps)
                          end if

                          if (upper_left > 0) then
                             grid_pts(:,iii  ,jjj+1) = pts(:,upper_left,sps)
                          end if

                          if (upper_right > 0) then
                             grid_pts(:,iii+1,jjj+1) = pts(:,upper_right,sps)
                          end if

                       end do
                    end do

                    ! Copy over the states

                    select case (BCFaceID(mm))
                    case (iMin)
                       fact = -1_realType
                       do kkk=1,2
                          wadj(kkk,1,1,:) = w(kkk+1,i  ,j  ,:)
                          wadj(kkk,2,1,:) = w(kkk+1,i+1,j  ,:)
                          wadj(kkk,1,2,:) = w(kkk+1,i  ,j+1,:)
                          wadj(kkk,2,2,:) = w(kkk+1,i+1,j+1,:)
                       end do
                    case (iMax)
                       fact = 1_realType
                       do kkk=1,2
                          wadj(kkk,1,1,:) = w(ib-kkk-1,i  ,j  ,:)
                          wadj(kkk,2,1,:) = w(ib-kkk-1,i+1,j  ,:)
                          wadj(kkk,1,2,:) = w(ib-kkk-1,i  ,j+1,:)
                          wadj(kkk,2,2,:) = w(ib-kkk-1,i+1,j+1,:)
                       end do
                    case (jMin)
                       fact = 1_realType
                       do kkk=1,2
                          wadj(kkk,1,1,:) = w(i  ,kkk+1,j  ,:)
                          wadj(kkk,2,1,:) = w(i+1,kkk+1,j  ,:)
                          wadj(kkk,1,2,:) = w(i  ,kkk+1,j+1,:)
                          wadj(kkk,2,2,:) = w(i+1,kkk+1,j+1,:)
                       end do
                    case (jMax)
                       fact = -1_realType
                       do kkk=1,2
                          wadj(kkk,1,1,:) = w(i  ,jb-kkk-1,j  ,:)
                          wadj(kkk,2,1,:) = w(i+1,jb-kkk-1,j  ,:)
                          wadj(kkk,1,2,:) = w(i  ,jb-kkk-1,j+1,:)
                          wadj(kkk,2,2,:) = w(i+1,jb-kkk-1,j+1,:)
                       end do
                    case (kMin)
                       fact = -1_realType
                       do kkk=1,2
                          wadj(kkk,1,1,:) = w(i  ,j  ,kkk+1,:)
                          wadj(kkk,2,1,:) = w(i+1,j  ,kkk+1,:)
                          wadj(kkk,1,2,:) = w(i  ,j+1,kkk+1,:)
                          wadj(kkk,2,2,:) = w(i+1,j+1,kkk+1,:)
                       end do
                    case (kMax)
                       fact = 1_realType
                       do kkk=1,2
                          wadj(kkk,1,1,:) = w(i  ,j  ,kb-kkk-1,:)
                          wadj(kkk,2,1,:) = w(i+1,j  ,kb-kkk-1,:)
                          wadj(kkk,1,2,:) = w(i  ,j+1,kb-kkk-1,:)
                          wadj(kkk,2,2,:) = w(i+1,j+1,kb-kkk-1,:)
                       end do
                    end select

                    call computeForcesAdj(force,moment,grid_pts,wAdj,&
                         refPoint,fact,iBeg,iEnd,jBeg,jEnd,i,j,righthanded)
                    ii = ii + 1
                    forcesAdj(:,ii,sps) = force

                 end do
              end do
           end if
        end do bocos
     end do domains
  end do

  ! Check The forces --- Thse should be "exact"
  tol = 5e-13
  do sps = 1,nTimeIntervalsSpectral
     do i=1,size(forces0,2)
        diff = norm(forces0(:,i,sps)-forcesAdj(:,i,sps),3)
        if (diff > tol) then
           print *,'myid,i,diff:',myid,i,diff,forces0(1,i,sps),forcesAdj(1,i,sps)
        end if
     end do
  end do
  ! Also check that the the sum matches the orignal forces and moments
  ! function:
  
  origForceSum = 0.0
  getForceSum = 0.0
  ADForceSum = 0.0
  do sps =1,ntimeintervalsspectral
     do nn=1,nDom
        call setPointers(nn,1_intType,sps)
        call forcesAndMoments(cFp, cFv, cMp, cMv, yplusMax)
        origForceSum(:,sps) = origForceSum(:,sps) + cFp + cFv
     end do
  end do
  
  do sps = 1,ntimeintervalsspectral
     do i=1,npts
        getForceSum(:,sps) = getForceSum(:,sps) + forces0(:,i,sps)
        ADForceSum(:,sps) = ADForceSum(:,sps) + forcesAdj(:,i,sps)
     end do
  end do

  fact = two/(gammaInf*pInf*MachCoef*MachCoef*surfaceRef*LRef*LRef*Pref) 
  origForceSum = origForceSum /fact

  ! MPI_REDUCE them and display
  
  call MPI_Reduce(origForceSum, origForceTotal, 3, SUMB_REAL, MPI_SUM, 0, &
       SUMB_COMM_WORLD,ierr)
  call EChk(ierr,__file__,__line__)
  call MPI_Reduce(ADForceSum, ADForceTotal, 3, SUMB_REAL, MPI_SUM, 0, &
       SUMB_COMM_WORLD,ierr)
  call EChk(ierr,__file__,__line__)
  call MPI_Reduce(getForceSum, getForceTotal, 3, SUMB_REAL, MPI_SUM, 0, &
       SUMB_COMM_WORLD,ierr)
  call EChk(ierr,__file__,__line__)

  if (myid == 0) then
     print *, 'Sum Check:'
     print *,'Orig forcesAndMoments        AD Routine              GetForces Routine'
     do sps = 1,ntimeintervalsspectral
        print *,'X:',origForceTotal(1,sps),ADForceTotal(1,sps),getForceTotal(1,sps),sps
        print *,'Y:',origForceTotal(2,sps),ADForceTotal(2,sps),getForceTotal(2,sps),sps
        print *,'Z:',origForceTotal(3,sps),ADForceTotal(3,sps),getForceTotal(3,sps),sps
     end do
  end if

  ! Now Check the dfdx matrix. First call setupCouplingMatrixStruct
  ! to generate the two matrices
 
!  call setupCouplingMatrixStruct(pts,npts,nTS)

  ! -----------------------------
  !           Check dFdx
  ! -----------------------------

  ! Currently Broken

  call MatGetOwnershipRange(dFdx,rowStart,rowEnd,ierr)
  call EChk(ierr,__file__,__line__)
  call MatGetOwnershipRangeColumn(dFdx,colStart,colEnd,ierr)
  call EChk(ierr,__file__,__line__)

  ! if (myid == 0) then 
!      print *,'------------ dfdx verification ----------'
!   end if
!   tol = 1e-6
!   h = 1e-7
  
!   do ipt = 1,npts   ! ----> Loop over Columns
!      do idim = 1,3  ! 
!         pts(idim,ipt) = pts(idim,ipt) + h
!         call getForces(forces,pts,npts)
!         deriv = (forces-forces0)/h
        
!         ! Now check this deriv with the local COLUMN from dFdx
!         icol = colStart + (ipt-1)*3 + idim -1

!         do jpt =1,npts   ! ----> Loop over the Rows 
!            do jdim =1,3  !
!               irow = rowStart+(jpt-1)*3+jdim-1
!               call MatGetValues(dFdx,1,irow,1,icol,vval,ierr)
!               call EChk(ierr,__file__,__line__)

!               diff = abs(vval-deriv(jdim,jpt))
!               if (myid == 0) then
!               if (diff > tol) then
!                  print *,'proc,id,diff:',myid,irow,icol,diff
!                  print *,vval,deriv(jdim,jpt)
!               end if
!            end if
!            end do
!         end do
   
!         pts(idim,ipt) = pts(idim,ipt) - h
!      end do
!   end do

!   call mpi_barrier(sumb_comm_world,ierr)
!   call EChk(ierr,__file__,__line__)
  ! -----------------------------
  !           Check dFdw
  ! -----------------------------

  if (myid == 0) then 
     print *,'------------ dfdw verification ----------'
  end if
!!$  tol = 1e-5
!!$  h = 1e-8
!!$  call MatGetOwnershipRange(dFdw,rowStart,rowEnd,ierr)
!!$  call EChk(ierr,__file__,__line__)
!!$  do nn=1,ndom
!!$     max_rel_err = 0.0
!!$     rel_err     = 0.0
!!$     i_err = 2
!!$     j_err = 2
!!$     k_err = 2
!!$     l_err = 1
!!$     err_count = 0
!!$     ad_val = 0.0
!!$     fd_val = 0.0
!!$     call setPointersAdj(nn,1_intType,1_intType)!-| 
!!$     do l=1,nw                                  ! |
!!$        do k=2,kl                               ! |->Loop Over Columns
!!$           do j=2,jl                            ! | 
!!$              do i=2,il                         ! |
!!$               
!!$                 w(i,j,k,l) = w(i,j,k,l) + h
!!$
!!$                 call computePressureAdjFullBlock(w,p) ! Full block
!!$                                                    ! version
!!$                 call applyAllBC(.True.)
!!$                 call setPointersAdj(nn,1_intType,1_intType) 
!!$
!!$                 call getForces(forces,pts,npts)
!!$                 call setPointersAdj(nn,1_intType,1_intType) 
!!$                                                             
!!$                 deriv = 0.0
!!$                 deriv = (forces-forces0)/h
!!$                 
!!$                 ! Now check this deriv with the local COLUMN from dFdw
!!$
!!$                 icol = globalCell(i,j,k)*nw+l-1
!!$              
!!$                 do jpt =1,npts   ! ----> Loop over the Rows 
!!$                       do jdim =1,3  !
!!$                          irow = rowStart+(jpt-1)*3+jdim-1
!!$                          call MatGetValues(dFdw,1,irow,1,icol,vec_value,ierr)
!!$                          call EChk(ierr,__file__,__line__)
!!$                          diff = abs(vec_value-deriv(jdim,jpt,sps))
!!$                          if (diff > 1e-16) then
!!$                             rel_err = diff/((vec_value+deriv(jdim,jpt,sps)))
!!$                          else
!!$                             rel_err = 0.0_realType
!!$                          end if
!!$
!!$                          if (rel_err > tol) then
!!$                             err_count = err_count + 1
!!$                          end if
!!$
!!$                          if (rel_err > tol .and. rel_err > max_rel_err) then
!!$                             max_rel_err = rel_err
!!$                             max_err     = diff
!!$                             fd_val      = deriv(jdim,jpt,sps)
!!$                             ad_val      = vec_value
!!$                             i_err = i
!!$                             j_err = j
!!$                             k_err = k
!!$                             l_err = l
!!$                          end if
!!$                       end do
!!$                    end do
!!$
!!$                 w(i,j,k,l) = w(i,j,k,l)-h
!!$                 ! Must call this to reset the pressures
!!$                 call computePressureAdjFullBlock(w,p)
!!$              end do ! state loop
!!$           end do !i loop
!!$        end do ! j loop
!!$     end do ! k loop
!!$
!!$     write(*,900)'P:',myid,' B:',nn,' Nerr:',err_count,' Max Rel:',max_rel_err,&
!!$          ' Max Err:',max_err,' AD:',ad_val,' FD:',fd_val,' Pos:',i_err,' ',j_err,' ',k_err,' ',l_err
!!$900  format (A,I2,A,I2,A,I5,A,G10.5,A,G10.5,A,G10.5,A,G10.5,A,I2,A,I2,A,I2,A,I2)
!!$
!!$  end do ! ndom loop

end subroutine verifyForces

function norm(X,n)
  ! Compute the L2 nomr of X
  implicit none
  double precision       :: X(n)
  double precision       :: norm
  integer                :: i,n
  norm = 0.0
  do i=1,n
     norm = norm + X(i)**2
  end do
  norm = sqrt(norm)
end function norm
