subroutine verifyForces(pts,npts,nTS)

  ! This routine does three thing:
  ! 1. Check that getForces and computeForcesAdj give the same results
  ! 2. Check that dFdw matrix is correct
  ! 3. Check that dFdx matrix is correct
#ifndef USE_NO_PETSC
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
  integer(kind=intType) :: i,j,k,ii,jj,kk,iii,jjj,kkk,l,sps,sps2
  integer(kind=intType) :: iBeg,jBeg,iEnd,jEnd,iStride,jStride
  integer(kind=intType) :: ipt,idim,jpt,jdim
  integer(kind=intType) :: lower_left,lower_right,upper_left,upper_right
  real(kind=realType)   :: fact,diff,tol,norm,h,vec_value,rel_err
  real(kind=realTYpe)   :: relDiff,adiff
  real(kind=realType)   :: vval
  real(kind=realType)   :: max_rel_err,max_err,ad_val,fd_val
  integer(kind=intType) :: i_err,j_err,k_err,l_err,err_count,pt_high,ind

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

                    ! This iii,jjj loop is over the nodes surrounding
                    ! the pt of interst
                    do iii =-1,1
                       do jjj=-1,1

                          ! ii is the node we're after...ACTUALLY is 1
                          ! less than the node we're intersted in
                          ! since ii starts at 0 above. The two loops
                          ! set the 9 coordindates into grid_pts

                          ind = jjj*istride + iii + 1

                          if (i+iii > 0 .and. i + iii <= iEnd .and. &
                              j+jjj > 0 .and. j + jjj <= jend) then
                             grid_pts(:,iii+2,jjj+2) = pts(:,ind+ii,sps)
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
  tol = 1e-12
  do sps = 1,nTimeIntervalsSpectral
     do i=1,npts
        diff = norm(forces0(:,i,sps)-forcesAdj(:,i,sps),3)
        if (diff > tol) then
           print *,'myid,i,diff:',myid,i,diff
           print *,'coord1:',forces0(1,i,sps),forcesAdj(1,i,sps)
           print *,'coord2:',forces0(2,i,sps),forcesAdj(2,i,sps)
           print *,'coord3:',forces0(3,i,sps),forcesAdj(3,i,sps)
           
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
  call EChk(ierr,__FILE__,__LINE__)
  call MPI_Reduce(ADForceSum, ADForceTotal, 3, SUMB_REAL, MPI_SUM, 0, &
       SUMB_COMM_WORLD,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call MPI_Reduce(getForceSum, getForceTotal, 3, SUMB_REAL, MPI_SUM, 0, &
       SUMB_COMM_WORLD,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  if (myid == 0) then
     print *, 'Sum Check:'
     print *,'Orig forcesAndMoments        AD Routine              GetForces Routine'
     do sps = 1,ntimeintervalsspectral
        print *,'SPS:',sps
        print *,'X:',origForceTotal(1,sps),ADForceTotal(1,sps),getForceTotal(1,sps)
        print *,'Y:',origForceTotal(2,sps),ADForceTotal(2,sps),getForceTotal(2,sps)
        print *,'Z:',origForceTotal(3,sps),ADForceTotal(3,sps),getForceTotal(3,sps)
     end do
  end if
 
  ! Now Check the dfdx matrix. First call setupCouplingMatrixStruct
  ! to generate the two matrices
  call setupCouplingMatrixStruct(pts,npts,nTS)


  !  -----------------------------
  !            Check dFdx
  !  -----------------------------

  call MatGetOwnershipRange(dFdx,rowStart,rowEnd,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call MatGetOwnershipRangeColumn(dFdx,colStart,colEnd,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  if (myid == 0) then 
     print *,'------------ dfdx verification ----------'
  end if
  tol = 1e-5
  h = 1e-7
  
  do sps =1,nTimeIntervalsSpectral
     do ipt = 1,npts   ! ----> Loop over Columns
        do idim = 1,3  ! 
           pts(idim,ipt,sps) = pts(idim,ipt,sps) + h
           call getForces(forces,pts,npts,nTS)
           deriv = (forces-forces0)/h
        
           ! Now check this deriv with the local COLUMN from dFdx
           icol = colStart + (ipt-1)*3 + idim -1

           do sps2=1,nTimeIntervalsSpectral
              do jpt =1,npts   ! ----> Loop over the Rows 
                 do jdim =1,3  !
                    irow = rowStart+(jpt-1)*3+jdim-1+(sps2-1)*npts*3
                    call MatGetValues(dFdx,1,irow,1,icol,vval,ierr)
                    call EChk(ierr,__FILE__,__LINE__)
                    
                    adiff = abs(vval-deriv(jdim,jpt,sps2))
                    reldiff = abs(vval-deriv(jdim,jpt,sps2))/(vval + 1e-16)
                    if (myid == 0) then
                       if (reldiff > tol) then
                          print *,'proc,id,diff:',myid,irow,icol,reldiff
                          print *,vval,deriv(jdim,jpt,sps2)
                       end if
                    end if
                 end do
              end do
           end do
           pts(idim,ipt,sps) = pts(idim,ipt,sps) - h
        end do
     end do
  end do
  call mpi_barrier(sumb_comm_world,ierr)
  call EChk(ierr,__FILE__,__LINE__)
 
  ! -----------------------------
  !           Check dFdw
  ! -----------------------------

  if (myid == 0) then 
     print *,'------------ dfdw verification ----------'
  end if
 
  tol = 1e-5
  h = 1e-8
  call MatGetOwnershipRange(dFdw,rowStart,rowEnd,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  do sps=1,nTimeIntervalsSpectral
     do nn=1,ndom
        max_rel_err = 0.0
        rel_err     = 0.0
        i_err = 2
        j_err = 2
        k_err = 2
        l_err = 1
        err_count = 0
        ad_val = 0.0
        fd_val = 0.0
        call setPointersAdj(nn,1_intType,sps)!-| 
        do l=1,nw                            ! |
           do k=2,kl                         ! |->Loop Over Columns
              do j=2,jl                      ! | 
                 do i=2,il                   ! |
               
                    w(i,j,k,l) = w(i,j,k,l) + h

                    call computePressureAdjFullBlock(w,p) ! Full block
                    ! version
                    call applyAllBC(.True.)
                    call setPointersAdj(nn,1_intType,sps)
                    
                    call getForces(forces,pts,npts,nTS)
                    call setPointersAdj(nn,1_intType,sps)
                    
                    deriv = 0.0
                    deriv = (forces-forces0)/h
                 
                    ! Now check this deriv with the local COLUMN from dFdw

                    icol = globalCell(i,j,k)*nw+l-1
              
                    do sps2=1,ntimeIntervalsSpectral
                       do jpt =1,npts   ! ----> Loop over the Rows 
                          do jdim =1,3  !
                             irow = rowStart+(jpt-1)*3+jdim-1 + (sps2-1)*npts*3
                             call MatGetValues(dFdw,1,irow,1,icol,vec_value,ierr)
                             call EChk(ierr,__FILE__,__LINE__)
                             diff = abs(vec_value-deriv(jdim,jpt,sps2))
                             if (diff > 1e-16) then
                                rel_err = diff/((vec_value+deriv(jdim,jpt,sps2)))
                             else
                                rel_err = 0.0_realType
                             end if
                             
                             if (rel_err > tol) then
                                err_count = err_count + 1
                             end if

                             if (rel_err > tol .and. rel_err > max_rel_err) then
                                max_rel_err = rel_err
                                max_err     = diff
                                fd_val      = deriv(jdim,jpt,sps2)
                                ad_val      = vec_value
                                i_err = i
                                j_err = j
                                k_err = k
                                l_err = l
                             end if
                          end do ! jdim loop
                       end do ! Pt loop
                    end do ! sps loop
                    w(i,j,k,l) = w(i,j,k,l)-h
                    ! Must call this to reset the pressures
                    call computePressureAdjFullBlock(w,p)
                 end do !i loop
              end do ! j loop
           end do ! k loop
        end do ! State Loop
        write(*,900)'P:',myid,' B:',nn,' Nerr:',err_count,' Max Rel:',max_rel_err,&
             ' Max Err:',max_err,' AD:',ad_val,' FD:',fd_val,' Pos:',i_err,' ',j_err,' ',k_err,' ',l_err
900     format (A,I2,A,I2,A,I5,A,G10.5,A,G10.5,A,G10.5,A,G10.5,A,I2,A,I2,A,I2,A,I2)
        
     end do ! ndom loop
  end do ! sps loop
#endif
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
