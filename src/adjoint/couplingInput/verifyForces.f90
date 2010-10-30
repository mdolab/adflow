subroutine verifyForces(pts,npts)

  ! This routine does three thing:
  ! 1. Check that getForces and computeForceCouplingAdj give the same results
  ! 2. Check that dSdw matrix is correct
  ! 3. Check that dSdx matrix is correct

  use adjointpetsc
  use BCTypes
  use blockPointers
  use flowvarrefstate ! nw
  use communication
  ! Subroutine Arguments
  implicit none

  real(kind=realType)  :: pts(3,npts)
  integer(kind=intType), intent(in) :: npts
  
  ! Local Variables
  real(kind=realType)        :: forces0(3,npts)
  real(kind=realType)        :: forces(3,npts)
  real(kind=realType)        :: forcesAdj(3,npts)
  real(kind=realType)        :: deriv(3,npts)
  real(kind=realType)        :: grid_pts(3,3,3)
  real(kind=realType)        :: force(3)
  real(kind=realType)        :: wAdj(2,2,2,nw),wRef
  integer(kind=intType) :: nn,mm,rowStart,rowEnd,ierr,icol,irow,colStart,colEnd
  integer(kind=intType) :: i,j,k,ii,jj,kk,iii,jjj,kkk,l
  integer(kind=intType) :: iBeg,jBeg,iEnd,jEnd,iStride,jStride
  integer(kind=intType) :: ipt,idim,jpt,jdim
  integer(kind=intType) :: lower_left,lower_right,upper_left,upper_right
  real(kind=realType)   :: fact,diff,tol,norm,h,vec_value,rel_err
  real(kind=realType)   :: max_rel_err,max_err,ad_val,fd_val
  integer(kind=intType) :: i_err,j_err,k_err,l_err,err_count
  !logical :: rightHanded
  ! Get the reference set of forces 

  call getForces(forces0,pts,npts)


  ! Now compute the forces using the computeForceCouplingAdj routine

  ii = 0

  domains: do nn=1,nDom
     call setPointers(nn,1_intType,1_intType)
     rightHanded = flowDoms(nn,1_intType,1_intType)%rightHanded
     
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
                          grid_pts(:,iii  ,jjj  ) = pts(:,lower_left)
                       end if

                       if (lower_right > 0) then
                          grid_pts(:,iii+1,jjj  ) = pts(:,lower_right)
                       end if
                       
                       if (upper_left > 0) then
                          grid_pts(:,iii  ,jjj+1) = pts(:,upper_left)
                       end if

                       if (upper_right > 0) then
                          grid_pts(:,iii+1,jjj+1) = pts(:,upper_right)
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
                    fact = -1_realType
                    do kkk=1,2
                       wadj(kkk,1,1,:) = w(i  ,kkk+1,j  ,:)
                       wadj(kkk,2,1,:) = w(i+1,kkk+1,j  ,:)
                       wadj(kkk,1,2,:) = w(i  ,kkk+1,j+1,:)
                       wadj(kkk,2,2,:) = w(i+1,kkk+1,j+1,:)
                    end do
                 case (jMax)
                    fact = 1_realType
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

                 call computeForceCouplingAdj(force,grid_pts,wAdj,fact,&
                      iBeg,iEnd,jBeg,jEnd,i,j,righthanded)
                 ii = ii + 1
                 forcesAdj(:,ii) = force

              end do
           end do
        end if
     end do bocos
  end do domains

  ! Check The forces --- Thse should be "exact"
  tol = 5e-13
  do i=1,size(forces0,2)
     diff = norm(forces0(:,i)-forcesAdj(:,i),3)
     if (diff > tol) then
        print *,'myid,i,diff:',myid,i,diff
     end if
  end do

  ! Now Check the dsdx matrix. First call setupCouplingMatrixStruct
  ! to generate the two matrices
 
  call setupCouplingMatrixStruct(pts,npts)

  ! -----------------------------
  !           Check dSdx
  ! -----------------------------

  call MatGetOwnershipRange(dSdx,rowStart,rowEnd,ierr)
  call MatGetOwnershipRangeColumn(dSdx,colStart,colEnd,ierr)
  if (myid == 0) then 
     print *,'------------ dsdx verification ----------'
  end if
  tol = 1e-6
  h = 1e-7
  do ipt = 1,npts   ! ----> Loop over Columns
     do idim = 1,3  ! 
        pts(idim,ipt) = pts(idim,ipt) + h
        call getForces(forces,pts,npts)
        deriv = (forces-forces0)/h
        
        ! Now check this deriv with the local COLUMN from dSdx
        icol = colStart + (ipt-1)*3 + idim -1

        do jpt =1,npts   ! ----> Loop over the Rows 
           do jdim =1,3  !
              irow = rowStart+(jpt-1)*3+jdim-1
              call MatGetValues(dSdx,1,irow,1,icol,vec_value,ierr)

              diff = abs(vec_value-deriv(jdim,jpt))
              if (myid == 0) then
              if (diff > tol) then
                 print *,'proc,id,diff:',myid,irow,icol,diff
              end if
           end if
           end do
        end do
   
        pts(idim,ipt) = pts(idim,ipt) - h
     end do
  end do

  call mpi_barrier(sumb_comm_world,ierr)
  ! -----------------------------
  !           Check dSdw
  ! -----------------------------

  if (myid == 0) then 
     print *,'------------ dsdw verification ----------'
  end if
  tol = 1e-5
  h = 1e-8
  call MatGetOwnershipRange(dSdw,rowStart,rowEnd,ierr)
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
     call setPointersAdj(nn,1_intType,1_intType)!-| 
     do l=1,nw                                  ! |
        do k=2,kl                               ! |->Loop Over Columns
           do j=2,jl                            ! | 
              do i=2,il                         ! |
               
                 w(i,j,k,l) = w(i,j,k,l) + h

                 call computeForcesPressureAdj(w,p) ! Full block
                                                    ! version
                 call applyAllBC(.True.)
                 call setPointersAdj(nn,1_intType,1_intType) 

                 call getForces(forces,pts,npts)
                 call setPointersAdj(nn,1_intType,1_intType) 
                                                             
                 deriv = 0.0
                 deriv = (forces-forces0)/h
                 
                 ! Now check this deriv with the local COLUMN from dSdw

                 icol = globalCell(i,j,k)*nw+l-1
              
                 do jpt =1,npts   ! ----> Loop over the Rows 
                       do jdim =1,3  !
                          irow = rowStart+(jpt-1)*3+jdim-1
                          call MatGetValues(dSdw,1,irow,1,icol,vec_value,ierr)
                          
                          diff = abs(vec_value-deriv(jdim,jpt))
                          if (diff > 1e-16) then
                             rel_err = diff/((vec_value+deriv(jdim,jpt)))
                          else
                             rel_err = 0.0_realType
                          end if

                          if (rel_err > tol) then
                             err_count = err_count + 1
                          end if

                          if (rel_err > tol .and. rel_err > max_rel_err) then
                             max_rel_err = rel_err
                             max_err     = diff
                             fd_val      = deriv(jdim,jpt)
                             ad_val      = vec_value
                             i_err = i
                             j_err = j
                             k_err = k
                             l_err = l
                          end if
                       end do
                    end do

                 w(i,j,k,l) = w(i,j,k,l)-h
                 ! Must call this to reset the pressures
                 call computeForcesPressureAdj(w,p)
              end do ! state loop
           end do !i loop
        end do ! j loop
     end do ! k loop

     write(*,900)'P:',myid,' B:',nn,' Nerr:',err_count,' Max Rel:',max_rel_err,&
          ' Max Err:',max_err,' AD:',ad_val,' FD:',fd_val,' Pos:',i_err,' ',j_err,' ',k_err,' ',l_err
900  format (A,I2,A,I2,A,I5,A,G10.5,A,G10.5,A,G10.5,A,G10.5,A,I2,A,I2,A,I2,A,I2)

  end do ! ndom loop

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
