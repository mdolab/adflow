
subroutine computeForceAndMomentAdj(Force,cForce,Lift,Drag,Cl,Cd,&
     moment,cMoment,alphaAdj,betaAdj,liftIndex,MachCoefAdj,&
     pointRefAdj,pts,npts,wBlock,rightHandedAdj,faceID,&
     iBeg,iEnd,jBeg,jEnd,ii_start)
     

  !     ******************************************************************
  !     *                                                                *
  !     * Compute the sum of the forces and moments on all blocks on     *
  !     * this processor. This function can be AD'd                      *
  !     *                                                                *
  !     ******************************************************************
  !
  use blockPointers
  use inputPhysics        ! equations
  use flowVarRefState     ! nw
  use inputDiscretization ! spaceDiscr, useCompactDiss
  use bcTypes             !imin,imax,jmin,jmax,kmin,kmax

  implicit none

  ! Subroutine Arguments
  ! Output
  real(kind=realType), intent(out) :: force(3),cForce(3)
  real(kind=realType), intent(out) :: Lift,Drag,Cl,Cd
  real(kind=realType), intent(out) :: moment(3),cMoment(3)

  ! Input
  real(kind=realType), intent(in) :: alphaAdj,betaAdj
  integer(kind=intType), intent(in):: liftIndex
  real(kind=realType), intent(in) :: MachCoefAdj
  
  real(kind=realType), intent(in) :: pointRefAdj(3)
  real(kind=realType), intent(in) :: pts(3,npts)
  integer(kind=intType),intent(in):: npts
  real(kind=realType), intent(in) :: wBlock(0:ib,0:jb,0:kb,nw)
  logical, intent(in) :: rightHandedAdj
  integer(kind=intType), intent(in):: faceID
  integer(kind=intType), intent(in):: iBeg,iEnd,jBeg,jEnd,ii_start
  ! Local Variables
  integer(kind=intType) :: ii
  real(kind=realType) :: addForce(3),addMoment(3),refPoint(3)
 
  real(kind=realType) :: liftDir(3),dragDir(3),freeStreamDir(3)
  real(kind=realType) :: grid_pts(3,3,3),wAdj(2,2,2,nw)
  
  integer(kind=intType) :: iStride,jStride,i,j
  integer(kind=intType) :: iii,jjj,kkk
  integer(kind=intTYpe) :: lower_left,lower_right,upper_left,upper_right
  real(kind=realType) :: fact
  real(kind=realType) :: velDirFreestreamAdj(3)
  ! Only need to zero force and moment -> these are summed again
  force = 0.0
  moment = 0.0
  
  iStride = iEnd-iBeg+1
  jStride = jEnd-jBeg+1
  ii = ii_start 

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

              if (lower_left > 0 .and. lower_left <= npts) then
                 grid_pts(:,iii  ,jjj  ) = pts(:,lower_left)
              end if
              
              if (lower_right > 0 .and. lower_right <=npts ) then
                 grid_pts(:,iii+1,jjj  ) = pts(:,lower_right)
              end if

              if (upper_left > 0 .and. upper_left <= npts) then
                 grid_pts(:,iii  ,jjj+1) = pts(:,upper_left)
              end if

              if (upper_right > 0 .and. upper_right <= npts) then
                 grid_pts(:,iii+1,jjj+1) = pts(:,upper_right)
              end if
              
           end do
        end do

        !Copy over the states
        
        select case (faceID)
        case (iMin)
           fact = -1_realType
           do kkk=1,2
              wadj(kkk,1,1,:) = wBlock(kkk+1,i  ,j  ,:)
              wadj(kkk,2,1,:) = wBlock(kkk+1,i+1,j  ,:)
              wadj(kkk,1,2,:) = wBlock(kkk+1,i  ,j+1,:)
              wadj(kkk,2,2,:) = wBlock(kkk+1,i+1,j+1,:)
           end do
        case (iMax)
           fact = 1_realType
           do kkk=1,2
              wadj(kkk,1,1,:) = wBlock(ib-kkk-1,i  ,j  ,:)
              wadj(kkk,2,1,:) = wBlock(ib-kkk-1,i+1,j  ,:)
              wadj(kkk,1,2,:) = wBlock(ib-kkk-1,i  ,j+1,:)
              wadj(kkk,2,2,:) = wBlock(ib-kkk-1,i+1,j+1,:)
           end do
        case (jMin)
           fact = 1_realType
           do kkk=1,2
              wadj(kkk,1,1,:) = wBlock(i  ,kkk+1,j  ,:)
              wadj(kkk,2,1,:) = wBlock(i+1,kkk+1,j  ,:)
              wadj(kkk,1,2,:) = wBlock(i  ,kkk+1,j+1,:)
              wadj(kkk,2,2,:) = wBlock(i+1,kkk+1,j+1,:)
           end do
        case (jMax)
           fact = -1_realType
           do kkk=1,2
              wadj(kkk,1,1,:) = wBlock(i  ,jb-kkk-1,j  ,:)
              wadj(kkk,2,1,:) = wBlock(i+1,jb-kkk-1,j  ,:)
              wadj(kkk,1,2,:) = wBlock(i  ,jb-kkk-1,j+1,:)
              wadj(kkk,2,2,:) = wBlock(i+1,jb-kkk-1,j+1,:)
           end do
        case (kMin)
           fact = -1_realType
           do kkk=1,2
              wadj(kkk,1,1,:) = wBlock(i  ,j  ,kkk+1,:)
              wadj(kkk,2,1,:) = wBlock(i+1,j  ,kkk+1,:)
              wadj(kkk,1,2,:) = wBlock(i  ,j+1,kkk+1,:)
              wadj(kkk,2,2,:) = wBlock(i+1,j+1,kkk+1,:)
           end do
        case (kMax)
           fact = 1_realType
           do kkk=1,2
              wadj(kkk,1,1,:) = wBlock(i  ,j  ,kb-kkk-1,:)
              wadj(kkk,2,1,:) = wBlock(i+1,j  ,kb-kkk-1,:)
              wadj(kkk,1,2,:) = wBlock(i  ,j+1,kb-kkk-1,:)
              wadj(kkk,2,2,:) = wBlock(i+1,j+1,kb-kkk-1,:)
           end do
        end select
        
        call computeForcesAdj(addForce,addMoment,&
             grid_pts,wAdj,pointRefAdj,fact,iBeg,iEnd,jBeg,jEnd,&
             i,j,righthandedAdj)

        ii = ii + 1
        force = force + addForce
        moment = moment + addMoment
        
     end do
  end do
  
  ! Now we know the sum of the force and moment contribution from this block
  
  ! First get cForce -> Coefficient of FOrce
  fact = two/(gammaInf*pInf*pRef*MachCoefAdj*MachCoefAdj*surfaceRef*LRef*LRef)
 
  cForce = fact*Force

  ! To get Lift,Drag,Cl and Cd get lift and drag directions


  call adjustInflowAngleForcesAdj(alphaAdj,betaAdj,velDirFreestreamAdj,&
       liftDir,dragDir,liftIndex)

  ! Take Dot Products ... this won't AD properly so we will write explictly
  !Lift = dot_product(Force,liftDir)
  !Drag = dot_product(Force,dragDir)

  Lift = Force(1)*liftDir(1) + Force(2)*liftDir(2) + Force(3)*liftDir(3)
  Drag = Force(1)*dragDir(1) + Force(2)*dragDir(2) + Force(3)*dragDir(3)

  Cl = Lift * fact
  Cd = Drag * fact

  ! Update fact for moment normalization
  fact = fact/(lengthRef*LRef)
  
  cMoment = moment * fact

end subroutine computeForceAndMomentAdj

