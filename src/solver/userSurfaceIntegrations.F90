module userSurfaceIntegrations

  use constants
  use communication, only : commType, internalCommType
  use userSurfaceIntegrationData
contains

  subroutine integrateUserSurfaces(localValues, famList, sps, withGathered, funcValues)

    use constants
    use block, onlY : flowDoms, nDom
    use flowVarRefState, only : pRef, rhoRef, pRef, timeRef, LRef, TRef
    use communication, only : myid, adflow_comm_world
    use utils, only : EChk, mynorm2
    use flowUtils, only : computePtot, computeTtot
    use sorting, only : famInList
    use zipperIntegrations, only : flowIntegrationZipper
    use utils, only : terminate
    implicit none

    ! Input Parameters
    real(kind=realType), dimension(nLocalValues), intent(inout) :: localValues
    integer(kind=intType), dimension(:), intent(in) :: famList
    integer(kind=intType), intent(in) :: sps
    logical, intent(in) :: withGathered
    real(kind=realType), intent(in), dimension(:) :: funcValues

    ! Working parameters
    integer(kind=intType) :: iSurf, i, j, k, jj, ierr, nn, iDim, nPts
    real(kind=realType), dimension(:), allocatable :: recvBuffer1, recvBuffer2
    real(kind=realType), dimension(:, :), allocatable :: vars
    integer(kind=intType), dimension(:), allocatable :: fams
    logical, dimension(:), allocatable :: ptValid
    type(userIntSurf), pointer :: surf

    ! Set the pointers for the required communication variables
    domainLoop:do nn=1, nDom
       if (flowDoms(nn, 1, sps)%addGridVelocities) then 
          call terminate("userSurfaceIntegrations", "Cannot use user-supplied surface integrations"& 
               &"on with moving grids")
       end if
       
       flowDoms(nn, 1, sps)%realCommVars(iRho)%var => flowDoms(nn, 1, sps)%w(:, :, :, iRho)
       flowDoms(nn, 1, sps)%realCommVars(iVx)%var => flowDoms(nn, 1, sps)%w(:, :, :, iVx)
       flowDoms(nn, 1, sps)%realCommVars(iVy)%var => flowDoms(nn, 1, sps)%w(:, :, :, iVy)
       flowDoms(nn, 1, sps)%realCommVars(iVz)%var => flowDoms(nn, 1, sps)%w(:, :, :, iVz)
       flowDoms(nn, 1, sps)%realCommVars(iZippFlowP)%var => flowDoms(nn, 1, sps)%P(:, :, :)
       flowDoms(nn, 1, sps)%realCommVars(iZippFlowGamma)%var => flowDoms(nn, 1, sps)%gamma(:, :, :)
       ! flowDoms(nn, 1, sps)%realCommVars(iZippFlowSface)%var => Not Implemented
    
       flowDoms(nn, 1, sps)%realCommVars(iZippFlowX)%var => flowDoms(nn, 1, sps)%x(:, :, :, 1)
       flowDoms(nn, 1, sps)%realCommVars(iZippFlowY)%var => flowDoms(nn, 1, sps)%x(:, :, :, 2)
       flowDoms(nn, 1, sps)%realCommVars(iZippFlowZ)%var => flowDoms(nn, 1, sps)%x(:, :, :, 3)
    end do domainLoop

    masterLoop: do iSurf=1, nUserIntSurfs

       ! Pointer for easier reading
       surf => userIntSurfs(iSurf)

       ! We will make a short-cut here: By definition user supplied
       ! surfaces have a fixed family, we won't do anything if we
       ! are not told to deal with this surface. 

       famInclude: if (famInList(surf%famID, famList)) then 
          nPts = size(surf%pts, 2)
             
          ! Communicate the face values and the nodal values
          if (myid == 0) then 
             allocate(recvBuffer1(6*nPts), recvBuffer2(3*nPts))
          else
             allocate(recvBuffer1(0), recvBuffer2(0))
          end if

          call commUserIntegrationSurfaceVars(recvBuffer1, iRho, iZippFlowGamma, surf%flowComm)
          call commUserIntegrationSurfaceVars(recvBuffer2, iZippFlowX, iZippFlowZ, surf%nodeComm)

          ! *Finally* we can do the actual integrations
          if (myid == 0)  then 

             ! Allocate some temporary data needed to supply to the
             ! zipper integration routine. 
             allocate(ptValid(npts), vars(npts, nZippFlowComm), fams(size(surf%conn, 2)))

             ! Prepare for the "zipper" integration call. We have to
             ! re-order the data according to the "inv" array in each
             ! of the two comms. 
             do i=1, nPts

                ! Flow Variables
                j = surf%flowComm%inv(i)
                vars(j, iRho:iZippFlowGamma) = recvBuffer1(6*(i-1) + iRho : 6*(i-1) + iZippFlowGamma)

                ! Sface is not implemented. To correctly do this,
                ! interpolate the three components of 's', do the dot
                ! product with the local normal to get the sFace value. 
                vars(j, iZippFlowSface)  = zero

                ! Node Comm Values
                j = surf%nodeComm%inv(i)
                vars(j, iZippFlowX:iZippFlowZ) = recvBuffer2(3*i-2:3*i)
                      
                ! The additional pt-valid array
                ptValid(j) = surf%nodeComm%valid(i)
             end do
             
             ! The family array is all the same value:
             fams = surf%famID

             ! Perform the actual integration
             call flowIntegrationZipper(.True., surf%conn, fams, vars, localValues, famList, sps, &
                  withGathered, funcValues, surf%nodeComm%Valid)
             deallocate(ptValid, vars, fams)
          end if
          deallocate(recvBuffer1, recvBuffer2)
       end if famInclude
    end do masterLoop
  end subroutine integrateUserSurfaces
#ifndef USE_COMPLEX
  subroutine integrateUserSurfaces_d(localValues, localValuesd, famList, sps, withGathered, &
       funcValues, funcValuesd)

    use constants
    use block, onlY : flowDoms, flowDomsd, nDom
    use flowVarRefState, only : pRef, rhoRef, pRef, timeRef, LRef, TRef
    use communication, only : myid, adflow_comm_world
    use utils, only : EChk, mynorm2
    use flowUtils, only : computePtot, computeTtot
    use sorting, only : famInList
    use zipperIntegrations_d, only : flowIntegrationZipper_d
    use utils, only : terminate
    implicit none

    ! Input Parameters
    real(kind=realType), dimension(nLocalValues), intent(inout) :: localValues, localValuesd
    integer(kind=intType), dimension(:), intent(in) :: famList
    integer(kind=intType), intent(in) :: sps
    logical, intent(in) :: withGathered
    real(kind=realType), intent(in), dimension(:) :: funcValues, funcValuesd

    ! Working parameters
    integer(kind=intType) :: iSurf, i, j, k, jj, ierr, nn, iDim, nPts
    real(kind=realType), dimension(:), allocatable :: recvBuffer1, recvBuffer2
    real(kind=realType), dimension(:), allocatable :: recvBuffer1d, recvBuffer2d
    real(kind=realType), dimension(:, :), allocatable :: vars, varsd 
    integer(kind=intType), dimension(:), allocatable :: fams
    logical, dimension(:), allocatable :: ptValid
    type(userIntSurf), pointer :: surf

    ! Set the pointers for the required communication variables
    domainLoop:do nn=1, nDom
       if (flowDoms(nn, 1, sps)%addGridVelocities) then 
          call terminate("userSurfaceIntegrations", "Cannot use user-supplied surface integrations"& 
               &"on with moving grids")
       end if
       
       flowDoms(nn, 1, sps)%realCommVars(iRho)%var => flowDoms(nn, 1, sps)%w(:, :, :, iRho)
       flowDoms(nn, 1, sps)%realCommVars(iVx)%var => flowDoms(nn, 1, sps)%w(:, :, :, iVx)
       flowDoms(nn, 1, sps)%realCommVars(iVy)%var => flowDoms(nn, 1, sps)%w(:, :, :, iVy)
       flowDoms(nn, 1, sps)%realCommVars(iVz)%var => flowDoms(nn, 1, sps)%w(:, :, :, iVz)
       flowDoms(nn, 1, sps)%realCommVars(iZippFlowP)%var => flowDoms(nn, 1, sps)%P(:, :, :)
       flowDoms(nn, 1, sps)%realCommVars(iZippFlowGamma)%var => flowDoms(nn, 1, sps)%gamma(:, :, :)
       ! flowDoms(nn, 1, sps)%realCommVars(iZippFlowSface)%var => Not Implemented
    
       flowDoms(nn, 1, sps)%realCommVars(iZippFlowX)%var => flowDoms(nn, 1, sps)%x(:, :, :, 1)
       flowDoms(nn, 1, sps)%realCommVars(iZippFlowY)%var => flowDoms(nn, 1, sps)%x(:, :, :, 2)
       flowDoms(nn, 1, sps)%realCommVars(iZippFlowZ)%var => flowDoms(nn, 1, sps)%x(:, :, :, 3)

       flowDoms(nn, 1, sps)%realCommVars(iRho+nZippFlowComm)%var => flowDomsd(nn, 1, sps)%w(:, :, :, iRho)
       flowDoms(nn, 1, sps)%realCommVars(iVx+nZippFlowComm)%var => flowDomsd(nn, 1, sps)%w(:, :, :, iVx)
       flowDoms(nn, 1, sps)%realCommVars(iVy+nZippFlowComm)%var => flowDomsd(nn, 1, sps)%w(:, :, :, iVy)
       flowDoms(nn, 1, sps)%realCommVars(iVz+nZippFlowComm)%var => flowDomsd(nn, 1, sps)%w(:, :, :, iVz)
       flowDoms(nn, 1, sps)%realCommVars(iZippFlowP+nZippFlowComm)%var => flowDomsd(nn, 1, sps)%P(:, :, :)
       flowDoms(nn, 1, sps)%realCommVars(iZippFlowGamma+nZippFlowComm)%var => flowDomsd(nn, 1, sps)%gamma(:, :, :)
       ! flowDoms(nn, 1, sps)%realCommVars(iZippFlowSface+nZippFlowComm)%var => Not Implemented
    
       flowDoms(nn, 1, sps)%realCommVars(iZippFlowX+nZippFlowComm)%var => flowDomsd(nn, 1, sps)%x(:, :, :, 1)
       flowDoms(nn, 1, sps)%realCommVars(iZippFlowY+nZippFlowComm)%var => flowDomsd(nn, 1, sps)%x(:, :, :, 2)
       flowDoms(nn, 1, sps)%realCommVars(iZippFlowZ+nZippFlowComm)%var => flowDomsd(nn, 1, sps)%x(:, :, :, 3)

    end do domainLoop

    masterLoop: do iSurf=1, nUserIntSurfs

       ! Pointer for easier reading
       surf => userIntSurfs(iSurf)

       ! We will make a short-cut here: By definition user supplied
       ! surfaces have a fixed family, we won't do anything if we
       ! are not told to deal with this surface. 

       famInclude: if (famInList(surf%famID, famList)) then 
          nPts = size(surf%pts, 2)
             
          ! Communicate the face values and the nodal values
          if (myid == 0) then 
             allocate(recvBuffer1(6*nPts), recvBuffer2(3*nPts))
             allocate(recvBuffer1d(6*nPts), recvBuffer2d(3*nPts))
          else
             allocate(recvBuffer1(0), recvBuffer2(0))
          end if

          call commUserIntegrationSurfaceVars_d(recvBuffer1, recvBuffer1d, iRho, iZippFlowGamma, surf%flowComm)
          call commUserIntegrationSurfaceVars_d(recvBuffer2, recvBuffer2d, iZippFlowX, iZippFlowZ, surf%nodeComm)

          ! *Finally* we can do the actual integrations
          if (myid == 0)  then 

             ! Allocate some temporary data needed to supply to the
             ! zipper integration routine. 
             allocate(ptValid(npts), vars(npts, nZippFlowComm), &
                  varsd(npts, nZippFlowComm), fams(size(surf%conn, 2)))

             ! Prepare for the "zipper" integration call. We have to
             ! re-order the data according to the "inv" array in each
             ! of the two comms. 
             do i=1, nPts

                ! Flow Variables
                j = surf%flowComm%inv(i)

                vars(j, iRho:iZippFlowGamma) = recvBuffer1(6*(i-1) + iRho : 6*(i-1) + iZippFlowGamma)
                varsd(j, iRho:iZippFlowGamma) = recvBuffer1d(6*(i-1) + iRho : 6*(i-1) + iZippFlowGamma)

                ! Sface is not implemented. To correctly do this,
                ! interpolate the three components of 's', do the dot
                ! product with the local normal to get the sFace value. 
                vars(j, iZippFlowSface)  = zero
                varsd(j, iZippFlowSface)  = zero

                ! Node Comm Values
                j = surf%nodeComm%inv(i)
                vars(j, iZippFlowX:iZippFlowZ) = recvBuffer2(3*i-2:3*i)
                varsd(j, iZippFlowX:iZippFlowZ) = recvBuffer2d(3*i-2:3*i)
                      
                ! The additional pt-valid array
                ptValid(j) = surf%nodeComm%valid(i)
             end do
             
             ! The family array is all the same value:
             fams = surf%famID

             ! Perform the actual integration
             call flowIntegrationZipper_d(.True., surf%conn, fams, vars, varsd, localValues, localValuesd, & 
                  famList, sps, withGathered, funcValues, funcValuesd, surf%nodeComm%Valid)
             deallocate(ptValid, vars, varsd, fams)
          end if
          deallocate(recvBuffer1, recvBuffer2, recvBuffer1d, recvBuffer2d)
       end if famInclude
    end do masterLoop
  end subroutine integrateUserSurfaces_d


#endif 

  subroutine addIntegrationSurface(pts, conn, famName, famID, nPts, nConn)
    ! Add a user-supplied integration surface. 

    use communication, only : myID
    use constants

    implicit none

    ! Input variables
    real(kind=realType), dimension(3, nPts), intent(in) :: pts
    integer(kind=intType), dimension(3, nConn), intent(in) :: conn
    integer(kind=intType), intent(in) :: nPts, nConn, famID
    character(len=*) :: famName
    type(userIntSurf), pointer :: surf

    ! Not really much to do here...we just have to save the data
    ! into the data structure untilly we actual have to do the
    ! search.
    nUserIntSurfs = nUserIntSurfs + 1
    if (nUserIntSurfs > nUserIntSurfsMax) then 
       print *,"Error: Exceeded the maximum number of user-supplied "&
            &"integration slices. Increase nUserIntSurfsMax"
       stop
    end if

    if (myid == 0) then 
       surf => userIntSurfs(nUserIntSurfs)

       allocate(surf%pts(3, nPts), surf%conn(3, nConn))
       surf%pts = pts
       surf%conn = conn
       surf%famName = famName
       surf%famID = famID
    end if

  end subroutine addIntegrationSurface

  subroutine buildVolumeADTs(oBlocks, useDual)

    ! This builds volume ADTs for the the owned blocks. It will build
    ! either the dual mesh or the primal mesh depening on the flag
    ! useDual. 

    use constants
    use overset, only : oversetBlock
    use blockPointers, only : nDom, x, ie, je, ke, il, jl, kl,  vol, ib, jb, kb,& 
         iBlank, BCData, nBocos, BCFaceID, BCType
    use adtBuild, only : buildSerialHex
    use utils, only : setPointers, EChk
    implicit none

    ! Input/Output Parameters
    type(oversetBlock), dimension(:), target :: oBlocks
    logical :: useDual

    ! Working Parameters
    integer(kind=intType) :: nInterpol, nn, i, j, k, iii, jjj, kkk
    integer(kind=intType) :: iStart, jStart, kStart, iEnd, jEnd, kEnd
    integer(kind=intType) :: ii, jj, kk, mm, nADT, nHexa, planeOffset
    type(oversetBlock), pointer :: oBlock

    nInterpol = 1 ! we get the ADT to compute the interpolated volume for us. 

    domainLoop: do nn=1, nDom

       call setPointers(nn, 1, 1)
       oBlock => oBlocks(nn)

       primalOrDual: if (useDual) then 

          ! Now setup the data for the ADT
          nHexa = il * jl * kl
          nADT = ie * je * ke
          oBlock%il = il
          oBlock%jl = jl
          oBlock%kl = kl

          allocate(oBlock%xADT(3, nADT), oBlock%hexaConn(8, nHexa), &
               oBlock%qualDonor(1, nADT))
          ! Fill up the xADT using cell centers (dual mesh)
          mm = 0
          do k=1, ke
             do j=1, je
                do i=1, ie
                   mm = mm + 1
                   oBlock%xADT(:, mm) = eighth*(&
                        x(i-1, j-1, k-1, :) + &
                        x(i  , j-1, k-1, :) + &
                        x(i-1, j  , k-1, :) + &
                        x(i  , j  , k-1, :) + &
                        x(i-1, j-1, k  , :) + &
                        x(i  , j-1, k  , :) + &
                        x(i-1, j  , k  , :) + &
                        x(i  , j  , k  , :))
                   oBlock%qualDonor(1, mm) = vol(i, j, k)
                end do
             end do
          end do
          
          do mm=1,nBocos
             ! We need to make sure that the interpolation does not
             ! use a halo behind an overset outer bound. This could
             ! happen because the last cell is interpolated (-1) and
             ! the iblank on the halos are still 1. Just set the
             ! quality very high so it is not accepted. 
          
             select case (BCFaceID(mm))
             case (iMin)
                iStart=1; iEnd=1;
                jStart=BCData(mm)%icBeg; jEnd=BCData(mm)%icEnd
                kStart=BCData(mm)%jcBeg; kEnd=BCData(mm)%jcEnd
             case (iMax)
                iStart=ie; iEnd=ie;
                jStart=BCData(mm)%icBeg; jEnd=BCData(mm)%icEnd
                kStart=BCData(mm)%jcBeg; kEnd=BCData(mm)%jcEnd
             case (jMin)
                iStart=BCData(mm)%icBeg; iEnd=BCData(mm)%icEnd
                jStart=1; jEnd=1
                kStart=BCData(mm)%jcBeg; kEnd=BCData(mm)%jcEnd
             case (jMax)
                iStart=BCData(mm)%icBeg; iEnd=BCData(mm)%icEnd
                jStart=je; jEnd=je;
                kStart=BCData(mm)%jcBeg; kEnd=BCData(mm)%jcEnd
             case (kMin)
                iStart=BCData(mm)%icBeg; iEnd=BCData(mm)%icEnd
                jStart=BCData(mm)%jcBeg; jEnd=BCData(mm)%jcEnd
                kStart=1; kEnd=1;
             case (kMax)
                iStart=BCData(mm)%icBeg; iEnd=BCData(mm)%icEnd
                jStart=BCData(mm)%jcBeg; jEnd=BCData(mm)%jcEnd
                kStart=ke; kEnd=ke;
             end select
             
             if (BCType(mm) == OversetOuterBound) then
                do k=kStart, kEnd
                   do j=jStart, jEnd
                      do i=iStart, iEnd
                         ! recompute the index
                         kk = (k-1)*ie*je + (j-1)*ie + i
                         oBlock%qualDonor(1, kk) = large
                      end do
                   end do
                end do
             end if
          end do
       

          mm = 0
          ! These are the 'elements' of the dual mesh.
          planeOffset = ie * je
          do k=2, ke
             do j=2, je
                do i=2, ie
                   mm = mm + 1
                   oBlock%hexaConn(1, mm) = (k-2)*planeOffset + (j-2)*ie + (i-2) + 1
                   oBlock%hexaConn(2, mm) = oBlock%hexaConn(1, mm) + 1 
                   oBlock%hexaConn(3, mm) = oBlock%hexaConn(2, mm) + ie
                   oBlock%hexaConn(4, mm) = oBlock%hexaConn(3, mm) - 1 

                   oBlock%hexaConn(5, mm) = oBlock%hexaConn(1, mm) + planeOffset
                   oBlock%hexaConn(6, mm) = oBlock%hexaConn(2, mm) + planeOffset
                   oBlock%hexaConn(7, mm) = oBlock%hexaConn(3, mm) + planeOffset
                   oBlock%hexaConn(8, mm) = oBlock%hexaConn(4, mm) + planeOffset
                end do
             end do
          end do
       else
          ! Note that we will be including the halo primal cells. This
          ! should slightly increase robusness for viscous off-wall
          ! spacing. This means the primal mesh has 1 MORE node/cell
          ! in each direction. 

          ! Now setup the data for the ADT
          nHexa = ie * je * ke
          nADT = ib * jb * kb
          oBlock%il = ie
          oBlock%jl = je
          oBlock%kl = ke

          allocate(oBlock%xADT(3, nADT), oBlock%hexaConn(8, nHexa), &
               oBlock%qualDonor(1, nADT))

          oBlock%qualDonor = zero
          ! Fill up the xADT using the primal nodes
          mm = 0
          do k=0, ke
             do j=0, je
                do i=0, ie
                   mm = mm + 1
                   oBlock%xADT(:, mm) = x(i, j, k, :)

                   ! Since we don't have all 8 volumes surrounding the
                   ! halo nodes, clip the volumes to be between 0 and ib etc.
                   do iii=0,1
                      do jjj=0,1
                         do kkk=0,1
                            ii = min(max(0, iii+i), ib)
                            jj = min(max(0, jjj+j), jb)
                            kk = min(max(0, kkk+k), kb)

                            oBlock%qualDonor(1, mm) = oBlock%qualDonor(1, mm) + &
                                 vol(ii, jj, kk)
                         end do
                      end do
                   end do

                   ! Dividing by 8 isn't strictly necessary but we'll
                   ! do it anyway.
                   oBlock%qualDonor(1, mm) = oBlock%qualDonor(1, mm) * eighth
                end do
             end do
          end do

          mm = 0
          ! These are the 'elements' of the dual mesh.
          planeOffset = ib * jb
          do k=1, ke
             do j=1, je
                do i=1, ie
                   mm = mm + 1
                   oBlock%hexaConn(1, mm) = (k-1)*planeOffset + (j-1)*ib + (i-1) + 1
                   oBlock%hexaConn(2, mm) = oBlock%hexaConn(1, mm) + 1 
                   oBlock%hexaConn(3, mm) = oBlock%hexaConn(2, mm) + ib
                   oBlock%hexaConn(4, mm) = oBlock%hexaConn(3, mm) - 1 

                   oBlock%hexaConn(5, mm) = oBlock%hexaConn(1, mm) + planeOffset
                   oBlock%hexaConn(6, mm) = oBlock%hexaConn(2, mm) + planeOffset
                   oBlock%hexaConn(7, mm) = oBlock%hexaConn(3, mm) + planeOffset
                   oBlock%hexaConn(8, mm) = oBlock%hexaConn(4, mm) + planeOffset
                end do
             end do
          end do
       end if primalOrDual

       ! Call the custom build routine -- Serial only, only Hexa volumes,
       ! we supply our own ADT Type

       call buildSerialHex(nHexa, nADT, oBlock%xADT, oBlock%hexaConn, oBlock%ADT)
    end do domainLoop

  end subroutine buildVolumeADTs

  subroutine performInterpolation(pts, oBlocks, useDual, comm)

    ! This routine performs the actual searches for the slices. It is
    ! generic in the sense that it will search an arbtitrary set of
    ! points on either the primal or dual meshes. The final required
    ! communication data is then written into the supplied comm. 

    use constants
    use block, only : interpPtType
    use communication, only : adflow_comm_world, myid, nProc
    use overset, only : oversetBlock
    use blockPointers, only : nDom, x, ie, je, ke, il, jl, kl, x, iBlank, vol
    use adtLocalSearch, only :  mindistancetreesearchsinglepoint, &
         containmenttreesearchsinglepoint
    use adtUtils, only : stack
    use adtData, only : adtBBoxTargetType
    use utils, only : setPointers, mynorm2, EChk
    use inputOverset, only : oversetProjTol
    use oversetUtilities, only : fracToWeights2, getCumulativeForm

    implicit none

    ! Input parameters:
    real(kind=realType), dimension(:, :), intent(in) :: pts
    type(userIntSurf) :: surf
    type(oversetBlock), dimension(:), target, intent(in) :: oBlocks
    logical, intent(in) :: useDual
    !type(oversetSurf), dimension(:), intent(in) :: oSurfs
    type(userSurfCommType) :: comm

    ! Working parameters
    type(oversetBlock), pointer :: oBlock
    type(interpPtType), dimension(:), allocatable :: surfFringes

    integer(Kind=intType) :: i, j, k, ii, jj, kk, iii, jjj, kkk, nn, mm 
    integer(kind=intType) :: iSurf, ierr, nInterpol, iProc
    integer(kind=intType) :: nHexa, nAdt, planeOffset, elemID, nPts
    real(kind=realType) :: xc(4), weight(8)
    integer mpiStatus(MPI_STATUS_SIZE) 

    real(kind=realType) :: uvw(5), uvw2(5), donorQual, xcheck(3)
    integer(kind=intType) :: intInfo(3), intInfo2(3)
    logical :: failed, invalid

    integer(kind=intType), dimension(:, :), allocatable :: donorInfo, intSend
    integer(kind=intType), dimension(:), allocatable :: procSizes
    real(kind=realType), dimension(:, :), allocatable :: donorFrac, realSend

    ! Variables we have to pass the ADT search routine
    integer(kind=intType), dimension(:), pointer :: BB
    type(adtBBoxTargetType), dimension(:), pointer :: BB2
    integer(kind=intType), dimension(:), pointer :: frontLeaves
    integer(kind=intType), dimension(:), pointer :: frontLeavesNew

    ! Data for the search
    allocate(BB(20), frontLeaves(25), frontLeavesNew(25), stack(100))

    nPts = size(pts, 2)

    ! Allocate donor information arrays
    allocate(donorFrac(4, nPts), donorInfo(5, nPts))
    donorInfo = -1
    donorFrac(4, :) = large
    nInterpol = 1
    domainSearch: do nn=1, nDom
       oBlock => oBlocks(nn)
       call setPointers(nn, 1, 1)

       ! Search the supplied pts one at a time
       elemLoop: do i=1, nPts

          xc(1:3) = pts(:, i)

          ! Call the standard tree search
          call containmentTreeSearchSinglePoint(oBlock%ADT, xc, intInfo, uvw, &
               oBlock%qualDonor, nInterpol, BB, frontLeaves, frontLeavesNew, failed)

          ! Make sure this point is not garbage.
          if (intInfo(1) >= 0) then 
             call fracToWeights2(uvw(1:3), weight)
             xcheck = zero
             do j=1,8
                xcheck = xcheck + weight(j)*oBlock%xADT(:, oBlock%hexaConn(j, intInfo(3)))
             end do

             if (mynorm2(xcheck - xc(1:3)) > oversetProjTol) then 
                failed = .True.
             end if
          end if

          if (intInfo(1) >= 0 .and. failed) then 
             ! we "found" a point but it is garbage. Do the failsafe search
             xc(4) = large
             call minDistanceTreeSearchSinglePoint(oBlock%ADT, xc, intInfo, uvw, &
                  oBlock%qualDonor, nInterpol, BB2, frontLeaves, frontLeavesNew)

             ! Check this one:
             call fracToWeights2(uvw(1:3), weight)
             xcheck = zero
             do j=1,8
                xcheck = xcheck + weight(j)*oBlock%xADT(:, oBlock%hexaConn(j, intInfo(3)))
             end do

             ! Since this is the last line of defence, relax the tolerance a bit
             if (mynorm2(xcheck - xc(1:3)) > 100*oversetProjTol) then 
                ! This fringe has not found a donor
                intInfo(1) = -1
             else
                ! This one has now passed.

                ! Important! uvw(4) is the distance squared for this search
                ! not interpolated value
                uvw(4) = uvw(5)
             end if
          end if

          elemFound: if (intInfo(1) >= 0) then 

             ! Donor and block and index information for this donor. 
             donorQual = uvw(4)
             elemID = intInfo(3) - 1 ! Make it zero based

             ! The dual mesh needs an offset of 1 becuase it only used
             ! 1:ie values. This is not necessary for the primal. 
             if (useDual) then 
                ii = mod(elemID, oBlock%il) + 1
                jj = mod(elemID/oBlock%il, oBlock%jl) + 1
                kk = elemID/(oBlock%il*oBlock%jl) + 1
             else
                ii = mod(elemID, oBlock%il) 
                jj = mod(elemID/oBlock%il, oBlock%jl) 
                kk = elemID/(oBlock%il*oBlock%jl) 
             end if
             ! Rememebr donorFrac(4, i) is the current best quality 
             if ( donorQual < donorFrac(4, i)) then 

                invalid = .False.

                ! For the dual mesh search, we have to make sure the
                ! potential donors are valid. Such a check is not
                ! necessary for the primal search since all nodes are
                ! considered valid. 

                if (useDual) then 
                   ! Check if the point is invalid. We can do this
                   ! with i-blank array. We can only accept compute
                   ! cells (iblank=1) or interpolated
                   ! cells(iblank=-1).
                   do kkk=0,1
                      do jjj=0,1
                         do iii=0,1
                            if (.not. (iblank(ii+iii, jj+jjj, kk+kkk) == 1 .or. &
                                 iblank(ii+iii, jj+jjj, kk+kkk) == -1)) then 
                               invalid = .True.
                            end if
                         end do
                      end do
                   end do
                end if

                if (.not. invalid) then 

                   ! Set the quality of the donor to the one we
                   ! just found. Save the rest of the necessary
                   ! information. 
                   donorInfo(1, i) = myid
                   donorInfo(2, i) = nn
                   donorInfo(3, i) = ii
                   donorInfo(4, i) = jj
                   donorInfo(5, i) = kk
                   donorFrac(1:3, i) = uvw(1:3)
                   donorFrac(4, i) = donorQual
                end if
             end if
          end if elemFound
       end do elemLoop
    end do domainSearch

    ! Next count up the number of valid donors we've found and compact
    ! the info back to that length.
    if (myid /=0) then 
       j = 0
       do i=1, nPts
          if (donorInfo(1, i) /= -1) then 
             j = j + 1
          end if
       end do
       allocate(intSend(6, j), realSend(4, i))
       if (j > 0) then 
          j = 0
          do i=1, nPts
             if (donorInfo(1, i) /= -1) then 
                j = j + 1
                intSend(1:5, j) = donorInfo(:, i)
                intSend(6, j) = i
                realSend(:, j) = donorFrac(:, i)
             end if
          end do
       end if
    else
       ! On the root proc, use intSend and realSend as the receiver
       ! buffer. These can be at most nPts sized.
       allocate(intSend(6, nPts), realSend(4, nPts))
    end if

    ! Gather up the sizes (j) to the root processor so he know who to
    ! expect data from.
    allocate(procSizes(0:nProc-1))

    call mpi_gather(j, 1, adflow_integer, procSizes, 1, &
         adflow_integer, 0, adflow_comm_world, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Next all the procs need to send all the information back to the
    ! root processor where we will determine the proper donors for
    ! each of the cells

    ! All procs except root fire off their data.
    if (myid >= 1) then 
       if (j > 0) then 
          call mpi_send(intSend, j*6, adflow_integer, 0, myid, &
               adflow_comm_world, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          call mpi_send(realSend, j*4, adflow_real, 0, myid, &
               adflow_comm_world, ierr)
          call EChk(ierr,__FILE__,__LINE__)
       end if
    end if

    ! And the root processor recieves it...
    if (myid == 0) then 
       do iProc=1, nProc-1
          ! Determine if this proc has sent anything:
          if (procSizes(iProc) /= 0) then 

             call MPI_recv(intSend, 6*nPts, adflow_integer, iProc, iProc,&
                  adflow_comm_world, mpiStatus, ierr)
             call EChk(ierr,__FILE__,__LINE__)

             call MPI_recv(realSend, 4*nPts, adflow_real, iProc, iProc,&
                  adflow_comm_world, mpiStatus, ierr)
             call EChk(ierr,__FILE__,__LINE__)

             ! Now process the data (intSend and realSend) that we
             ! just received. We don't need to check the status for
             ! the sizes becuase we already know the sizes from the
             ! initial gather we did. 

             do i=1, procSizes(iProc)
                ii = intSend(6, i)

                if (realSend(4, i) < donorFrac(4, ii)) then 
                   ! The incoming quality is better. Accept it. 
                   donorInfo(1:5, ii) = intSend(1:5, i)
                   donorFrac(:, ii) = realSend(:, i)
                end if
             end do
          end if
       end do

       ! To make this easier, convert the information we have to a
       ! 'fringeType' array so we can use the pre-existing sorting
       ! routine.
       allocate(surfFringes(nPts))
       do i=1, nPts
          surfFringes(i)%donorProc = donorInfo(1, i)
          surfFringes(i)%donorBlock = donorInfo(2, i)
          surfFringes(i)%dI = donorInfo(3, i)
          surfFringes(i)%dJ = donorInfo(4, i)
          surfFringes(i)%dK = donorInfo(5, i)
          surfFringes(i)%donorFrac = donorFrac(1:3, i)
          ! Use the myBlock attribute to keep track of the original
          ! index. When we sort the fringes, they will no longer be
          ! in the same order
          surfFringes(i)%myBlock = i
       end do

       ! Perform the actual sort. 
       call qsortInterpPtType(surfFringes, nPts)

       ! We will reuse-proc sizes to now mean the number of elements
       ! that the processor *actually* has to send. We will include
       ! the root proc itself in the calc becuase that will tell us
       ! the size of the internal comm structure. 

       procSizes = 0
       allocate(comm%valid(nPts))
       comm%valid = .True.
       do i=1, nPts
          if (surfFringes(i)%donorProc < 0) then 
             ! We dont have a donor. Flag this point as invalid
             comm%valid(i) = .False.
          end if

          ! Dump the points without donors on the root proc by making
          ! sure j is at least 0 for the root proc. These will just
          ! simply be ignored during the comm. 
          j = max(surfFringes(i)%donorProc, 0)
          procSizes(j) = procSizes(j) + 1
       end do
    end if

    ! Simply broadcast out the the proc sizes back to everyone so all
    ! processors know if they are to receive anything back. 
    call mpi_bcast(procSizes, nProc, adflow_Integer, 0, adflow_comm_world, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! We can now save some of the final comm data required on the
    ! root proc for this surf. 
    allocate(comm%procSizes(0:nProc-1), comm%procDisps(0:nProc))

    ! Copy over procSizes and generate the cumulative form of
    ! the size array, procDisps
    comm%procSizes = procSizes
    call getCumulativeForm(comm%procSizes, nProc, comm%procDisps)

    ! Record the elemInverse which is necessary to index into
    ! the original conn array. 
    if (myid == 0) then 
       allocate(comm%Inv(nPts))
       do i=1, nPts
          comm%Inv(i) = surfFringes(i)%myBlock
       end do
    end if

    ! Now we can send out the final donor information to the
    ! processors that must supply it. 
    comm%nDonor = procSizes(myID)
    allocate(comm%frac(3, comm%nDonor), comm%donorInfo(4, comm%nDonor))

    if (myid >= 1) then 
       if (comm%nDonor > 0) then 
          ! We are responible for at least 1 donor. We have to make
          ! use of the intSend and realSend buffers again (which
          ! are guaranteed to be big enough). The reason we can't
          ! dump the data in directlyis that intSend and realSend
          ! have a different leading index than we need on the
          ! final data structure. 

          call MPI_recv(intSend, 6*comm%nDonor, adflow_integer, 0, myid, &
               adflow_comm_world, mpiStatus, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          call MPI_recv(realSend, 4*comm%nDonor, adflow_real, 0, myID, &
               adflow_comm_world, mpiStatus, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! Copy into final structure
          do i=1, comm%nDonor
             comm%donorInfo(:, i) = intSend(1:4, i)
             comm%frac(:, i) = realSend(1:3, i)
          end do
       end if
    else
       ! We are the root processor.
       if (comm%nDonor > 0) then 
          ! We need to copy out our donor info on the root proc if we have any
          do i= comm%procDisps(myID)+1, comm%procDisps(myID+1)
             comm%donorInfo(1, i) = surfFringes(i)%donorBlock
             comm%donorInfo(2, i) = surfFringes(i)%dI
             comm%donorInfo(3, i) = surfFringes(i)%dJ
             comm%donorInfo(4, i) = surfFringes(i)%dK
             comm%frac(1:3, i)    = surfFringes(i)%donorFrac
          end do
       end if

       ! Now loop over the rest of the procs and send out the info we
       ! need. We have to temporarily copy the data back out of
       ! fringes to the intSend and realSend arrays
       do iProc=1, nProc-1

          if (comm%procSizes(iProc) > 0) then 
             ! Have something to send here:
             j = 0
             do i=comm%procDisps(iProc)+1, comm%procDisps(iProc+1)
                j = j + 1

                intSend(1, j) = surfFringes(i)%donorBlock
                intSend(2, j) = surfFringes(i)%dI      
                intSend(3, j) = surfFringes(i)%dJ
                intSend(4, j) = surfFringes(i)%dK
                realSend(1:3, j) = surfFringes(i)%donorFrac
             end do

             call mpi_send(intSend, j*6, adflow_integer, iProc, iProc, &
                  adflow_comm_world, ierr)
             call EChk(ierr,__FILE__,__LINE__)

             call mpi_send(realSend, j*4, adflow_real, iProc, iProc, &
                  adflow_comm_world, ierr)
             call EChk(ierr,__FILE__,__LINE__)
          end if
       end do

       ! Deallocate data allocatd only on root proc
       deallocate(surfFringes)
    end if

    ! Nuke rest of allocated on all procs
    deallocate(intSend, realSend, procSizes, donorInfo, donorFrac)

  end subroutine performInterpolation

  subroutine interpolateIntegrationSurfaces

    ! This routine performs the actual searches for the slices. We
    ! reuse much of the same machinery as is used in the overset code.

    use constants
    use communication, only : adflow_comm_world, myid
    use overset, only : oversetBlock
    use blockPointers, only : nDom, ie, je, ke, il, jl, kl
    use adtBuild, only : buildSerialHex, destroySerialHex
    use utils, only : setPointers, EChk

    implicit none

    ! Working parameters
    type(oversetBlock), dimension(nDom), target :: oBlocks
    type(userIntSurf), pointer :: surf

    integer(Kind=intType) :: iSurf, ii, i, nn, nPts, ierr
    real(kind=realType), dimension(:, :), allocatable :: pts
    logical :: useDual
    if (nUserIntSurfs == 0) then 
       return
    end if
    primalDualLoop: do ii=1, 2
       if (ii==1) then 
          useDual = .True.
       else
          useDual = .False.
       end if

       call buildVolumeADTs(oBlocks, useDual)

       masterLoop: do iSurf=1, nUserIntSurfs

          surf => userIntSurfs(iSurf)

          if (myid == 0) then 

             ! We are interpolating the nodal values for both the
             ! nodes and the solution variables. 
             nPts = size(surf%pts, 2)
             allocate(pts(3, nPts))
             nodeLoop: do i=1, nPts
                pts(:, i) = surf%pts(:, i)
             end do nodeLoop
          end if

          ! Send the number of points back to all procs:
          call mpi_bcast(nPts, 1, adflow_integer, 0, adflow_comm_world, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! All other procs except the root allocate space and receive
          ! the pt array. 
          if (myid /= 0) then 
             allocate(pts(3, nPts))
          end if

          call mpi_bcast(pts, 3*nPts, adflow_real, 0, adflow_comm_world, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! Call the actual interpolation routine
          if (ii==1) then 
             call performInterpolation(pts, oBlocks, .True., surf%flowComm)
          else
             call performInterpolation(pts, oBlocks, .False., surf%nodeComm)
          end if

          deallocate(pts)
       end do masterLoop

       ! Destroy the ADT Data and allocated values
       do nn=1, nDom
          call destroySerialHex(oBlocks(nn)%ADT)
          deallocate(oBlocks(nn)%xADT, oBlocks(nn)%hexaConn, oBlocks(nn)%qualDonor)
       end do
    end do primalDualLoop
  end subroutine interpolateIntegrationSurfaces

  subroutine commUserIntegrationSurfaceVars(recvBuffer, varStart, varEnd, comm)

    use constants
    use block, onlY : flowDoms, nDom
    use communication, only : myid, adflow_comm_world
    use utils, only : EChk
    use oversetUtilities, only :fracToWeights

    implicit none

    ! Input/Output
    real(kind=realType), dimension(:) :: recvBuffer
    integer(kind=intType), intent(in) :: varStart, varEnd
    type(userSurfCommType) :: comm

    ! Working
    real(kind=realType), dimension(:), allocatable :: sendBuffer
    integer(Kind=intType) :: d1, i1, j1, k1, jj, k, nvar, i, ierr
    real(kind=realType), dimension(8) :: weight

    ! The number of variables we are transferring:
    nVar = varEnd - varStart + 1

    ! We assume that the pointers to the realCommVars have already been set. 

    allocate(sendBuffer(nVar*comm%nDonor))

    ! First generate the interpolated data necessary
    jj = 0
    donorLoop: do i=1, comm%nDonor
       ! Convert the frac to weights
       call fracToWeights(comm%frac(:, i), weight)

       ! Block and indices for easier reading. The +1 is due to the
       ! pointer offset on realCommVars.

       d1 = comm%donorInfo(1, i) ! Block Index
       i1 = comm%donorInfo(2, i)+1 ! donor I index
       j1 = comm%donorInfo(3, i)+1 ! donor J index
       k1 = comm%donorInfo(4, i)+1 ! donor K index

       ! We are interpolating nVar variables
       do k=varStart, varEnd
          jj = jj + 1
          if (d1 > 0) then ! Is this pt valid?
             sendBuffer(jj) = &
                  weight(1)*flowDoms(d1,1,1)%realCommVars(k)%var(i1  ,j1  ,k1  ) + &
                  weight(2)*flowDoms(d1,1,1)%realCommVars(k)%var(i1+1,j1  ,k1  ) + &
                  weight(3)*flowDoms(d1,1,1)%realCommVars(k)%var(i1  ,j1+1,k1  ) + &
                  weight(4)*flowDoms(d1,1,1)%realCommVars(k)%var(i1+1,j1+1,k1  ) + &
                  weight(5)*flowDoms(d1,1,1)%realCommVars(k)%var(i1  ,j1  ,k1+1) + &
                  weight(6)*flowDoms(d1,1,1)%realCommVars(k)%var(i1+1,j1  ,k1+1) + &
                  weight(7)*flowDoms(d1,1,1)%realCommVars(k)%var(i1  ,j1+1,k1+1) + &
                  weight(8)*flowDoms(d1,1,1)%realCommVars(k)%var(i1+1,j1+1,k1+1)
          end if
       end do
    end do donorLoop

    ! Now we can do an mpi_gatherv to the root proc:
    call mpi_gatherv(sendBuffer, nVar*comm%nDonor, adflow_real, recvBuffer, &
         nVar*comm%procSizes, nVar*comm%procDisps, adflow_real, 0, adflow_comm_world, ierr)
    call EChk(ierr,__FILE__,__LINE__)
    deallocate(sendBuffer)

  end subroutine commUserIntegrationSurfaceVars

 subroutine commUserIntegrationSurfaceVars_d(recvBuffer, recvBufferd, varStart, varEnd, comm)

    use constants
    use block, onlY : flowDoms, nDom
    use communication, only : myid, adflow_comm_world
    use utils, only : EChk
    use oversetUtilities, only :fracToWeights

    implicit none

    ! Input/Output
    real(kind=realType), dimension(:) :: recvBuffer, recvBufferd
    integer(kind=intType), intent(in) :: varStart, varEnd
    type(userSurfCommType) :: comm

    ! Working
    real(kind=realType), dimension(:), allocatable :: sendBuffer, sendBufferd
    integer(Kind=intType) :: d1, i1, j1, k1, jj, k, nvar, i, ierr
    real(kind=realType), dimension(8) :: weight

    ! The number of variables we are transferring:
    nVar = varEnd - varStart + 1

    ! We assume that the pointers to the realCommVars have already been set. 

    allocate(sendBuffer(nVar*comm%nDonor), &
         sendBufferd(nVar*comm%nDonor))

    ! First generate the interpolated data necessary
    jj = 0
    donorLoop: do i=1, comm%nDonor
       ! Convert the frac to weights
       call fracToWeights(comm%frac(:, i), weight)

       ! Block and indices for easier reading. The +1 is due to the
       ! pointer offset on realCommVars.

       d1 = comm%donorInfo(1, i) ! Block Index
       i1 = comm%donorInfo(2, i)+1 ! donor I index
       j1 = comm%donorInfo(3, i)+1 ! donor J index
       k1 = comm%donorInfo(4, i)+1 ! donor K index

       ! We are interpolating nVar variables
       do k=varStart, varEnd
          jj = jj + 1
          if (d1 > 0) then ! Is this pt valid?
             sendBuffer(jj) = &
                  weight(1)*flowDoms(d1,1,1)%realCommVars(k)%var(i1  ,j1  ,k1  ) + &
                  weight(2)*flowDoms(d1,1,1)%realCommVars(k)%var(i1+1,j1  ,k1  ) + &
                  weight(3)*flowDoms(d1,1,1)%realCommVars(k)%var(i1  ,j1+1,k1  ) + &
                  weight(4)*flowDoms(d1,1,1)%realCommVars(k)%var(i1+1,j1+1,k1  ) + &
                  weight(5)*flowDoms(d1,1,1)%realCommVars(k)%var(i1  ,j1  ,k1+1) + &
                  weight(6)*flowDoms(d1,1,1)%realCommVars(k)%var(i1+1,j1  ,k1+1) + &
                  weight(7)*flowDoms(d1,1,1)%realCommVars(k)%var(i1  ,j1+1,k1+1) + &
                  weight(8)*flowDoms(d1,1,1)%realCommVars(k)%var(i1+1,j1+1,k1+1)
             sendBufferd(jj) = &
                  weight(1)*flowDoms(d1,1,1)%realCommVars(k+varEnd)%var(i1  ,j1  ,k1  ) + &
                  weight(2)*flowDoms(d1,1,1)%realCommVars(k+varEnd)%var(i1+1,j1  ,k1  ) + &
                  weight(3)*flowDoms(d1,1,1)%realCommVars(k+varEnd)%var(i1  ,j1+1,k1  ) + &
                  weight(4)*flowDoms(d1,1,1)%realCommVars(k+varEnd)%var(i1+1,j1+1,k1  ) + &
                  weight(5)*flowDoms(d1,1,1)%realCommVars(k+varEnd)%var(i1  ,j1  ,k1+1) + &
                  weight(6)*flowDoms(d1,1,1)%realCommVars(k+varEnd)%var(i1+1,j1  ,k1+1) + &
                  weight(7)*flowDoms(d1,1,1)%realCommVars(k+varEnd)%var(i1  ,j1+1,k1+1) + &
                  weight(8)*flowDoms(d1,1,1)%realCommVars(k+varEnd)%var(i1+1,j1+1,k1+1)
          end if
       end do
    end do donorLoop

    ! Now we can do an mpi_gatherv to the root proc:
    call mpi_gatherv(sendBuffer, nVar*comm%nDonor, adflow_real, recvBuffer, &
         nVar*comm%procSizes, nVar*comm%procDisps, adflow_real, 0, adflow_comm_world, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call mpi_gatherv(sendBufferd, nVar*comm%nDonor, adflow_real, recvBufferd, &
         nVar*comm%procSizes, nVar*comm%procDisps, adflow_real, 0, adflow_comm_world, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    deallocate(sendBuffer, sendBufferd)

  end subroutine commUserIntegrationSurfaceVars_d

subroutine commUserIntegrationSurfaceVars_b(recvBuffer, recvBufferd, varStart, varEnd, comm)

    use constants
    use block, onlY : flowDoms, nDom
    use communication, only : myid, adflow_comm_world
    use utils, only : EChk
    use oversetUtilities, only :fracToWeights

    implicit none

    ! Input/Output
    real(kind=realType), dimension(:) :: recvBuffer, recvBufferd
    integer(kind=intType), intent(in) :: varStart, varEnd
    type(userSurfCommType) :: comm

    ! Working
    real(kind=realType), dimension(:), allocatable :: sendBuffer, sendBufferd
    integer(Kind=intType) :: d1, i1, j1, k1, jj, k, nvar, i, ierr
    real(kind=realType), dimension(8) :: weight

    ! The number of variables we are transferring:
    nVar = varEnd - varStart + 1

    allocate(sendBufferd(nVar*comm%nDonor))

    ! Adjoint of a gatherv is a scatterv. Flip the send/recv relative
    ! to the forward call. 
    call mpi_scatterv(recvBufferd, nVar*comm%procSizes, nVar*comm%procDisps, adflow_real, &
         sendBufferd, nVar*comm%nDonor, adflow_real, 0, adflow_comm_world, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! First generate the interpolated data necessary
    jj = 0
    donorLoop: do i=1, comm%nDonor
       ! Convert the frac to weights
       call fracToWeights(comm%frac(:, i), weight)

       ! Block and indices for easier reading. The +1 is due to the
       ! pointer offset on realCommVars.

       d1 = comm%donorInfo(1, i) ! Block Index
       i1 = comm%donorInfo(2, i)+1 ! donor I index
       j1 = comm%donorInfo(3, i)+1 ! donor J index
       k1 = comm%donorInfo(4, i)+1 ! donor K index

       ! We are interpolating nVar variables
       do k=varStart, varEnd
          jj = jj + 1
          if (d1 > 0) then ! Is this pt valid?

             ! Accumulate back onto the derivative variables. 
             flowDoms(d1,1,1)%realCommVars(k+nVar)%var(i1  , j1  , k1) = & 
                  flowDoms(d1,1,1)%realCommVars(k+nVar)%var(i1  , j1  , k1) + weight(1)*sendBufferd(jj)

             flowDoms(d1,1,1)%realCommVars(k+nVar)%var(i1+1, j1  , k1) = & 
                  flowDoms(d1,1,1)%realCommVars(k+nVar)%var(i1+1, j1  , k1) + weight(2)*sendBufferd(jj)

             flowDoms(d1,1,1)%realCommVars(k+nVar)%var(i1  , j1+1, k1) = & 
                  flowDoms(d1,1,1)%realCommVars(k+nVar)%var(i1  , j1+1, k1) + weight(3)*sendBufferd(jj)

             flowDoms(d1,1,1)%realCommVars(k+nVar)%var(i1+1, j1+1, k1) = & 
                  flowDoms(d1,1,1)%realCommVars(k+nVar)%var(i1+1, j1+1, k1) + weight(4)*sendBufferd(jj)

             flowDoms(d1,1,1)%realCommVars(k+nVar)%var(i1  , j1  , k1+1) = & 
                  flowDoms(d1,1,1)%realCommVars(k+nVar)%var(i1  , j1  , k1+1) + weight(5)*sendBufferd(jj)

             flowDoms(d1,1,1)%realCommVars(k+nVar)%var(i1+1, j1  , k1+1) = & 
                  flowDoms(d1,1,1)%realCommVars(k+nVar)%var(i1+1, j1  , k1+1) + weight(6)*sendBufferd(jj)

             flowDoms(d1,1,1)%realCommVars(k+nVar)%var(i1  , j1+1, k1+1) = & 
                  flowDoms(d1,1,1)%realCommVars(k+nVar)%var(i1  , j1+1, k1+1) + weight(7)*sendBufferd(jj)

             flowDoms(d1,1,1)%realCommVars(k+nVar)%var(i1+1, j1+1, k1+1) = & 
                  flowDoms(d1,1,1)%realCommVars(k+nVar)%var(i1+1, j1+1, k1+1) + weight(8)*sendBufferd(jj)
          end if
       end do
    end do donorLoop
   
    deallocate(sendBuffer, sendBufferd)

  end subroutine commUserIntegrationSurfaceVars_b



  subroutine qsortInterpPtType(arr, nn)

    use constants
    use block ! Cannot use-only becuase of <= operator
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nn

    type(interpPtType), dimension(*), intent(inout) :: arr
    !
    !      Local variables.
    !
    integer(kind=intType), parameter :: m = 7

    integer(kind=intType) :: nStack
    integer(kind=intType) :: i, j, k, r, l, jStack, ii

    integer :: ierr

    type(interpPtType) :: a, tmp

    integer(kind=intType), allocatable, dimension(:) :: stack
    integer(kind=intType), allocatable, dimension(:) :: tmpStack

    ! Allocate the memory for stack.

    nStack = 100
    allocate(stack(nStack), stat=ierr)
    if(ierr /= 0)                         &
         call terminate("qsortinterpPt", &
         "Memory allocation failure for stack")

    ! Initialize the variables that control the sorting.

    jStack = 0
    l      = 1
    r      = nn

    ! Start of the algorithm

    do

       ! Check for the size of the subarray.

       if((r-l) < m) then

          ! Perform insertion sort

          do j=l+1,r
             a = arr(j)
             do i=(j-1),l,-1
                if(arr(i) <= a) exit
                arr(i+1) = arr(i)
             enddo
             arr(i+1) = a
          enddo

          ! In case there are no more elements on the stack, exit from
          ! the outermost do-loop. Algorithm has finished.

          if(jStack == 0) exit

          ! Pop stack and begin a new round of partitioning.

          r = stack(jStack)
          l = stack(jStack-1)
          jStack = jStack - 2

       else

          ! Subarray is larger than the threshold for a linear sort.
          ! Choose median of left, center and right elements as
          ! partitioning element a.
          ! Also rearrange so that (l) <= (l+1) <= (r).

          k = (l+r)/2
          tmp      = arr(k)             ! Swap the elements
          arr(k)   = arr(l+1)           ! k and l+1.
          arr(l+1) = tmp

          if(arr(r) < arr(l)) then
             tmp    = arr(l)             ! Swap the elements
             arr(l) = arr(r)             ! r and l.
             arr(r) = tmp
          endif

          if(arr(r) < arr(l+1)) then
             tmp      = arr(l+1)         ! Swap the elements
             arr(l+1) = arr(r)           ! r and l+1.
             arr(r)   = tmp
          endif

          if(arr(l+1) < arr(l)) then
             tmp      = arr(l+1)         ! Swap the elements
             arr(l+1) = arr(l)           ! l and l+1.
             arr(l)   = tmp
          endif

          ! Initialize the pointers for partitioning.

          i = l+1
          j = r
          a = arr(l+1)

          ! The innermost loop

          do

             ! Scan up to find element >= a.
             do
                i = i+1
                if(a <= arr(i)) exit
             enddo

             ! Scan down to find element <= a.
             do
                j = j-1
                if(arr(j) <= a) exit
             enddo

             ! Exit the loop in case the pointers i and j crossed.

             if(j < i) exit

             ! Swap the element i and j.

             tmp    = arr(i)
             arr(i) = arr(j)
             arr(j) = tmp
          enddo

          ! Swap the entries j and l+1. Remember that a equals
          ! arr(l+1).

          arr(l+1) = arr(j)
          arr(j)   = a

          ! Push pointers to larger subarray on stack,
          ! process smaller subarray immediately.

          jStack = jStack + 2
          if(jStack > nStack) then

             ! Storage of the stack is too small. Reallocate.

             allocate(tmpStack(nStack), stat=ierr)
             if(ierr /= 0)                         &
                  call terminate("qsortinterpPt", &
                  "Memory allocation error for tmpStack")
             tmpStack = stack

             ! Free the memory of stack, store the old value of nStack
             ! in tmp and increase nStack.

             deallocate(stack, stat=ierr)
             if(ierr /= 0)                         &
                  call terminate("qsortinterpPt", &
                  "Deallocation error for stack")
             ii = nStack
             nStack = nStack + 100

             ! Allocate the memory for stack and copy the old values
             ! from tmpStack.

             allocate(stack(nStack), stat=ierr)
             if(ierr /= 0)                         &
                  call terminate("qsortinterpPt", &
                  "Memory reallocation error for stack")
             stack(1:ii) = tmpStack(1:ii)

             ! And finally release the memory of tmpStack.

             deallocate(tmpStack, stat=ierr)
             if(ierr /= 0)                         &
                  call terminate("qsortinterpPt", &
                  "Deallocation error for tmpStack")
          endif

          if((r-i+1) >= (j-l)) then
             stack(jStack)   = r
             r               = j-1
             stack(jStack-1) = j
          else
             stack(jStack)   = j-1
             stack(jStack-1) = l
             l               = j
          endif

       endif
    enddo

    ! Release the memory of stack.

    deallocate(stack, stat=ierr)
    if(ierr /= 0)                         &
         call terminate("qsortinterpPt", &
         "Deallocation error for stack")

    ! Check in debug mode whether the array is really sorted.

    if( debug ) then
       do i=1,(nn-1)
          if(arr(i+1) < arr(i))                 &
               call terminate("qsortinterpPt", &
               "Array is not sorted correctly")
       enddo
    endif

  end subroutine qsortInterpPtType

end module userSurfaceIntegrations
