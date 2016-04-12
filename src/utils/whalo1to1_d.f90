subroutine whalo1to1_d(level, start, end, commPressure,       &
                            commVarGamma, commLamVis, commEddyVis, &
                            commPattern, internal)
!
!      ******************************************************************
!      *                                                                *
!      * whalo1to1 exchanges the 1 to 1 internal halo's derivatives     *
!      *                                                                *
!      ******************************************************************
!
       use block
       use communication
       use inputTimeSpectral
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, start, end
       logical, intent(in) :: commPressure, commVarGamma
       logical, intent(in) :: commLamVis, commEddyVis

       type(commType), dimension(*), intent(in)         :: commPattern
       type(internalCommType), dimension(*), intent(in) :: internal

       integer(kind=intType) :: nVar, nn, k, sps

       spectralModes: do sps=1,nTimeIntervalsSpectral

          ! Set the pointers for the required variables
          domainLoop:do nn=1, nDom
             nVar = 0

             do k=start, end
                nVar = nVar + 1 
                flowDoms(nn, level, sps)%realCommVars(nVar)%var => &
                     flowDomsd(nn, level, sps)%w(:, :, :, k)
             end do

             if( commPressure )  then 
                nVar = nVar + 1
                flowDoms(nn, level, sps)%realCommVars(nVar)%var => &
                     flowDomsd(nn, level, sps)%P(:, :, :)
             end if

             if( commVarGamma ) then 
                nVar = nVar + 1
                flowDoms(nn, level, sps)%realCommVars(nVar)%var => &
                     flowDomsd(nn, 1, sps)%gamma(:, :, :)
             end if

             if( commLamVis ) then 
                nVar = nVar + 1
                flowDoms(nn, level, sps)%realCommVars(nVar)%var => &
                     flowDomsd(nn, 1, sps)%rlv(:, :, :)
             end if

             if( commEddyVis ) then 
                nVar = nVar + 1
                flowDoms(nn, level, sps)%realCommVars(nVar)%var => &
                     flowDomsd(nn, level, sps)%rev(:, :, :)
             end if

             if(nVar == 0) return
          end do domainLoop

          ! Run the generic exchange
          call wHalo1to1RealGeneric(nVar, level, sps, commPattern, internal)

       end do spectralModes

     end subroutine whalo1to1_d
     
