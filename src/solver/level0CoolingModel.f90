!
!      ******************************************************************
!      *                                                                *
!      * File:          level0CoolingModel.f90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-26-2005                                      *
!      * Last modified: 08-15-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine level0CoolingModel
!
!      ******************************************************************
!      *                                                                *
!      * The routine level0CoolingModel computes the source terms       *
!      * corresponding to the so called level 0 cooling model for       *
!      * turbines. It is a very simple model which injects uniformly    *
!      * the mass, momentum and energy source terms in a grid plane     *
!      * choosen by the user.                                           *
!      *                                                                *
!      * This model has been developed by Pratt and Whitney and should  *
!      * not be given to third parties. This implementation assumes     *
!      * the x-direction is the axial direction.                        *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use block
       use cgnsGrid
       use coolingModelLevel0
       use flowVarRefState
       use inputTimeSpectral
       use iteration
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: nn, mm, i, j, sps, nbl, nbg
       integer(kind=intType) :: indSol, indNorm, indX1, indX2

       real(kind=realType) :: dMadd, dPadd, dTadd, massFlowWheel
       real(kind=realType) :: wx, wy, wz, rx, ry, rz, dx, dy, dz
       real(kind=realType) :: velrx, velry, velrz
       real(kind=realType) :: rho, u, v, w, E, p
       real(kind=realType) :: flwn, advtd, adhod, fitm

       real(kind=realType), dimension(:,:),   pointer :: pp, gamma
       real(kind=realType), dimension(:,:,:), pointer :: xx1, xx2
       real(kind=realType), dimension(:,:,:), pointer :: ss, dww, ww
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if no cooling model is used.

       if(nPlanesLevel0CoolingModel == 0) return

       ! Loop over the number of spectral solutions.

       spectralLoop: do sps=1,nTimeIntervalsSpectral

         ! Compute the massflow of the entire wheel. The last argument
         ! indicates that the computeInletMassFlowFullWheel should
         ! perform the allreduce such that the global quantity is
         ! computed.

         call computeInletMassFlowFullWheel(sps, massFlowWheel, .true.)

         ! Loop over the number of cooling planes.

         planeLoop: do nn=1,nPlanesLevel0CoolingModel
!
!          **************************************************************
!          *                                                            *
!          * Step 1: Determine the mass flow, total pressure loss and   *
!          *         total temperature loss coefficients for this       *
!          *         cooling plane.                                     *
!          *                                                            *
!          **************************************************************
!
           dMadd = level0Cooling(nn,currentLevel)%mDotRatio &
                 * massFlowWheel/level0Cooling(nn,currentLevel)%area
           dPadd = log(one + level0Cooling(nn,currentLevel)%dpLog)
           dTadd = log(one + level0Cooling(nn,currentLevel)%dTLog)
!
!          **************************************************************
!          *                                                            *
!          * Step 2: Loop over the local number of subfaces that        *
!          *         contribute to this cooling plane and add the       *
!          *         appropriate source terms to the downstream cells.  *
!          *                                                            *
!          **************************************************************
!
           subfaceLoop: do mm=1,level0Cooling(nn,currentLevel)%nSubfaces

             ! Determine the rotational speed and center of rotation
             ! of the block to which this subface belongs.

             nbl = level0Cooling(nn,currentLevel)%blockID(mm)
             nbg = flowDoms(nbl,currentLevel,sps)%cgnsBlockID

             wx = timeRef*cgnsDoms(nbg)%rotRate(1)
             wy = timeRef*cgnsDoms(nbg)%rotRate(2)
             wz = timeRef*cgnsDoms(nbg)%rotRate(3)

             rx = cgnsDoms(nbg)%rotCenter(1)
             ry = cgnsDoms(nbg)%rotCenter(2)
             rz = cgnsDoms(nbg)%rotCenter(3)

             ! Set the pointers for the solution, pressure, residual
             ! and face normals, depending on the situation. Note that
             ! iMin, iMax, etc. are a bit abused here, because these
             ! planes may lie in the interior of the blocks.
             ! Also note that gamma and dw are only allocated on the
             ! fine grid and therefore the pointer must be set to the
             ! fine grid arrays.

             indSol  = level0Cooling(nn,currentLevel)%indSol(mm)
             indNorm = level0Cooling(nn,currentLevel)%indNorm(mm)
             indX1   = level0Cooling(nn,currentLevel)%indX1(mm)
             indX2   = level0Cooling(nn,currentLevel)%indX2(mm)

             select case (level0Cooling(nn,currentLevel)%indexDir(mm))
               case (iMin, iMax)
                 ww    => flowDoms(nbl,currentLevel,sps)%w(indSol,1:,1:,:)
                 pp    => flowDoms(nbl,currentLevel,sps)%p(indSol,1:,1:)
                 gamma => flowDoms(nbl,           1,sps)%gamma(indSol,1:,1:)
                 dww   => flowDoms(nbl,           1,sps)%dw(indSol,1:,1:,:)
                 ss    => flowDoms(nbl,currentLevel,sps)%sI(indNorm,:,:,:)
                 xx1   => flowDoms(nbl,currentLevel,sps)%x(indX1,1:,1:,:)
                 xx2   => flowDoms(nbl,currentLevel,sps)%x(indX2,1:,1:,:)

               !=========================================================

               case (jMin, jMax)
                 ww    => flowDoms(nbl,currentLevel,sps)%w(1:,indSol,1:,:)
                 pp    => flowDoms(nbl,currentLevel,sps)%p(1:,indSol,1:)
                 gamma => flowDoms(nbl,           1,sps)%gamma(1:,indSol,1:)
                 dww   => flowDoms(nbl,           1,sps)%dw(1:,indSol,1:,:)
                 ss    => flowDoms(nbl,currentLevel,sps)%sJ(:,indNorm,:,:)
                 xx1   => flowDoms(nbl,currentLevel,sps)%x(1:,indX1,1:,:)
                 xx2   => flowDoms(nbl,currentLevel,sps)%x(1:,indX2,1:,:)

               !=========================================================

               case (kMin, kMax)
                 ww    => flowDoms(nbl,currentLevel,sps)%w(1:,1:,indSol,:)
                 pp    => flowDoms(nbl,currentLevel,sps)%p(1:,1:,indSol)
                 gamma => flowDoms(nbl,           1,sps)%gamma(1:,1:,indSol)
                 dww   => flowDoms(nbl,           1,sps)%dw(1:,1:,indSol,:)
                 ss    => flowDoms(nbl,currentLevel,sps)%sJ(:,:,indNorm,:)
                 xx1   => flowDoms(nbl,currentLevel,sps)%x(1:,1:,indX1,:)
                 xx2   => flowDoms(nbl,currentLevel,sps)%x(1:,1:,indX2,:)

             end select

             ! Loop over the range of cells to which the source terms
             ! must be added.

             do j=level0Cooling(nn,currentLevel)%jcBeg(mm), &
                  level0Cooling(nn,currentLevel)%jcEnd(mm)
               do i=level0Cooling(nn,currentLevel)%icBeg(mm), &
                    level0Cooling(nn,currentLevel)%icEnd(mm)

                 ! Determine the coordinates relative to the rotation
                 ! center.

                 dx = eighth*(xx1(i-1,j-1,1) + xx1(i,j-1,1) &
                    +         xx1(i-1,j,  1) + xx1(i,j,  1) &
                    +         xx2(i-1,j-1,1) + xx2(i,j-1,1) &
                    +         xx2(i-1,j,  1) + xx2(i,j,  1)) - rx
                 dy = eighth*(xx1(i-1,j-1,2) + xx1(i,j-1,2) &
                    +         xx1(i-1,j,  2) + xx1(i,j,  2) &
                    +         xx2(i-1,j-1,2) + xx2(i,j-1,2) &
                    +         xx2(i-1,j,  2) + xx2(i,j,  2)) - ry
                 dz = eighth*(xx1(i-1,j-1,3) + xx1(i,j-1,3) &
                    +         xx1(i-1,j,  3) + xx1(i,j,  3) &
                    +         xx2(i-1,j-1,3) + xx2(i,j-1,3) &
                    +         xx2(i-1,j,  3) + xx2(i,j,  3)) - rz

                 ! Determine the variables in the rotating frame.
                 ! The static variables are the same as in the intertial
                 ! frame, but the velocity components as well as the
                 ! energy must be must be corrected, i.e. the wheel
                 ! speed must be substracted.

                 velrx = wy*rz - wz*ry
                 velry = wz*rx - wx*rz
                 velrz = wx*ry - wy*rx

                 rho = ww(i,j,irho)
                 u   = ww(i,j,ivx) - velrx
                 v   = ww(i,j,ivy) - velry
                 w   = ww(i,j,ivz) - velrz
                 E   = ww(i,j,irhoE) &
                     - half*rho*sqrt(velrx**2 + velry**2 + velrz**2)
                 p   = pp(i,j)

                 ! Update the residuals. Keep in mind that minus the
                 ! time update (times volume) is stored in ddw.

                 flwn  = dMadd*ss(i,j,1)
                 advtd = dPadd*p + dTadd*(E - p/(gamma(i,j)-one))
                 adhod = dTadd*(E + p)
                 fitm  = u*ss(i,j,1) + v*ss(i,j,2) + w*ss(i,j,3)

                 dww(i,j,irho)  = dww(i,j,irho)  -   flwn
                 dww(i,j,imx)   = dww(i,j,imx)   - u*flwn - advtd*ss(i,j,1)
                 dww(i,j,imy)   = dww(i,j,imy)   - v*flwn - advtd*ss(i,j,2)
                 dww(i,j,imz)   = dww(i,j,imz)   - w*flwn - advtd*ss(i,j,3)
                 dww(i,j,irhoE) = dww(i,j,irhoE) - (E+p)*flwn/rho &
                                - adhod*fitm
               enddo
             enddo

           enddo subfaceLoop

         enddo planeLoop
       enddo spectralLoop

       end subroutine level0CoolingModel
