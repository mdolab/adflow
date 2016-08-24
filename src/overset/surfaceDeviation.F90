!
!     **********************************************************************
!     *                                                                    *
!     * surfaceDeviation computes an approximation of the maximum          *       
!     * deviation a surface could be as compared to an underlying "exact"  *
!     * surface. The purpose is to compute an adaptive "near wall distance"*
!     * value that can be used to determine if a point is "close" to a     *
!     * wall. 
!     *                                                                    *
!     **********************************************************************
!

subroutine surfaceDeviation(level, sps)

  use communication 
  use blockPointers
  use overset
  use inputTimeSpectral
  use utils, only : setPointers
  implicit none

  ! Input Parameters
  integer(kind=intType), intent(in) :: level, sps

  ! Local Variables
  integer(kind=intType) :: i, j, k, ii, jj, kk, nn, iBeg, iEnd, jBeg, jEnd, mm
  real(kind=realType) :: checkDeviation, deviation
  real(kind=realType), dimension(:, :, :), pointer :: xx
  logical :: isWallType

  ! Loop over blocks
  do nn=1, nDom
     call setPointers(nn, level, sps)

     bocoLoop: do mm=1, nBocos
        wallType: if (isWallType(BCType(mm))) then
           
           ! Extract pointers for the primal wall mesh

           select case (BCFaceID(mm))
           case (iMin)
              xx => x(1, :, :, :)
           case (iMax)
              xx => x(il, :, :, :)
           case (jMin)
              xx => x(:, 1, :, :)
           case (jMax)
              xx => x(:, jl, :, :)
           case (kMin)
              xx => x(:, :, 1, :)
           case (kMax)
              xx => x(:, :, kl, :)
           end select

           jBeg = BCdata(mm)%jnBeg; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg; iEnd = BCData(mm)%inEnd  

           ! The procedure goes in 2 passes. The first pass checks all
           ! the i-direction edges, and the second all the j-direction
           ! edges. For every edge, we estimate the max deviation
           ! along that edge and then the surface will use the maximum
           ! deviation from each of the 4 edges. Ie we scatter the
           ! edge deviation to the two cells next to it. We only do
           ! the real cells here. Boundary halos get -one set (below)
           ! and then actual compute halos are set with an exchange. 

           bcData(mm)%delta = -one

           ! ------------------
           ! Check the i-edges
           ! ------------------
           do j=jBeg, jEnd       ! <------- Node loop
              do i=iBeg+1, iEnd  ! <------- Face Loop

                 ! We will creating a local cubic approximation of the
                 ! local edge. This will use node i-2, i-1, i, and
                 ! i+1. However, due to the pointer offset, these are
                 ! all shifted by 1 to get: i-1, i, i+1, i+2

                 deviation = checkDeviation(xx(i-1, j, :), xx(i, j, :), xx(i+1, j, :), &
                      xx(i+2, j, :))

                 ! Cell to the bottom:
                 if (j-1 >= jBeg+1) then 
                    bcData(mm)%delta(i, j-1) = max(bcData(mm)%delta(i, j-1), deviation)
                 end if

                 ! Cell to the top:
                 if (j+1 <= jEnd) then 
                    bcData(mm)%delta(i, j+1) = max(bcData(mm)%delta(i, j+1), deviation)
                 end if
              end do
           end do

           ! -----------------
           ! Check the j-edges
           ! -----------------
           do j=jBeg+1, jEnd   ! <------- Face loop
              do i=iBeg, iEnd  ! <------- Node Loop

                 ! We will creating a local cubic approximation of the
                 ! local edge. This will use node j-2, j-1, j, and
                 ! j+1. However, due to the pointer offset, these are
                 ! all shifted by 1 to get: j-1, j, j+1, j+2

                 deviation = checkDeviation(xx(i, j-1, :), xx(i, j, :), xx(i, j+1, :), &
                      xx(i, j+2, :))

                 ! Cell to the left:
                 if (i-1 >= iBeg+1) then 
                    bcData(mm)%delta(i-1, j) = max(bcData(mm)%delta(i-1, j), deviation)
                 end if

                 ! Cell to the right:
                 if (i+1 <= iEnd) then 
                    bcData(mm)%delta(i+1, j) = max(bcData(mm)%delta(i+1, j), deviation)
                 end if

              end do
           end do
        end if wallType
     end do bocoLoop
  end do

  ! Exchange so that halos get correct values set as well. THIS IS BROKEN FIX IT!
  !call exchangeSurfaceDelta(level, sps, commPatternCell_1st, internalCell_1st)

  ! Now make one pass back and compute a delta for the nodes. Of
  ! course, this technically makes no sense: The nodes should
  ! *exactly* match the surface by definition. However, since we are
  ! using this as a surrogate for what is near a surface, it make a
  ! little sense. Essentially we go through the nodes, and take the
  ! max deviation from the cells surrpounding it. 

  ! Loop over blocks
  do nn=1, nDom
     call setPointers(nn, level, sps)

     bocoLoop2: do mm=1, nBocos
        wallType2: if (isWallType(BCType(mm))) then
           
           jBeg = BCdata(mm)%jnBeg; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg; iEnd = BCData(mm)%inEnd  
           
           do j=jBeg, jEnd
              do i=iBeg, iEnd
                 
                 ! Since we are taking the max and the boundary halos
                 ! have a value of -one it's ok to blindy just take
                 ! the max from each of the 4 cells surrounding each node. 

                 BCData(mm)%deltaNode(i, j) = max(&
                      BCData(mm)%delta(i  , j  ), &
                      BCData(mm)%delta(i+1, j  ), &
                      BCData(mm)%delta(i  , j+1), &
                      BCData(mm)%delta(i+1, j+1))
              end do
           end do
        end if wallType2
     end do bocoLoop2
  end do

end subroutine surfaceDeviation

function checkDeviation(P0, P1, P2, P3)

  ! Find the maximum deviation between a local cubic approximation
  ! formed by nodes P0, P1, P2 and P3, with the linear approximation
  ! formed by nodes P1 and P2. 

  ! See this article for the implementation. 
  ! https://en.wikipedia.org/wiki/Centripetal_Catmull-Rom_spline

  use constants
  implicit none

  ! Input Parameters
  real(kind=realType), intent(in), dimension(3) :: P0, P1, P2, P3

  ! Function value
  real(kind=realType) ::  checkDeviation
  
  ! Working Parameters
  real(kind=realType) :: t0, t1, t2, t3
  real(kind=realType), dimension(3) :: A1, A2, A3, B1, B2
  real(kind=realType), parameter :: alpha=half
  integer(kind=intType), parameter :: N=20
  integer(kind=intType) :: i
  real(kind=realType) :: t, P(3), Q(3), s

  t0 = zero
  t1 = t0 + norm2(P1-P0)**alpha
  t2 = t1 + norm2(P2-P1)**alpha
  t3 = t2 + norm2(P3-P2)**alpha

  ! Normalize
  t1 = t1/t3
  t2 = t2/t3
  t3 = one

  ! Loop over the number of points to check. We need to go between t2
  ! and t3. No need to check the first and last since the devaition
  ! there is zero by construction.
  checkDeviation = zero

  do i=1, N
     s = (i-one)/(N-one)
     t = (one-s)*t1 + s*t2


     ! Spline pt
     A3 = (t3-t)/(t3-t2)*P2 + (t-t2)/(t3-t2)*P3
     A2 = (t2-t)/(t2-t1)*P1 + (t-t1)/(t2-t1)*P2
     A1 = (t1-t)/(t1-t0)*P0 + (t-t0)/(t1-t0)*P1
     
     B2 = (t3-t)/(t3-t1)*A2 + (t-t1)/(t3-t1)*A3
     B1 = (t2-t)/(t2-t0)*A1 + (t-t0)/(t2-t0)*A2
     
     P = (t2-t)/(t2-t1)*B1 + (t-t1)/(t2-t1)*B2

     ! Now project the cubic point onto the line to get point Q

     Q = P1 + dot_product(P-P1, P2-P1)/dot_product(P2-P1, P2-P1) * (P2 - P1)

     ! Just get the distance between the two points.
     checkDeviation = max(checkDeviation, norm2(Q-P))
     
 end do

 ! if (checkDeviation > .0055) then 
 !    print *, 'VARIABLES = X, Y, Z'
 !    print *,'ZONE I=2'
 !    print *,'DATAPACKING=POINT'
 
 !    do i=1, N
 !       s = (i-one)/(N-one)
 !       t = (one-s)*t1 + s*t2
       
 !       ! Spline pt
 !       A3 = (t3-t)/(t3-t2)*P2 + (t-t2)/(t3-t2)*P3
 !       A2 = (t2-t)/(t2-t1)*P1 + (t-t1)/(t2-t1)*P2
 !       A1 = (t1-t)/(t1-t0)*P0 + (t-t0)/(t1-t0)*P1
       
 !       B2 = (t3-t)/(t3-t1)*A2 + (t-t1)/(t3-t1)*A3
 !       B1 = (t2-t)/(t2-t0)*A1 + (t-t0)/(t2-t0)*A2
       
 !       P = (t2-t)/(t2-t1)*B1 + (t-t1)/(t2-t1)*B2
       
 !       print *, P
 !    end do
 ! end if

end function checkDeviation
