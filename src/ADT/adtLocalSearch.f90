!
!     ******************************************************************
!     *                                                                *
!     * File:          adtLocalSearch.f90                              *
!     * Author:        Edwin van der Weide                             *
!     * Starting date: 02-10-2006                                      *
!     * Last modified: 03-17-2006                                      *
!     *                                                                *
!     ******************************************************************
!
      module adtLocalSearch
!
!     ******************************************************************
!     *                                                                *
!     * Module which contains the subroutines to perform the local     *
!     * searches, i.e. the tree traversals.                            *
!     *                                                                *
!     ******************************************************************
!
      use adtUtils
      implicit none

      !=================================================================

      contains

        !===============================================================

        subroutine containmentTreeSearch(ADT,       coor,   &
                                         intInfo,   uvw,    &
                                         arrDonor,  nCoor,  &
                                         nInterpol)
!
!       ****************************************************************
!       *                                                              *
!       * This routine performs the actual containment search in the   *
!       * local tree. It is a local routine in the sense that no       *
!       * communication is involved.                                   *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * ADT:       ADT type whose ADT must be searched               *
!       * nCoor:     Number of coordinates for which the element must  *
!       *            be determined.                                    *
!       * coor:      The coordinates of these points.                  *
!       * nInterpol: Number of variables to be interpolated.           *
!       * arrDonor:  Array with the donor data; needed to obtain the   *
!       *            interpolated data.                                *
!       *                                                              *
!       * Subroutine intent(out) arguments.                            *
!       * ---------------------------------                            *
!       * intInfo: 2D integer array, in which the following output     *
!       *          be be stored:                                       *
!       *          intInfo(1,:): processor ID of the processor where   *
!       *                        the element is stored. This of course *
!       *                        is myID. If no element is found this  *
!       *                        value is set to -1.                   *
!       *          intInfo(2,:): The element type of the element.      *
!       *          intInfo(3,:): The element ID of the element in the  *
!       *                        connectivity.                         *
!       * uvw:     2D floating point array to store the parametric     *
!       *          coordinates of the point in the transformed element *
!       *          and the interpolated data:                          *
!       *          uvw(1, :): Parametric u-weight.                     *
!       *          uvw(2, :): Parametric v-weight.                     *
!       *          uvw(3, :): Parametric w-weight.                     *
!       *          uvw(4:,:): Interpolated solution, if desired. It is *
!       *                     possible to call this routine with       *
!       *                     nInterpol == 0.                          *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments.
!
        type(adtType), intent(inout) :: ADT
        integer(kind=intType), intent(in) :: nCoor
        integer(kind=intType), intent(in) :: nInterpol

        real(kind=realType), dimension(:,:), intent(in) :: coor
        real(kind=realType), dimension(:,:), intent(in) :: arrDonor

        integer(kind=intType), dimension(:,:), intent(out) :: intInfo
        real(kind=realType),   dimension(:,:), intent(out) :: uvw
!
!       Local variables.
!
        integer(kind=intType), dimension(:), pointer :: BB
        integer(kind=intType), dimension(:), pointer :: frontLeaves
        integer(kind=intType), dimension(:), pointer :: frontLeavesNew
        integer(kind=intType) :: nAllocBB, nAllocFront
        integer(kind=intType) :: ierr, nn
!
!       ****************************************************************
!       *                                                              *
!       * Begin execution.                                             *
!       *                                                              *
!       ****************************************************************
!

        ! Initial allocation of the arrays for the tree traversal.

        nAllocBB    = 10
        nAllocFront = 25

        allocate(BB(nAllocBB), frontLeaves(nAllocFront), &
                 frontLeavesNew(nAllocFront), stat=ierr)
        if(ierr /= 0)                                    &
          call adtTerminate(ADT, "containmentTreeSearch", &
                            "Memory allocation failure for BB, &
                            &frontLeaves and frontLeavesNew.")

        ! Loop over the number of coordinates to be treated.

        coorLoop: do nn=1, nCoor

           call containmentTreeSearchSinglePoint(ADT, coor(:, nn), &
           intInfo(:, nn), uvw(:, nn), arrDonor, nInterpol, BB, &
           frontLeaves, frontLeavesNew)

        end do coorLoop

        ! Release the memory allocated in this routine.

        deallocate(BB, frontLeaves, frontLeavesNew, stat=ierr)
        if(ierr /= 0)                                    &
          call adtTerminate(ADT, "containmentTreeSearch", &
                            "Deallocation failure for BB, &
                            &frontLeaves and frontLeavesNew.")

        end subroutine containmentTreeSearch

      subroutine containmentTreeSearchSinglePoint(ADT, coor, &
           intInfo, uvw, arrDonor, nInterpol, BB, &
           frontLeaves, frontLeavesNew)
!
!       ****************************************************************
!       *                                                              *
!       * This routine is replaces the original containment            *
!       * tree search, however, it has been optimized for searching a  *
!       * single query point at a time. Specifically, this means that  *
!       * there are no allocatable arrays inside this routine; they    *
!       * must be supplied exterally. Also since only a single point   *
!       * is searched we can fix some of the dimensions. This routine  *
!       * requires pointers for BB, frontLeaves and frontLeavesNew to  *
!       * be passed to the routine. Since this routine is called a     *
!       * very large number of times, these cannot be allocated        *
!       * dynamically inside for speed purposes.                       *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * ADT:       ADT type whose ADT must be searched               *
!       * coor:      The coordinate of the point to be searched.       *
!       * nInterpol: Number of variables to be interpolated.           *
!       * arrDonor:  Array with the donor data; needed to obtain the   *
!       *            interpolated data.                                *
!       *                                                              *
!       * Subroutine intent(out) arguments.                            *
!       * ---------------------------------                            *
!       * intInfo: 1D integer array of length three , in which the     *
!       *        following output                                      *
!       *          be be stored:                                       *
!       *          intInfo(1): processor ID of the processor where    
!       *                        the element is stored. This of course *
!       *                        is myID. If no element is found this  *
!       *                        value is set to -1.                   *
!       *          intInfo(2): The element type of the element.        
!       *          intInfo(3): The element ID of the element in the   
!       *                        connectivity.                         *
!       * uvw:     1D floating point array to store the parametric     *
!       *          coordinates of the point in the transformed element *
!       *          and the interpolated data:                          *
!       *          uvw(1): Parametric u-weight.                        *
!       *          uvw(2): Parametric v-weight.                        *
!       *          uvw(3): Parametric w-weight.                        *
!       *          uvw(4:): Interpolated solution(s), if desired. It is*
!       *                     possible to call this routine with       *
!       *                     nInterpol == 0.                          *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments.
!
        type(adtType), intent(inout) :: ADT
        integer(kind=intType), intent(in) :: nInterpol

        real(kind=realType), dimension(3), intent(in) :: coor
        real(kind=realType), dimension(:,:), intent(in) :: arrDonor

        integer(kind=intType), dimension(3), intent(out) :: intInfo
        real(kind=realType),   dimension(:), intent(out) :: uvw

        integer(kind=intType), dimension(:), pointer :: BB
        integer(kind=intType), dimension(:), pointer :: frontLeaves
        integer(kind=intType), dimension(:), pointer :: frontLeavesNew
!
!       Local parameters used in the Newton algorithm.
!
        integer(kind=intType), parameter :: iterMax   = 15
        real(kind=realType),   parameter :: adtEps    = 1.e-25_realType
        real(kind=realType),   parameter :: thresConv = 1.e-10_realType
!
!       Local variables.
!
        integer :: ierr

        integer(kind=intType) :: ii, kk, ll, mm, nn
        integer(kind=intType) :: nBB, nFrontLeaves, nFrontLeavesNew
        integer(kind=intType) :: nAllocBB, nAllocFront
        integer(kind=intType) :: i, nNodeElement

        integer(kind=intType), dimension(8) :: n

        real(kind=realType) :: u, v, w, uv, uw, vw, wvu, du, dv, dw
        real(kind=realType) :: oneMinusU, oneMinusV, oneMinusW
        real(kind=realType) :: oneMinusUMinusV
        real(kind=realType) :: a11, a12, a13, a21, a22, a23
        real(kind=realType) :: a31, a32, a33, val

        real(kind=realType), dimension(3)     :: x, f
        real(kind=realType), dimension(8)     :: weight
        real(kind=realType), dimension(3,2:8) :: xn

        real(kind=realType), dimension(:,:), pointer :: xBBox

        logical :: elementFound

        type(adtLeafType), dimension(:), pointer :: ADTree
!
!       ****************************************************************
!       *                                                              *
!       * Begin execution.                                             *
!       *                                                              *
!       ****************************************************************
!
        ! Set some pointers to make the code more readable.

        xBBox  => ADT%xBBox
        ADTree => ADT%ADTree

        ! Determine the sizes from the arrays we have been passed

        nAllocBB    = size(BB)
        nAllocFront = size(frontLeaves)
        
        ! Initialize the processor ID to -1 to indicate that no
        ! corresponding volume element is found.

        intInfo(1) = -1
!
!         **************************************************************
!         *                                                            *
!         * Part 1. Traverse the tree and determine the target         *
!         *         bounding boxes, which may contain the element.     *
!         *                                                            *
!         **************************************************************
!
          ! Start at the root, i.e. set the front leaf to the root leaf.
          ! Also initialize the number of possible bounding boxes to 0.

          nBB = 0

          nFrontLeaves   = 1
          frontLeaves(1) = 1

          treeTraversalLoop: do

            ! Initialize the number of leaves for the new front, i.e.
            ! the front of the next round, to 0.

            nFrontLeavesNew = 0

            ! Loop over the leaves of the current front.

            currentFrontLoop: do ii=1,nFrontLeaves

              ! Store the ID of the leaf a bit easier and loop over
              ! its two children.

              ll = frontLeaves(ii)

              childrenLoop: do mm=1,2

                ! Determine whether this child contains a bounding box
                ! or a leaf of the next level.

                kk = ADTree(ll)%children(mm)
               
                terminalTest: if(kk < 0) then

                  ! Child contains a bounding box. Check if the
                  ! coordinate is inside the bounding box.

                  kk = -kk
                  if(coor(1) >= xBBox(1,kk) .and. &
                     coor(1) <= xBBox(4,kk) .and. &
                     coor(2) >= xBBox(2,kk) .and. &
                     coor(2) <= xBBox(5,kk) .and. &
                     coor(3) >= xBBox(3,kk) .and. &
                     coor(3) <= xBBox(6,kk)) then

                    ! Coordinate is inside the bounding box. Store the
                    ! bounding box in the list of possible candidates.

                    if(nBB == nAllocBB) &
                      call reallocPlus(BB, nAllocBB, 10, ADT)

                    nBB     = nBB + 1
                    BB(nBB) = kk
                  endif

                else terminalTest

                  ! Child contains a leaf. Check if the coordinate is
                  ! inside the bounding box of the leaf.

                  if(coor(1) >= ADTree(kk)%xMin(1) .and. &
                     coor(1) <= ADTree(kk)%xMax(4) .and. &
                     coor(2) >= ADTree(kk)%xMin(2) .and. &
                     coor(2) <= ADTree(kk)%xMax(5) .and. &
                     coor(3) >= ADTree(kk)%xMin(3) .and. &
                     coor(3) <= ADTree(kk)%xMax(6)) then

                    ! Coordinate is inside the leaf. Store the leaf in
                    ! the list for the new front.

                    if(nFrontLeavesNew == nAllocFront) then
                      i = nAllocFront
                      call reallocPlus(frontLeavesNew, i, 25, ADT)
                      call reallocPlus(frontLeaves, nAllocFront, 25, ADT)
                    endif

                    nFrontLeavesNew = nFrontLeavesNew + 1
                    frontLeavesNew(nFrontLeavesNew) = kk

                  endif

                endif terminalTest

              enddo childrenLoop

            enddo currentFrontLoop

            ! End of the loop over the current front. If the new front
            ! is empty the entire tree has been traversed and an exit is
            ! made from the corresponding loop.

            if(nFrontLeavesNew == 0) exit treeTraversalLoop

            ! Copy the data of the new front leaves into the current
            ! front for the next round.

            nFrontLeaves = nFrontLeavesNew
            do ll=1,nFrontLeaves
              frontLeaves(ll) = frontLeavesNew(ll)
            enddo

          enddo treeTraversalLoop
!
!         **************************************************************
!         *                                                            *
!         * Part 2: Loop over the selected bounding boxes and check if *
!         *         the corresponding elements contain the point.      *
!         *                                                            *
!         **************************************************************
!
          elementFound = .false.

          BBoxLoop: do mm=1,nBB

            ! Determine the element type stored in this bounding box.

            kk = BB(mm)
            select case (ADT%elementType(kk))

              case (adtTetrahedron)

                ! Element is a tetrahedron.
                ! Compute the coordinates relative to node 1.

                ll   = ADT%elementID(kk)
                n(1) = ADT%tetraConn(1,ll)

                do i=2,4
                  n(i) = ADT%tetraConn(i,ll)

                  xn(1,i) = ADT%coor(1,n(i)) - ADT%coor(1,n(1))
                  xn(2,i) = ADT%coor(2,n(i)) - ADT%coor(2,n(1))
                  xn(3,i) = ADT%coor(3,n(i)) - ADT%coor(3,n(1))
                enddo

                x(1) = coor(1) - ADT%coor(1,n(1))
                x(2) = coor(2) - ADT%coor(2,n(1))
                x(3) = coor(3) - ADT%coor(3,n(1))

                ! Determine the matrix for the linear transformation
                ! from the standard element to the current element.

                a11 = xn(1,2); a12 = xn(1,3); a13 = xn(1,4)
                a21 = xn(2,2); a22 = xn(2,3); a23 = xn(2,4)
                a31 = xn(3,2); a32 = xn(3,3); a33 = xn(3,4)

                ! Compute the determinant. Make sure that it is not zero
                ! and invert the value.

                val = a11*(a22*a33 - a32*a23) + a21*(a13*a32 - a12*a33) &
                    + a31*(a12*a23 - a13*a22)
                val = sign(adtOne,val)/max(abs(val),adtEps)

                ! Compute the u, v, w weights for the given coordinate.

                u = val*((a22*a33 - a23*a32)*x(1) &
                  +      (a13*a32 - a12*a33)*x(2) &
                  +      (a12*a23 - a13*a22)*x(3))
                v = val*((a23*a31 - a21*a33)*x(1) &
                  +      (a11*a33 - a13*a31)*x(2) &
                  +      (a13*a21 - a11*a23)*x(3))
                w = val*((a21*a32 - a22*a31)*x(1) &
                  +      (a12*a31 - a11*a32)*x(2) &
                  +      (a11*a22 - a12*a21)*x(3))

                ! Check if the coordinate is inside the tetrahedron.
                ! If so, set elementFound to .true. and determine the
                ! interpolation weights.

                if(u >= adtZero .and. v >= adtZero .and. &
                   w >= adtZero .and. (u+v+w) <= adtOne) then
                  elementFound = .true.

                  ! Set the number of interpolation nodes to 4 and
                  ! determine the interpolation weights.

                  nNodeElement = 4

                  weight(1) = adtOne - u - v - w
                  weight(2) = u
                  weight(3) = v
                  weight(4) = w
                endif

              !=========================================================

              case (adtPyramid)

                ! Element is a pyramid.
                ! Compute the coordinates relative to node 1.

                ll   = ADT%elementID(kk)
                n(1) = ADT%pyraConn(1,ll)

                do i=2,5
                  n(i) = ADT%pyraConn(i,ll)

                  xn(1,i) = ADT%coor(1,n(i)) - ADT%coor(1,n(1))
                  xn(2,i) = ADT%coor(2,n(i)) - ADT%coor(2,n(1))
                  xn(3,i) = ADT%coor(3,n(i)) - ADT%coor(3,n(1))
                enddo

                x(1) = coor(1) - ADT%coor(1,n(1))
                x(2) = coor(2) - ADT%coor(2,n(1))
                x(3) = coor(3) - ADT%coor(3,n(1))

                ! Modify the coordinates of node 3, such that it
                ! corresponds to the weights of the u*v term in the
                ! transformation.

                xn(1,3) = xn(1,3) - xn(1,2) - xn(1,4)
                xn(2,3) = xn(2,3) - xn(2,2) - xn(2,4)
                xn(3,3) = xn(3,3) - xn(3,2) - xn(3,4)

                ! Set the starting values of u, v and w such that it is
                ! somewhere in the middle of the element. In this way the
                ! Jacobian matrix is always regular, even if the element
                ! is degenerate.

                u = adtHalf; v = adtHalf; w = adtHalf

                ! The Newton algorithm to determine the parametric
                ! weights u, v and w for the given coordinate.

                NewtonPyra: do ll=1,iterMax

                  ! Compute the RHS.

                  uv = u*v
                  oneMinusW = adtOne - w

                  f(1) = oneMinusW*(xn(1,2)*u + xn(1,4)*v + xn(1,3)*uv) &
                       + xn(1,5)*w - x(1)
                  f(2) = oneMinusW*(xn(2,2)*u + xn(2,4)*v + xn(2,3)*uv) &
                       + xn(2,5)*w - x(2)
                  f(3) = oneMinusW*(xn(3,2)*u + xn(3,4)*v + xn(3,3)*uv) &
                       + xn(3,5)*w - x(3)

                  ! Compute the Jacobian.

                  a11 = oneMinusW*(xn(1,2) + xn(1,3)*v)
                  a12 = oneMinusW*(xn(1,4) + xn(1,3)*u)
                  a13 = xn(1,5) - xn(1,2)*u - xn(1,4)*v - xn(1,3)*uv

                  a21 = oneMinusW*(xn(2,2) + xn(2,3)*v)
                  a22 = oneMinusW*(xn(2,4) + xn(2,3)*u)
                  a23 = xn(2,5) - xn(2,2)*u - xn(2,4)*v - xn(2,3)*uv

                  a31 = oneMinusW*(xn(3,2) + xn(3,3)*v)
                  a32 = oneMinusW*(xn(3,4) + xn(3,3)*u)
                  a33 = xn(3,5) - xn(3,2)*u - xn(3,4)*v - xn(3,3)*uv

                  ! Compute the determinant. Make sure that it is not zero
                  ! and invert the value. The cut off is needed to be able
                  ! to handle exceptional cases for degenerate elements.

                  val = a11*(a22*a33 - a32*a23) + a21*(a13*a32 - a12*a33) &
                      + a31*(a12*a23 - a13*a22)
                  val = sign(adtOne,val)/max(abs(val),adtEps)

                  ! Compute the new values of u, v and w.

                  du = val*((a22*a33 - a23*a32)*f(1) &
                     +      (a13*a32 - a12*a33)*f(2) &
                     +      (a12*a23 - a13*a22)*f(3))
                  dv = val*((a23*a31 - a21*a33)*f(1) &
                     +      (a11*a33 - a13*a31)*f(2) &
                     +      (a13*a21 - a11*a23)*f(3))
                  dw = val*((a21*a32 - a22*a31)*f(1) &
                     +      (a12*a31 - a11*a32)*f(2) &
                     +      (a11*a22 - a12*a21)*f(3))

                  u = u - du; v = v - dv; w = w - dw

                  ! Exit the loop if the update of the parametric
                  ! weights is below the threshold

                  val = sqrt(du*du + dv*dv + dw*dw)
                  if(val <= thresConv) exit NewtonPyra
                enddo NewtonPyra

                ! Check if the coordinate is inside the pyramid.
                ! If so, set elementFound to .true. and determine the
                ! interpolation weights.

                if(u     >= adtZero .and. v     >= adtZero .and. &
                   w     >= adtZero .and. (u+w) <= adtOne  .and. &
                   (v+w) <= adtOne) then
                  elementFound = .true.

                  ! Set the number of interpolation nodes to 5 and
                  ! determine the interpolation weights.

                  nNodeElement = 5

                  oneMinusU = adtOne - u
                  oneMinusV = adtOne - v
                  oneMinusW = adtOne - w

                  weight(1) = oneMinusU*oneMinusV*oneMinusW
                  weight(2) =         u*oneMinusV*oneMinusW
                  weight(3) =         u*        v*oneMinusW
                  weight(4) = oneMinusU*        v*oneMinusW
                  weight(5) =                             w
                endif

              !=========================================================

              case (adtPrism)

                ! Element is a prism.
                ! Compute the coordinates relative to node 1.

                 ll   = ADT%elementID(kk)
                 n(1) = ADT%prismsConn(1,ll)

                do i=2,6
                  n(i) = ADT%prismsConn(i,ll)

                  xn(1,i) = ADT%coor(1,n(i)) - ADT%coor(1,n(1))
                  xn(2,i) = ADT%coor(2,n(i)) - ADT%coor(2,n(1))
                  xn(3,i) = ADT%coor(3,n(i)) - ADT%coor(3,n(1))
                enddo

                x(1) = coor(1) - ADT%coor(1,n(1))
                x(2) = coor(2) - ADT%coor(2,n(1))
                x(3) = coor(3) - ADT%coor(3,n(1))

                ! Modify the coordinates of node 5 and 6, such that they
                ! correspond to the weights of the u*w and v*w term in the
                ! transformation respectively.

                xn(1,5) = xn(1,5) - xn(1,2) - xn(1,4)
                xn(2,5) = xn(2,5) - xn(2,2) - xn(2,4)
                xn(3,5) = xn(3,5) - xn(3,2) - xn(3,4)

                xn(1,6) = xn(1,6) - xn(1,3) - xn(1,4)
                xn(2,6) = xn(2,6) - xn(2,3) - xn(2,4)
                xn(3,6) = xn(3,6) - xn(3,3) - xn(3,4)

                ! Set the starting values of u, v and w such that it is
                ! somewhere in the middle of the element. In this way the
                ! Jacobian matrix is always regular, even if the element
                ! is degenerate.

                u = adtFourth; v = adtFourth; w = adtHalf

                ! The Newton algorithm to determine the parametric
                ! weights u, v and w for the given coordinate.

                NewtonPrisms: do ll=1,iterMax

                  ! Compute the RHS.

                  uw = u*w; vw = v*w

                  f(1) = xn(1,2)*u  + xn(1,3)*v  + xn(1,4)*w &
                       + xn(1,5)*uw + xn(1,6)*vw - x(1)
                  f(2) = xn(2,2)*u  + xn(2,3)*v  + xn(2,4)*w &
                       + xn(2,5)*uw + xn(2,6)*vw - x(2)
                  f(3) = xn(3,2)*u  + xn(3,3)*v  + xn(3,4)*w &
                       + xn(3,5)*uw + xn(3,6)*vw - x(3)

                  ! Compute the Jacobian.

                  a11 = xn(1,2) + xn(1,5)*w
                  a12 = xn(1,3) + xn(1,6)*w
                  a13 = xn(1,4) + xn(1,5)*u + xn(1,6)*v

                  a21 = xn(2,2) + xn(2,5)*w
                  a22 = xn(2,3) + xn(2,6)*w
                  a23 = xn(2,4) + xn(2,5)*u + xn(2,6)*v

                  a31 = xn(3,2) + xn(3,5)*w
                  a32 = xn(3,3) + xn(3,6)*w
                  a33 = xn(3,4) + xn(3,5)*u + xn(3,6)*v

                  ! Compute the determinant. Make sure that it is not zero
                  ! and invert the value. The cut off is needed to be able
                  ! to handle exceptional cases for degenerate elements.

                  val = a11*(a22*a33 - a32*a23) + a21*(a13*a32 - a12*a33) &
                      + a31*(a12*a23 - a13*a22)
                  val = sign(adtOne,val)/max(abs(val),adtEps)

                  ! Compute the new values of u, v and w.

                  du = val*((a22*a33 - a23*a32)*f(1) &
                     +      (a13*a32 - a12*a33)*f(2) &
                     +      (a12*a23 - a13*a22)*f(3))
                  dv = val*((a23*a31 - a21*a33)*f(1) &
                     +      (a11*a33 - a13*a31)*f(2) &
                     +      (a13*a21 - a11*a23)*f(3))
                  dw = val*((a21*a32 - a22*a31)*f(1) &
                     +      (a12*a31 - a11*a32)*f(2) &
                     +      (a11*a22 - a12*a21)*f(3))

                  u = u - du; v = v - dv; w = w - dw

                  ! Exit the loop if the update of the parametric
                  ! weights is below the threshold

                  val = sqrt(du*du + dv*dv + dw*dw)
                  if(val <= thresConv) exit NewtonPrisms
                enddo NewtonPrisms

                ! Check if the coordinate is inside the prism.
                ! If so, set elementFound to .true. and determine the
                ! interpolation weights.

                if(u     >= adtZero .and. v >= adtZero .and. &
                   w     >= adtZero .and. w <= adtOne  .and. &
                   (u+v) <= adtOne) then
                  elementFound = .true.

                  ! Set the number of interpolation nodes to 6 and
                  ! determine the interpolation weights.

                  nNodeElement = 6

                  oneMinusUminusV = adtOne - u - v
                  oneMinusW       = adtOne - w

                  weight(1) = oneMinusUminusV*oneMinusW
                  weight(2) =               u*oneMinusW
                  weight(3) =               v*oneMinusW
                  weight(4) = oneMinusUminusV*        w
                  weight(5) =               u*        w
                  weight(6) =               v*        w
                endif

              !=========================================================

              case (adtHexahedron)

                ! Element is a hexahedron.
                ! Compute the coordinates relative to node 1.

                ll   = ADT%elementID(kk)
                n(1) = ADT%hexaConn(1,ll)

                do i=2,8
                  n(i) = ADT%hexaConn(i,ll)

                  xn(1,i) = ADT%coor(1,n(i)) - ADT%coor(1,n(1))
                  xn(2,i) = ADT%coor(2,n(i)) - ADT%coor(2,n(1))
                  xn(3,i) = ADT%coor(3,n(i)) - ADT%coor(3,n(1))
                enddo

                x(1) = coor(1) - ADT%coor(1,n(1))
                x(2) = coor(2) - ADT%coor(2,n(1))
                x(3) = coor(3) - ADT%coor(3,n(1))

                ! Modify the coordinates of node 3, 6, 8 and 7 such that
                ! they correspond to the weights of the u*v, u*w, v*w and
                ! u*v*w term in the transformation respectively.

                xn(1,7) = xn(1,7) + xn(1,2) + xn(1,4) + xn(1,5) &
                        -           xn(1,3) - xn(1,6) - xn(1,8)
                xn(2,7) = xn(2,7) + xn(2,2) + xn(2,4) + xn(2,5) &
                        -           xn(2,3) - xn(2,6) - xn(2,8)
                xn(3,7) = xn(3,7) + xn(3,2) + xn(3,4) + xn(3,5) &
                        -           xn(3,3) - xn(3,6) - xn(3,8)

                xn(1,3) = xn(1,3) - xn(1,2) - xn(1,4)
                xn(2,3) = xn(2,3) - xn(2,2) - xn(2,4)
                xn(3,3) = xn(3,3) - xn(3,2) - xn(3,4)

                xn(1,6) = xn(1,6) - xn(1,2) - xn(1,5)
                xn(2,6) = xn(2,6) - xn(2,2) - xn(2,5)
                xn(3,6) = xn(3,6) - xn(3,2) - xn(3,5)

                xn(1,8) = xn(1,8) - xn(1,4) - xn(1,5)
                xn(2,8) = xn(2,8) - xn(2,4) - xn(2,5)
                xn(3,8) = xn(3,8) - xn(3,4) - xn(3,5)

                ! Set the starting values of u, v and w such that it is
                ! somewhere in the middle of the element. In this way the
                ! Jacobian matrix is always regular, even if the element
                ! is degenerate.

                u = adtHalf; v = adtHalf; w = adtHalf

                ! The Newton algorithm to determine the parametric
                ! weights u, v and w for the given coordinate.

                NewtonHexa: do ll=1,iterMax

                  ! Compute the RHS.

                  uv = u*v; uw = u*w; vw = v*w; wvu = u*v*w

                  f(1) = xn(1,2)*u   + xn(1,4)*v  + xn(1,5)*w  &
                       + xn(1,3)*uv  + xn(1,6)*uw + xn(1,8)*vw &
                       + xn(1,7)*wvu - x(1)
                  f(2) = xn(2,2)*u   + xn(2,4)*v  + xn(2,5)*w  &
                       + xn(2,3)*uv  + xn(2,6)*uw + xn(2,8)*vw &
                       + xn(2,7)*wvu - x(2)
                  f(3) = xn(3,2)*u   + xn(3,4)*v  + xn(3,5)*w  &
                       + xn(3,3)*uv  + xn(3,6)*uw + xn(3,8)*vw &
                       + xn(3,7)*wvu - x(3)

                  ! Compute the Jacobian.

                  a11 = xn(1,2) + xn(1,3)*v + xn(1,6)*w + xn(1,7)*vw
                  a12 = xn(1,4) + xn(1,3)*u + xn(1,8)*w + xn(1,7)*uw
                  a13 = xn(1,5) + xn(1,6)*u + xn(1,8)*v + xn(1,7)*uv

                  a21 = xn(2,2) + xn(2,3)*v + xn(2,6)*w + xn(2,7)*vw
                  a22 = xn(2,4) + xn(2,3)*u + xn(2,8)*w + xn(2,7)*uw
                  a23 = xn(2,5) + xn(2,6)*u + xn(2,8)*v + xn(2,7)*uv

                  a31 = xn(3,2) + xn(3,3)*v + xn(3,6)*w + xn(3,7)*vw
                  a32 = xn(3,4) + xn(3,3)*u + xn(3,8)*w + xn(3,7)*uw
                  a33 = xn(3,5) + xn(3,6)*u + xn(3,8)*v + xn(3,7)*uv

                  ! Compute the determinant. Make sure that it is not zero
                  ! and invert the value. The cut off is needed to be able
                  ! to handle exceptional cases for degenerate elements.

                  val = a11*(a22*a33 - a32*a23) + a21*(a13*a32 - a12*a33) &
                      + a31*(a12*a23 - a13*a22)
                  val = sign(adtOne,val)/max(abs(val),adtEps)

                  ! Compute the new values of u, v and w.

                  du = val*((a22*a33 - a23*a32)*f(1) &
                     +      (a13*a32 - a12*a33)*f(2) &
                     +      (a12*a23 - a13*a22)*f(3))
                  dv = val*((a23*a31 - a21*a33)*f(1) &
                     +      (a11*a33 - a13*a31)*f(2) &
                     +      (a13*a21 - a11*a23)*f(3))
                  dw = val*((a21*a32 - a22*a31)*f(1) &
                     +      (a12*a31 - a11*a32)*f(2) &
                     +      (a11*a22 - a12*a21)*f(3))

                  u = u - du; v = v - dv; w = w - dw

                  ! Exit the loop if the update of the parametric
                  ! weights is below the threshold

                  val = sqrt(du*du + dv*dv + dw*dw)
                  if(val <= thresConv) exit NewtonHexa
                enddo NewtonHexa

                ! Check if the coordinate is inside the hexahedron.
                ! If so, set elementFound to .true. and determine the
                ! interpolation weights.

                if(u >= adtZero .and. u <= adtOne .and. &
                   v >= adtZero .and. v <= adtOne .and. &
                   w >= adtZero .and. w <= adtOne) then
                  elementFound = .true.

                  ! Set the number of interpolation nodes to 8 and
                  ! determine the interpolation weights.

                  nNodeElement = 8

                  oneMinusU = adtOne - u
                  oneMinusV = adtOne - v
                  oneMinusW = adtOne - w

                  weight(1) = oneMinusU*oneMinusV*oneMinusW
                  weight(2) =         u*oneMinusV*oneMinusW
                  weight(3) =         u*        v*oneMinusW
                  weight(4) = oneMinusU*        v*oneMinusW
                  weight(5) = oneMinusU*oneMinusV*        w
                  weight(6) =         u*oneMinusV*        w
                  weight(7) =         u*        v*        w
                  weight(8) = oneMinusU*        v*        w
                endif

            end select

            ! If the coordinate is inside the element store all the
            ! necessary information and exit the loop over the target
            ! bounding boxes.

            if( elementFound ) then

              ! The processor, element type and local element ID.

              intInfo(1) = ADT%myID
              intInfo(2) = ADT%elementType(kk)
              intInfo(3) = ADT%elementID(kk)

              ! The parametric weights.

              uvw(1) = u
              uvw(2) = v
              uvw(3) = w

              ! The interpolated solution.

              do ll=1,nInterpol
                ii = 3+ll
                uvw(ii) = weight(1)*arrDonor(ll,n(1))
                do i=2,nNodeElement
                  uvw(ii) = uvw(ii) + weight(i)*arrDonor(ll,n(i))
                enddo
              enddo

              ! And exit the loop over the bounding boxes.

              exit BBoxLoop
            endif

          enddo BBoxLoop

        end subroutine containmentTreeSearchSinglePoint


        !***************************************************************
        !***************************************************************

        subroutine minDistanceTreeSearch(ADT,       coor,   &
                                         intInfo,   uvw,    &
                                         arrDonor,  nCoor,  &
                                         nInterpol)
!
!       ****************************************************************
!       *                                                              *
!       * This routine performs the actual minimum distance search in  *
!       * the local tree. It is a local routine in the sense that no   *
!       * communication is involved.                                   *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * ADT:       ADT type whose ADT must be searched               *
!       * nCoor:     Number of coordinates for which the element must  *
!       *            be determined.                                    *
!       * coor:      The coordinates and the currently stored minimum  *
!       *            distance squared of these points:                 *
!       *            coor(1,;): Coordinate 1.                          *
!       *            coor(2,;): Coordinate 2.                          *
!       *            coor(3,;): Coordinate 3.                          *
!       *            coor(4,;): The currently stored minimum distance  *
!       *            squared.                                          *
!       * nInterpol: Number of variables to be interpolated.           *
!       * arrDonor:  Array with the donor data; needed to obtain the   *
!       *            interpolated data.                                *
!       *                                                              *
!       * Subroutine intent(out) arguments.                            *
!       * ---------------------------------                            *
!       * intInfo: 2D integer array, in which the following output     *
!       *          will be stored:                                     *
!       *          intInfo(1,:): processor ID of the processor where   *
!       *                        the element is stored. This of course *
!       *                        is myID. If no element is found this  *
!       *                        value is set to -1.                   *
!       *          intInfo(2,:): The element type of the element.      *
!       *          intInfo(3,:): The element ID of the element in the  *
!       *                        the connectivity.                     *
!       * uvw:     2D floating point array to store the parametric     *
!       *          coordinates of the point in the transformed element *
!       *          as well as the new distance squared and the         *
!       *          interpolated solution:                              *
!       *          uvw(1, :): Parametric u-weight.                     *
!       *          uvw(2, :): Parametric v-weight.                     *
!       *          uvw(3, :): Parametric w-weight.                     *
!       *          uvw(4, :): The new distance squared.                *
!       *          uvw(5:,:): Interpolated solution, if desired. It is *
!       *                     possible to call this routine with       *
!       *                     nInterpol == 0.                          *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments.
!
        type(adtType), intent(inout) :: ADT
        integer(kind=intType), intent(in) :: nCoor
        integer(kind=intType), intent(in) :: nInterpol

        real(kind=realType), dimension(:,:), intent(in) :: coor
        real(kind=realType), dimension(:,:), intent(in) :: arrDonor

        integer(kind=intType), dimension(:,:), intent(out) :: intInfo
        real(kind=realType),   dimension(:,:), intent(out) :: uvw
!
!       Local parameters used in the Newton algorithm.
!
        integer(kind=intType), parameter :: iterMax   = 15
        real(kind=realType),   parameter :: adtEps    = 1.e-25_realType
        real(kind=realType),   parameter :: thresConv = 1.e-10_realType
!
!       Local variables.
!
        integer :: ierr, nn
        integer(kind=intType) :: nAllocBB, nAllocFront, nStack
        integer(kind=intType), dimension(:), pointer :: frontLeaves
        integer(kind=intType), dimension(:), pointer :: frontLeavesNew
        type(adtBBoxTargetType), dimension(:), pointer :: BB
!
!       ****************************************************************
!       *                                                              *
!       * Begin execution.                                             *
!       *                                                              *
!       ****************************************************************
!
   
        ! Initial allocation of the arrays for the tree traversal as well
        ! as the stack array used in the qsort routine. The latter is
        ! done, because the qsort routine is called for every coordinate
        ! and therefore it is more efficient to allocate the stack once
        ! rather than over and over again. The disadvantage of course is
        ! that an essentially local variable, stack, is now stored in
        ! adtData.

        nAllocBB    =  10
        nAllocFront =  25
        nStack      = 100

        allocate(stack(nStack), BB(nAllocBB), frontLeaves(nAllocFront), &
                 frontLeavesNew(nAllocFront), stat=ierr)
        if(ierr /= 0)                                    &
          call adtTerminate(ADT, "minDistanceTreeSearch", &
                            "Memory allocation failure for stack, BB, &
                            &etc.")

        ! Loop over the number of coordinates to be treated.

        coorLoop: do nn=1,nCoor

           call minDistanceTreeSearchSinglePoint(ADT, coor(:, nn), &
                intInfo(:, nn), uvw(:, nn), arrDonor, nInterpol, BB, &
                frontLeaves, frontLeavesNew)

        enddo coorLoop

        ! Release the memory allocated in this routine.

        deallocate(stack, BB, frontLeaves, frontLeavesNew, stat=ierr)
        if(ierr /= 0)                                    &
          call adtTerminate(ADT, "minDistanceTreeSearch", &
                            "Deallocation failure for stack, BB, etc.")

        end subroutine minDistanceTreeSearch


      subroutine minDistanceTreeSearchSinglePoint(ADT, coor, intInfo, &
           uvw, arrDonor, nInterpol, BB, frontLeaves, frontLeavesNew)
!
!       ****************************************************************
!       *                                                              *
!       * This routine performs the actual minimum distance search for *
!       * a single point on the local tree. It is local in the sens    *
!       * that no communication is involved. This routine does the     *
!       * actual search. The minDistanceTreeSearch is just a wrapper   *
!       * around this routine. The reason for the split is that the    *
!       * overset mesh connectivity requires efficient calling with    *
!       * a single coordinate. Therefore, this rouine does not         *
!       * allocate/deallocate any variables.                           *
!       *                                                              *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * ADT:       ADT type whose ADT must be searched               *
!       * coor:      The coordinates and the currently stored minimum  *
!       *            distance squared of these points:                 *
!       *            coor(1): Coordinate 1.                            *
!       *            coor(2): Coordinate 2.                            *
!       *            coor(3): Coordinate 3.                            *
!       *            coor(4): The currently stored minimum distance    *
!       *            squared.                                          *
!       * nInterpol: Number of variables to be interpolated.           *
!       * arrDonor:  Array with the donor data; needed to obtain the   *
!       *            interpolated data.                                *
!       *                                                              *
!       * Subroutine intent(out) arguments.                            *
!       * ---------------------------------                            *
!       * intInfo: 1D integer array, in which the following output     *
!       *          will be stored:                                     *
!       *          intInfo(1): processor ID of the processor where   *
!       *                        the element is stored. This of course *
!       *                        is myID. If no element is found this  *
!       *                        value is set to -1.                   *
!       *          intInfo(2): The element type of the element.        *
!       *          intInfo(3): The element ID of the element in the    *
!       *                        the connectivity.                     *
!       * uvw:     2D floating point array to store the parametric     *
!       *          coordinates of the point in the transformed element *
!       *          as well as the new distance squared and the         *
!       *          interpolated solution:                              *
!       *          uvw(1): Parametric u-weight.                        *
!       *          uvw(2): Parametric v-weight.                        *
!       *          uvw(3): Parametric w-weight.                        *
!       *          uvw(4): The new distance squared.                   *
!       *          uvw(5): Interpolated solution, if desired. It is    *
!       *                     possible to call this routine with       *
!       *                     nInterpol == 0.                          *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments.
!
        type(adtType), intent(inout) :: ADT
        integer(kind=intType), intent(in) :: nInterpol

        real(kind=realType), dimension(4), intent(in) :: coor
        real(kind=realType), dimension(:,:), intent(in) :: arrDonor

        integer(kind=intType), dimension(3), intent(out) :: intInfo
        real(kind=realType),   dimension(5), intent(out) :: uvw
        integer(kind=intType), dimension(:), pointer :: frontLeaves
        integer(kind=intType), dimension(:), pointer :: frontLeavesNew
        type(adtBBoxTargetType), dimension(:), pointer :: BB
!
!       Local parameters used in the Newton algorithm.
!
        integer(kind=intType), parameter :: iterMax   = 15
        real(kind=realType),   parameter :: adtEps    = 1.e-25_realType
        real(kind=realType),   parameter :: thresConv = 1.e-10_realType
!
!       Local variables.
!
        integer :: ierr

        integer(kind=intType) :: ii, kk, ll, mm, nn, activeLeaf
        integer(kind=intType) :: nBB, nFrontLeaves, nFrontLeavesNew
        integer(kind=intType) :: nAllocBB, nAllocFront, nNodeElement
        integer(kind=intType) :: i, kkk

        integer(kind=intType), dimension(8) :: n, m

        real(kind=realType) :: dx, dy, dz, d1, d2, invLen, val
        real(kind=realType) :: u, v, w, uv, uold, vold, vn, du, dv
        real(kind=realType) :: uu, vv, ww

        real(kind=realType), dimension(2) :: dd
        real(kind=realType), dimension(3) :: x1, x21, x41, x3142, xf
        real(kind=realType), dimension(3) :: vf, vt, a, b, norm, an, bn
        real(kind=realType), dimension(3) :: chi
        real(kind=realType), dimension(8) :: weight

        real(kind=realType), dimension(:,:), pointer :: xBBox

        logical :: elementFound
        type(adtLeafType),       dimension(:), pointer :: ADTree
!
!       ****************************************************************
!       *                                                              *
!       * Begin execution.                                             *
!       *                                                              *
!       ****************************************************************
!
        ! Set some pointers to make the code more readable.

        xBBox  => ADT%xBBox
        ADTree => ADT%ADTree

        ! Initial allocation of the arrays for the tree traversal as well
        ! as the stack array used in the qsort routine. The latter is
        ! done, because the qsort routine is called for every coordinate
        ! and therefore it is more efficient to allocate the stack once
        ! rather than over and over again. The disadvantage of course is
        ! that an essentially local variable, stack, is now stored in
        ! adtData.

        nAllocBB    =  size(BB)
        nAllocFront =  size(frontLeaves)
        nStack      =  size(stack)

        ! Initialize the processor ID to -1 to indicate that no
        ! corresponding volume element is found and the new minimum
        ! distance squared to the old value.
        
        intInfo(1)     = -1
        uvw(4) = coor(4)
        !
!         **************************************************************
!         *                                                            *
!         * Part 1. Determine the possible minimum distance squared to *
!         *         the root leaf. If larger than the current distance *
!         *         there is no need to search this tree.              *
!         *                                                            *
!         **************************************************************
!
          if(     coor(1) < ADTree(1)%xMin(1)) then
            dx =  coor(1) - ADTree(1)%xMin(1)
          else if(coor(1) > ADTree(1)%xMax(4)) then
            dx =  coor(1) - ADTree(1)%xMax(4)
          else
            dx = adtZero
          endif

          if(     coor(2) < ADTree(1)%xMin(2)) then
            dy =  coor(2) - ADTree(1)%xMin(2)
          else if(coor(2) > ADTree(1)%xMax(5)) then
            dy =  coor(2) - ADTree(1)%xMax(5)
          else
            dy = adtZero
          endif

          if(     coor(3) < ADTree(1)%xMin(3)) then
            dz =  coor(3) - ADTree(1)%xMin(3)
          else if(coor(3) > ADTree(1)%xMax(6)) then
            dz =  coor(3) - ADTree(1)%xMax(6)
          else
            dz = adtZero
          endif

          ! Continue with the next coordinate if the possible distance
          ! squared to the root leaf is larger than the currently stored
          ! value.

          if((dx*dx + dy*dy + dz*dz) >= uvw(4)) return
!
!         **************************************************************
!         *                                                            *
!         * Part 2. Find a likely bounding box, which minimizes the    *
!         *         guaranteed distance.                               *
!         *                                                            *
!         **************************************************************
!
          activeLeaf = 1

          ! Traverse the tree until a terminal leaf is found.

          treeTraversal1: do

            ! Exit the loop when a terminal leaf has been found.
            ! This is indicated by a negative value of activeLeaf.

            if(activeLeaf < 0) exit treeTraversal1

            ! Determine the guaranteed distance squared for both children.

            do mm=1,2

              ! Determine whether the child contains a bounding box or
              ! another leaf of the tree.

              ll = ADTree(activeLeaf)%children(mm)
              if(ll > 0) then

                ! Child contains a leaf. Determine the guaranteed distance
                ! vector to the leaf.

                d1 = abs(coor(1) - ADTree(ll)%xMin(1))
                d2 = abs(coor(1) - ADTree(ll)%xMax(4))
                dx = max(d1,d2)

                d1 = abs(coor(2) - ADTree(ll)%xMin(2))
                d2 = abs(coor(2) - ADTree(ll)%xMax(5))
                dy = max(d1,d2)

                d1 = abs(coor(3) - ADTree(ll)%xMin(3))
                d2 = abs(coor(3) - ADTree(ll)%xMax(6))
                dz = max(d1,d2)

              else

                ! Child contains a bounding box. Determine the guaranteed
                ! distance vector to it.

                ll = -ll

                d1 = abs(coor(1) - xBBox(1,ll))
                d2 = abs(coor(1) - xBBox(4,ll))
                dx = max(d1,d2)

                d1 = abs(coor(2) - xBBox(2,ll))
                d2 = abs(coor(2) - xBBox(5,ll))
                dy = max(d1,d2)

                d1 = abs(coor(3) - xBBox(3,ll))
                d2 = abs(coor(3) - xBBox(6,ll))
                dz = max(d1,d2)

              endif

              ! Compute the guaranteed distance squared for this child.

              dd(mm) = dx*dx + dy*dy + dz*dz

            enddo

            ! Determine which will be the next leaf in the tree traversal.
            ! This will be the leaf which has the minimum guaranteed
            ! distance. In case of ties take the left leaf, because this
            ! leaf may have more children.

            if(dd(1) <= dd(2)) then
              activeLeaf = ADTree(activeLeaf)%children(1)
            else
              activeLeaf = ADTree(activeLeaf)%children(2)
            endif

          enddo treeTraversal1

          ! Store the minimum of the just computed guaranteed distance
          ! squared and the currently stored value in uvw.

          uvw(4) = min(uvw(4),dd(1),dd(2))
!
!         **************************************************************
!         *                                                            *
!         * Part 3. Find the bounding boxes whose possible minimum     *
!         *         distance is less than the currently stored value.  *
!         *                                                            *
!         **************************************************************
!
          ! In part 1 it was already tested that the possible distance
          ! squared of the root leaf was less than the current value.
          ! Therefore initialize the current front to the root leaf and
          ! set the number of bounding boxes to 0.

          nBB = 0

          nFrontLeaves   = 1
          frontLeaves(1) = 1

          ! Second tree traversal. Now to find all possible bounding
          ! box candidates.

          treeTraversal2: do

            ! Initialize the number of leaves for the new front, i.e.
            ! the front of the next round, to 0.

            nFrontLeavesNew = 0

            ! Loop over the leaves of the current front.

            currentFrontLoop: do ii=1,nFrontLeaves

              ! Store the ID of the leaf a bit easier and loop over
              ! its two children.

              ll = frontLeaves(ii)

              childrenLoop: do mm=1,2

                ! Determine whether this child contains a bounding box
                ! or a leaf of the next level.

                kk = ADTree(ll)%children(mm)
                terminalTest: if(kk < 0) then

                  ! Child contains a bounding box. Determine the possible
                  ! minimum distance squared to this bounding box.

                  kk = -kk

                  if(     coor(1) < xBBox(1,kk)) then
                    dx =  coor(1) - xBBox(1,kk)
                  else if(coor(1) > xBBox(4,kk)) then
                    dx =  coor(1) - xBBox(4,kk)
                  else
                    dx = adtZero
                  endif

                  if(     coor(2) < xBBox(2,kk)) then
                    dy =  coor(2) - xBBox(2,kk)
                  else if(coor(2) > xBBox(5,kk)) then
                    dy =  coor(2) - xBBox(5,kk)
                  else
                    dy = adtZero
                  endif

                  if(     coor(3) < xBBox(3,kk)) then
                    dz =  coor(3) - xBBox(3,kk)
                  else if(coor(3) > xBBox(6,kk)) then
                    dz =  coor(3) - xBBox(6,kk)
                  else
                    dz = adtZero
                  endif

                  d2 = dx*dx + dy*dy + dz*dz

                  ! If this distance squared is less than the current
                  ! value, store this bounding box as a target.

                  testStoreBBox: if(d2 < uvw(4)) then

                    ! Check if the memory must be reallocated.

                    if(nBB == nAllocBB) &
                      call reallocBBoxTargetTypePlus(BB, nAllocBB, &
                                                     10, ADT)

                    ! Update the counter and store the data.

                    nBB              = nBB + 1
                    BB(nBB)%ID       = kk
                    BB(nBB)%posDist2 = d2

                    ! Although in step 2, i.e. the first tree traversal,
                    ! the guaranteed distance squared to a bounding box
                    ! has already been computed, this has been done only
                    ! for a likely candidate and not for all the possible
                    ! candidates. As this test is relatively cheap, do it
                    ! now for this bounding box.
              
                    d1 = abs(coor(1) - xBBox(1,kk))
                    d2 = abs(coor(1) - xBBox(4,kk))
                    dx = max(d1,d2)

                    d1 = abs(coor(2) - xBBox(2,kk))
                    d2 = abs(coor(2) - xBBox(5,kk))
                    dy = max(d1,d2)

                    d1 = abs(coor(3) - xBBox(3,kk))
                    d2 = abs(coor(3) - xBBox(6,kk))
                    dz = max(d1,d2)

                    d2 = dx*dx + dy*dy + dz*dz
                    uvw(4) = min(uvw(4),d2)

                  endif testStoreBBox

                else terminalTest

                  ! Child contains a leaf. Compute the possible minimum
                  ! distance squared to the current coordinate.

                  if(     coor(1) < ADTree(kk)%xMin(1)) then
                    dx =  coor(1) - ADTree(kk)%xMin(1)
                  else if(coor(1) > ADTree(kk)%xMax(4)) then
                    dx =  coor(1) - ADTree(kk)%xMax(4)
                  else
                    dx = adtZero
                  endif

                  if(     coor(2) < ADTree(kk)%xMin(2)) then
                    dy =  coor(2) - ADTree(kk)%xMin(2)
                  else if(coor(2) > ADTree(kk)%xMax(5)) then
                    dy =  coor(2) - ADTree(kk)%xMax(5)
                  else
                    dy = adtZero
                  endif

                  if(     coor(3) < ADTree(kk)%xMin(3)) then
                    dz =  coor(3) - ADTree(kk)%xMin(3)
                  else if(coor(3) > ADTree(kk)%xMax(6)) then
                    dz =  coor(3) - ADTree(kk)%xMax(6)
                  else
                    dz = adtZero
                  endif

                  d2 = dx*dx + dy*dy + dz*dz

                  ! If this distance squared is less than the current
                  ! value, store this leaf in the new front.

                  testStoreLeave: if(d2 < uvw(4)) then

                    ! Check if enough memory has been allocated and
                    ! store the leaf.

                    if(nFrontLeavesNew == nAllocFront) then
                      i = nAllocFront
                      call reallocPlus(frontLeavesNew, i, 25, ADT)
                      call reallocPlus(frontLeaves, nAllocFront, 25, ADT)
                    endif

                    nFrontLeavesNew = nFrontLeavesNew + 1
                    frontLeavesNew(nFrontLeavesNew) = kk

                    ! Compute the guaranteed distance squared to this leaf.
                    ! It may be less than the currently stored value.

                    d1 = abs(coor(1) - ADTree(kk)%xMin(1))
                    d2 = abs(coor(1) - ADTree(kk)%xMax(4))
                    dx = max(d1,d2)

                    d1 = abs(coor(2) - ADTree(kk)%xMin(2))
                    d2 = abs(coor(2) - ADTree(kk)%xMax(5))
                    dy = max(d1,d2)

                    d1 = abs(coor(3) - ADTree(kk)%xMin(3))
                    d2 = abs(coor(3) - ADTree(kk)%xMax(6))
                    dz = max(d1,d2)

                    d2 = dx*dx + dy*dy + dz*dz
                    uvw(4) = min(uvw(4),d2)

                  endif testStoreLeave

                endif terminalTest
              enddo childrenLoop
            enddo currentFrontLoop

            ! End of the loop over the current front. If the new front
            ! is empty the entire tree has been traversed and an exit is
            ! made from the corresponding loop.

            if(nFrontLeavesNew == 0) exit treeTraversal2

            ! Copy the data of the new front leaves into the current
            ! front for the next round.

            nFrontLeaves = nFrontLeavesNew
            do ll=1,nFrontLeaves
              frontLeaves(ll) = frontLeavesNew(ll)
            enddo

          enddo treeTraversal2

          ! Sort the target bounding boxes in increasing order such that
          ! the one with the smallest possible distance is first.

          call qsortBBoxTargets(BB, nBB, ADT)
!
!         **************************************************************
!         *                                                            *
!         * Part 4: Loop over the selected bounding boxes and check if *
!         *         the corresponding element minimizes the distance.  *
!         *                                                            *
!         **************************************************************
!
          elementFound = .false.

          BBoxLoop: do mm=1,nBB

            ! Exit the loop if the possible minimum distance of this
            ! bounding box is not smaller than the current value.
            ! Remember that BB has been sorted in increasing order.

            if(uvw(4) <= BB(mm)%posDist2) exit BBoxLoop

            ! Determine the element type stored in this bounding box.

            kk = BB(mm)%ID
            select case (ADT%elementType(kk))

              case (adtTriangle)
                call adtTerminate(ADT, "minDistanceTreeSearch", &
                                  "Minimum distance search for &
                                  &triangles not implemented yet")

              !=========================================================

              case (adtQuadrilateral)

                ! Temporary implementation. I'm waiting for Juan to come
                ! up with his more sophisticated algorithm.

                ! This is a surface element, so set the parametric weight
                ! w to zero.

                w = adtZero

                ! Determine the 4 vectors which completely describe
                ! the quadrilateral face

                ll   = ADT%elementID(kk)
                n(1) = ADT%quadsConn(1,ll)
                n(2) = ADT%quadsConn(2,ll)
                n(3) = ADT%quadsConn(3,ll)
                n(4) = ADT%quadsConn(4,ll)

                x1(1) = ADT%coor(1,n(1))
                x1(2) = ADT%coor(2,n(1))
                x1(3) = ADT%coor(3,n(1))

                x21(1) = ADT%coor(1,n(2)) - x1(1)
                x21(2) = ADT%coor(2,n(2)) - x1(2)
                x21(3) = ADT%coor(3,n(2)) - x1(3)

                x41(1) = ADT%coor(1,n(4)) - x1(1)
                x41(2) = ADT%coor(2,n(4)) - x1(2)
                x41(3) = ADT%coor(3,n(4)) - x1(3)

                x3142(1) = ADT%coor(1,n(3)) - x1(1) - x21(1) - x41(1)
                x3142(2) = ADT%coor(2,n(3)) - x1(2) - x21(2) - x41(2)
                x3142(3) = ADT%coor(3,n(3)) - x1(3) - x21(3) - x41(3)

                ! Initialize u and v to 0.5 and determine the
                ! corresponding coordinates on the face, which is the
                ! centroid.

                u  = adtHalf
                v  = adtHalf
                uv = u*v

                xf(1) = x1(1) + u*x21(1) + v*x41(1) + uv*x3142(1)
                xf(2) = x1(2) + u*x21(2) + v*x41(2) + uv*x3142(2)
                xf(3) = x1(3) + u*x21(3) + v*x41(3) + uv*x3142(3)

                ! Newton loop to determine the point on the surface,
                ! which minimizes the distance to the given coordinate.

                NewtonQuads: do ll=1,iterMax

                  ! Store the current values of u and v for a stop
                  ! criterion later on.

                  uold = u
                  vold = v

                  ! Determine the vector vf from xf to given coordinate.

                  vf(1) = coor(1) - xf(1)
                  vf(2) = coor(2) - xf(2)
                  vf(3) = coor(3) - xf(3)

                  ! Determine the tangent vectors in u- and v-direction.
                  ! Store these in a and b respectively.

                  a(1) = x21(1) + v*x3142(1)
                  a(2) = x21(2) + v*x3142(2)
                  a(3) = x21(3) + v*x3142(3)

                  b(1) = x41(1) + u*x3142(1)
                  b(2) = x41(2) + u*x3142(2)
                  b(3) = x41(3) + u*x3142(3)

                  ! Determine the normal vector of the face by taking the
                  ! cross product of a and b. Afterwards this vector will
                  ! be scaled to a unit vector.

                  norm(1) = a(2)*b(3) - a(3)*b(2)
                  norm(2) = a(3)*b(1) - a(1)*b(3)
                  norm(3) = a(1)*b(2) - a(2)*b(1)

                  invLen = adtOne/max(adtEps,sqrt(norm(1)*norm(1) &
                         +                        norm(2)*norm(2) &
                         +                        norm(3)*norm(3)))

                  norm(1) = norm(1)*invLen
                  norm(2) = norm(2)*invLen
                  norm(3) = norm(3)*invLen

                  ! Determine the projection of the vector vf onto
                  ! the face.

                  vn = vf(1)*norm(1) + vf(2)*norm(2) + vf(3)*norm(3)
                  vt(1) = vf(1) - vn*norm(1)
                  vt(2) = vf(2) - vn*norm(2)
                  vt(3) = vf(3) - vn*norm(3)

                  ! The vector vt points from the current point on the
                  ! face to the new point. However this new point lies on
                  ! the plane determined by the vectors a and b, but not
                  ! necessarily on the face itself. The new point on the
                  ! face is obtained by projecting the point in the a-b
                  ! plane onto the face. this can be done by determining
                  ! the coefficients du and dv, such that vt = du*a + dv*b.
                  ! To solve du and dv the vectors normal to a and b
                  ! inside the plane ab are needed.

                  an(1) = a(2)*norm(3) - a(3)*norm(2)
                  an(2) = a(3)*norm(1) - a(1)*norm(3)
                  an(3) = a(1)*norm(2) - a(2)*norm(1)

                  bn(1) = b(2)*norm(3) - b(3)*norm(2)
                  bn(2) = b(3)*norm(1) - b(1)*norm(3)
                  bn(3) = b(1)*norm(2) - b(2)*norm(1)

                  ! Solve du and dv. the clipping of vn should not be
                  ! active, as this would mean that the vectors a and b
                  ! are parallel. This corresponds to a quad degenerated
                  ! to a line, which should not occur in the surface mesh.

                  vn = a(1)*bn(1) + a(2)*bn(2) + a(3)*bn(3)
                  vn = sign(max(adtEps,abs(vn)),vn)
                  du = (vt(1)*bn(1) + vt(2)*bn(2) + vt(3)*bn(3))/vn

                  vn = b(1)*an(1) + b(2)*an(2) + b(3)*an(3)
                  vn = sign(max(adtEps,abs(vn)),vn)
                  dv = (vt(1)*an(1) + vt(2)*an(2) + vt(3)*an(3))/vn

                  ! Determine the new parameter values uu and vv. These
                  ! are limited to 0 <= (uu,vv) <= 1.

                  u = u + du; u = min(adtOne,max(adtZero,u))
                  v = v + dv; v = min(adtOne,max(adtZero,v))

                  ! Determine the final values of the corrections.

                  du = abs(u-uold)
                  dv = abs(v-vold)

                  ! Determine the new coordinates of the point xf.

                  uv  = u*v
                  xf(1) = x1(1) + u*x21(1) + v*x41(1) + uv*x3142(1)
                  xf(2) = x1(2) + u*x21(2) + v*x41(2) + uv*x3142(2)
                  xf(3) = x1(3) + u*x21(3) + v*x41(3) + uv*x3142(3)

                  ! Exit the loop if the update of the parametric
                  ! weights is below the threshold

                  val = sqrt(du*du + dv*dv)
                  if(val <= thresConv) exit NewtonQuads

                enddo NewtonQuads

                ! Compute the distance squared between the given
                ! coordinate and the point xf.

                dx = coor(1) - xf(1)
                dy = coor(2) - xf(2)
                dz = coor(3) - xf(3)

                val = dx*dx + dy*dy + dz*dz

                ! If the distance squared is less than the current value
                ! store the wall distance and interpolation info and
                ! indicate that an element was found.

                if(val < uvw(4)) then
                  uvw(4)    = val
                  nNodeElement = 4
                  elementFound = .true.

                  kkk  = kk;   uu   = u;    vv   = v;    ww   = w
                  m(1) = n(1); m(2) = n(2); m(3) = n(3); m(4) = n(4)

                  weight(1) = (adtOne - u)*(adtOne - v)
                  weight(2) =           u *(adtOne - v)
                  weight(3) =           u *          v
                  weight(4) = (adtOne - u)*          v
                endif

              !=========================================================

              case (adtTetrahedron)
                call adtTerminate(ADT, "minDistanceTreeSearch", &
                                  "Minimum distance search for &
                                  &tetrahedra not implemented yet")

              !===========================================================

              case (adtPyramid)
                call adtTerminate(ADT, "minDistanceTreeSearch", &
                                  "Minimum distance search for &
                                  &pyramids not implemented yet")

              !===========================================================

              case (adtPrism)
                call adtTerminate(ADT, "minDistanceTreeSearch", &
                                  "Minimum distance search for &
                                  &prisms not implemented yet")

              !===========================================================

              case (adtHexahedron)

                ! Determine the element ID and the corresponding
                ! 8 node ID's.

                ll = ADT%elementID(kk)
                n(1) = ADT%hexaConn(1,ll)
                n(2) = ADT%hexaConn(2,ll)
                n(3) = ADT%hexaConn(3,ll)
                n(4) = ADT%hexaConn(4,ll)
                n(5) = ADT%hexaConn(5,ll)
                n(6) = ADT%hexaConn(6,ll)
                n(7) = ADT%hexaConn(7,ll)
                n(8) = ADT%hexaConn(8,ll)

                ! Call the subroutine minD2Hexa to do the work.

                call minD2Hexa(coor(1:3),            &
                               ADT%coor(1:3,n(1)), &
                               ADT%coor(1:3,n(2)), &
                               ADT%coor(1:3,n(3)), &
                               ADT%coor(1:3,n(4)), &
                               ADT%coor(1:3,n(5)), &
                               ADT%coor(1:3,n(6)), &
                               ADT%coor(1:3,n(7)), &
                               ADT%coor(1:3,n(8)), &
                               val, chi, ierr)

                ! If the distance squared is less than the current value
                ! store the wall distance and interpolation info and
                ! indicate that an element was found.

                if(val < uvw(4)) then
                  uvw(4)    = val
                  nNodeElement = 8
                  elementFound = .true.

                  kkk = kk;
                  uu  = chi(1); vv = chi(2); ww = chi(3)

                  m(1) = n(1); m(2) = n(2); m(3) = n(3); m(4) = n(4)
                  m(5) = n(5); m(6) = n(6); m(7) = n(7); m(8) = n(8)

                  weight(1) = (adtOne - uu)*(adtOne - vv)*(adtOne - ww)
                  weight(2) =           uu *(adtOne - vv)*(adtOne - ww)
                  weight(3) =           uu *          vv *(adtOne - ww)
                  weight(4) = (adtOne - uu)*          vv *(adtOne - ww)
                  weight(5) = (adtOne - uu)*(adtOne - vv)*          ww
                  weight(6) =           uu *(adtOne - vv)*          ww
                  weight(7) =           uu *          vv *          ww
                  weight(8) = (adtOne - uu)*          vv *          ww
                endif

            end select

          enddo BBoxLoop

          ! Check if an element was found. As all the minimum distance
          ! searches are initialized by the calling routine (to support
          ! periodic searches) this is not always the case.

          if( elementFound ) then

            ! Store the interpolation information for this point.
            ! First the integer info, i.e. the processor ID, element type
            ! and local element ID.

            intInfo(1) = ADT%myID
            intInfo(2) = ADT%elementType(kkk)
            intInfo(3) = ADT%elementID(kkk)

            ! The parametric weights. Note that the wall distance
            ! squared, stored in the 4th position of uvw, already
            ! contains the correct value.

            uvw(1) = uu
            uvw(2) = vv
            uvw(3) = ww

            ! The interpolated solution, if needed.

            do ll=1,nInterpol
              ii = 4+ll
              uvw(ii) = weight(1)*arrDonor(ll,m(1))
              do i=2, nNodeElement
                uvw(ii) = uvw(ii) + weight(i)*arrDonor(ll,m(i))
              enddo
            enddo

          endif

        end subroutine minDistanceTreeSearchSinglePoint

        subroutine intersectionTreeSearchSinglePoint(ADT, coor, &
             intInfo, BB, frontLeaves, frontLeavesNew)
!
!       ****************************************************************
!       *                                                              *
!       * This routine is used in the ray casting approach to determine*
!       * if a given ray intersects any of the surface elements. The   *
!       * purpose is to determine if a point of interest is inside     *
!       * or outside the (closed) surface defined by the ADT.          *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * ADT:       ADT type whose ADT must be searched               *
!       * coor(3):   The coordinate of the point to be searched.       *
!       *                                                              *
!       * Subroutine intent(out) arguments.                            *
!       * ---------------------------------                            *
!       * intInfo: Intersection info. The number of intersections we   *
!       *          found. 
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments.
!
        type(adtType), intent(inout) :: ADT

        real(kind=realType), dimension(3), intent(in) :: coor
        integer(kind=intType), intent(out) :: intInfo

        integer(kind=intType), dimension(:), pointer :: BB
        integer(kind=intType), dimension(:), pointer :: frontLeaves
        integer(kind=intType), dimension(:), pointer :: frontLeavesNew
!
!       Local variables.
!
        integer :: ierr

        integer(kind=intType) :: ii, kk, ll, mm, nn
        integer(kind=intType) :: nBB, nFrontLeaves, nFrontLeavesNew
        integer(kind=intType) :: nAllocBB, nAllocFront
        integer(kind=intType) :: i, nNodeElement
        real(kind=realType), dimension(:,:), pointer :: xBBox
        logical :: elementFound
        type(adtLeafType), dimension(:), pointer :: ADTree
!
!       ****************************************************************
!       *                                                              *
!       * Begin execution.                                             *
!       *                                                              *
!       ****************************************************************
!
        ! Set some pointers to make the code more readable.

        xBBox  => ADT%xBBox
        ADTree => ADT%ADTree

        ! Determine the sizes from the arrays we have been passed

        nAllocBB    = size(BB)
        nAllocFront = size(frontLeaves)
        
        ! Initialize the number of possible intersections to 0

        intInfo = 0
!
!         **************************************************************
!         *                                                            *
!         * Part 1. Traverse the tree and determine the target         *
!         *         bounding boxes, which may contain the intersection *
!         *                                                            *
!         **************************************************************
!
          ! Start at the root, i.e. set the front leaf to the root leaf.
          ! Also initialize the number of possible bounding boxes to 0.

          nBB = 0

          nFrontLeaves   = 1
          frontLeaves(1) = 1

          treeTraversalLoop: do

            ! Initialize the number of leaves for the new front, i.e.
            ! the front of the next round, to 0.

            nFrontLeavesNew = 0

            ! Loop over the leaves of the current front.

            currentFrontLoop: do ii=1,nFrontLeaves

              ! Store the ID of the leaf a bit easier and loop over
              ! its two children.

              ll = frontLeaves(ii)

              childrenLoop: do mm=1,2

                ! Determine whether this child contains a bounding box
                ! or a leaf of the next level.

                kk = ADTree(ll)%children(mm)
                terminalTest: if(kk < 0) then

                  ! Child contains a bounding box. Check if the
                  ! ray is inside the bounding box.

                  kk = -kk
                  if(  coor(1) <= xBBox(4,kk) .and. &

                       coor(2) >= xBBox(2,kk) .and. &
                       coor(2) <= xBBox(5,kk) .and. &
                       coor(3) >= xBBox(3,kk) .and. &
                       coor(3) <= xBBox(6,kk)) then

                    ! Ray intersectst he bounding box. That's all we
                    ! wanted to know:
                     intInfo = 1
                     exit treeTraversalLoop
                  end if
                else terminalTest

                  ! Child contains a leaf. Check if the coordinate is
                  ! inside the bounding box of the leaf.

                  if(coor(1) <= ADTree(kk)%xMax(4) .and. &

                     coor(2) >= ADTree(kk)%xMin(2) .and. &
                     coor(2) <= ADTree(kk)%xMax(5) .and. &
                     coor(3) >= ADTree(kk)%xMin(3) .and. &
                     coor(3) <= ADTree(kk)%xMax(6)) then

                    ! Coordinate is inside the leaf. Store the leaf in
                    ! the list for the new front.

                    if(nFrontLeavesNew == nAllocFront) then
                      i = nAllocFront
                      call reallocPlus(frontLeavesNew, i, 25, ADT)
                      call reallocPlus(frontLeaves, nAllocFront, 25, ADT)
                    endif

                    nFrontLeavesNew = nFrontLeavesNew + 1
                    frontLeavesNew(nFrontLeavesNew) = kk
                  end if

                endif terminalTest

              enddo childrenLoop

            enddo currentFrontLoop

            ! End of the loop over the current front. If the new front
            ! is empty the entire tree has been traversed and an exit is
            ! made from the corresponding loop.

            if(nFrontLeavesNew == 0) then 
               exit treeTraversalLoop
            end if
            ! Copy the data of the new front leaves into the current
            ! front for the next round.

            nFrontLeaves = nFrontLeavesNew
            do ll=1,nFrontLeaves
              frontLeaves(ll) = frontLeavesNew(ll)
            enddo

          enddo treeTraversalLoop
!
!         **************************************************************
!         *                                                            *
!         * Part 2: Loop over the selected bounding boxes and check if *
!         *         the corresponding elements contain the point.      *
!         *                                                            *
!         **************************************************************
!
          !intInfo = nBB
          ! elementFound = .false.

          ! BBoxLoop: do mm=1,nBB

          !   ! Determine the element type stored in this bounding box.

          !   kk = BB(mm)
          !   select case (ADT%elementType(kk))

          !     case (adtQuadrilateral)


          !     !=========================================================

          !     case (adtTriangle)

         

          !   end select
     
          ! enddo BBoxLoop

        end subroutine intersectionTreeSearchSinglePoint


      end module adtLocalSearch
