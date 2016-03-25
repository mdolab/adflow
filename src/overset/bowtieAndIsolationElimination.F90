subroutine bowTieAndIsolationElimination(level, sps)

  use blockPointers
  use BCTypes
  use communication
  implicit none

  ! Input Parameters
  integer(kind=intType), intent(in) :: level, sps

  ! Local variables
  integer(kind=intType) :: mm, nn, i, j, k, e, iBeg, iEnd, jBeg, jEnd
  logical :: side(4), isWallType

  integer(kind=intType), dimension(:, :), pointer :: ibp, gcp
  integer(kind=intType), dimension(:, :), allocatable :: toFlip, nE, nC

  ! This routine initializes the surface cell iblank based on the
  ! volume iblank. It is not a straight copy since we a little
  ! preprocessing to eliminate a few particularly nasty cases.
  ! Three analysis are performed:
  ! 1. Bow-tie elimination
  ! 2. Single cell elmination

  bowTieLoop: do E=0, 2
     domainLoop1: do nn=1, nDom
        call setPointers(nn, level, sps)
        
        bocoLoop1: do mm=1, nBocos
           wallType1: if (isWallType(BCType(mm))) then
              
              select case (BCFaceID(mm))
              case (iMin)
                 ibp => iblank(2, :, :)
                 gcp => globalCell(2, :, :)
              case (iMax)
                 ibp => iblank(il, :, :)
                 gcp => globalCell(il, :, :)
              case (jMin)
                 ibp => iblank(:, 2, :)
                 gcp => globalCell(:, 2, :)
              case (jMax)
                 ibp => iblank(:, jl, :)
                 gcp => globalCell(:, jl, :)
              case (kMin)
                 ibp => iblank(:, :, 2)
                 gcp => globalCell(:, :, 2)
              case (kMax)
                 ibp => iblank(:, :, kl)
                 gcp => globalCell(:, :, kl)
              end select
              
              ! -------------------------------------------------
              ! Step 2: Bow-tie elimination: Elimiate cells
              ! that touch only at a corner.
              ! -------------------------------------------------
              
              ! Make bounds a little easier to read. Owned cells only
              ! from now on.
              jBeg = BCData(mm)%jnBeg+1 ; jEnd = BCData(mm)%jnEnd
              iBeg = BCData(mm)%inBeg+1 ; iEnd = BCData(mm)%inEnd

              ! Allocate two tmporary auxilary arrays 'eN'->
              ! edgeNeighbours and 'cN'-> cornerNeighbous. For every
              ! comute determine the number of compute neighbours
              ! connected along edges and at corners
              allocate(nE(iBeg:iEnd, jBeg:jEnd), nC(iBeg:iEnd, jBeg:jEnd))!, &
              
              call findBowTies()
              
              do j=jBeg, jEnd
                 do i=iBeg, iEnd
                    if (BCData(mm)%iBlank(i, j) > 0 .and. nC(i,j) >=1 .and. nE(i,j) <=E) then
                       BCData(mm)%iBlank(i, j) = 0
                    end if
                 end do
              end do
              
              deallocate(nC, nE)
           end if wallType1
        end do bocoLoop1
     end do domainLoop1

     ! Since we potentially changed iBlanks, we need to updated by
     ! performing an exchange.
     call exchangeSurfaceIBlanks(level, sps, commPatternCell_2nd, internalCell_2nd)
     
  end do bowTieLoop
     
  domainLoop2: do nn=1, nDom
     call setPointers(nn, level, sps)

     bocoLoop2: do mm=1, nBocos
        wallType2: if (isWallType(BCType(mm))) then

           select case (BCFaceID(mm))
           case (iMin)
              gcp => globalCell(2, :, :)
           case (iMax)
              gcp => globalCell(il, :, :)
           case (jMin)
              gcp => globalCell(:, 2, :)
           case (jMax)
              gcp => globalCell(:, jl, :)
           case (kMin)
              gcp => globalCell(:, :, 2)
           case (kMax)
              gcp => globalCell(:, :, kl)
           end select

           ! Make bounds a little easier to read. Owned cells only
           ! from now on.
           jBeg = BCData(mm)%jnBeg+1 ; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg+1 ; iEnd = BCData(mm)%inEnd

           ! -------------------------------------------------
           ! Step 3: Single-cell elimination: Elimiate cells
           ! that do not touch any other cells.
           ! -------------------------------------------------

           allocate(nE(iBeg:iEnd, jBeg:jEnd), nC(iBeg:iEnd, jBeg:jEnd))

           call setNeighbourCounts()

           ! This is easy, if a compute cell is stil around with no
           ! neighbours, kill it
           do j=jBeg, jEnd
              do i=iBeg, iEnd
                 if (BCData(mm)%iBlank(i, j) == 1 .and. nE(i,j) == 0 .and. nC(i,j) == 0) then
                    BCData(mm)%iBlank(i, j) = 0
                 end if
              end do
           end do

           deallocate(nE, nC)

        end if wallType2
     end do bocoLoop2
  end do domainLoop2

  ! Again, since we potentially changed iBlanks, we need to updated by
  ! performing an exchange.
  call exchangeSurfaceIBlanks(level, sps, commPatternCell_2nd, internalCell_2nd)
contains
  subroutine findBowTies

    implicit none
    ! For every compute determine the number of compute neighbours
    ! connected along edges (nE) and the number of bow-ties (in nC)

    integer(kind=intType) :: i, j, e(4)

    nE = 0
    nC = 0

    do j=jBeg, jEnd
       do i=iBeg, iEnd
          if (BCData(mm)%iBlank(i, j) >= 0 ) then
             e = 0

             !    |  e3  |
             !  --c4-----c3-
             !    |      |
             ! e4 |   x  | e2
             !    |      |
             ! -- c1-----c2-
             !    |  e1  |

             ! Set the status of each of the 4 edges:

             if (BCData(mm)%iBlank(i, j-1) == 1) &
                  e(1) = 1

             if (BCData(mm)%iBlank(i+1, j) == 1) &
                  e(2) = 1

             if (BCData(mm)%iBlank(i, j+1) == 1) &
                  e(3) = 1

             if (BCData(mm)%iBlank(i-1, j) == 1) &
                  e(4) = 1

             ! Check the 4 corner neighbours for bow-tie status
             if (BCData(mm)%iBlank(i-1, j-1) == 1 .and. e(4) == 0 .and. e(1) == 0) &
                  nC(i, j) = nC(i, j) + 1

             if (BCData(mm)%iBlank(i+1, j-1) == 1 .and. e(1) == 0 .and. e(2) == 0) &
                  nC(i, j) = nC(i, j) + 1

             if (BCData(mm)%iBlank(i+1, j+1) == 1 .and. e(2) == 0 .and. e(3) == 0) &
                  nC(i, j) = nC(i, j) + 1

             if (BCData(mm)%iBlank(i-1, j+1) == 1 .and. e(3) == 0 .and. e(4) == 0) &
                  nC(i, j) = nC(i, j) + 1
             
             nE(i, j) = sum(e)
          end if
       end do
    end do
  end subroutine findBowTies

 subroutine setNeighbourCounts

    implicit none
    ! For every comute determine the number of compute neighbours
    ! connected along edges and at corners

    integer(kind=intType) :: i, j

    nE = 0
    nC = 0

    do j=jBeg, jEnd
       do i=iBeg, iEnd
          if (BCData(mm)%iBlank(i, j) == 1) then

             !    |  e3  |
             !  --c4-----c3-
             !    |      |
             ! e4 |   x  | e2
             !    |      |
             ! -- c1-----c2-
             !    |  e1  |

             ! Set the status of each of the 4 edges:

             if (BCData(mm)%iBlank(i, j-1) == 1) &
                  nE(i, j) = nE(i, j) + 1

             if (BCData(mm)%iBlank(i+1, j) == 1) &
                  nE(i, j) = nE(i, j) + 1

             if (BCData(mm)%iBlank(i, j+1) == 1) &
                  nE(i, j) = nE(i, j) + 1

             if (BCData(mm)%iBlank(i-1, j) == 1) &
                  nE(i, j) = nE(i, j) + 1

             ! Check the 4 corner neighbours for compute neighbour
             if (BCData(mm)%iBlank(i-1, j-1) == 1) &
                  nC(i, j) = nC(i, j) + 1

             if (BCData(mm)%iBlank(i+1, j-1) == 1) &
                  nC(i, j) = nC(i, j) + 1

             if (BCData(mm)%iBlank(i+1, j+1) == 1) &
                  nC(i, j) = nC(i, j) + 1

             if (BCData(mm)%iBlank(i-1, j+1) == 1) &
                  nC(i, j) = nC(i, j) + 1
          end if
       end do
    end do
  end subroutine setNeighbourCounts
end subroutine bowTieAndIsolationElimination
