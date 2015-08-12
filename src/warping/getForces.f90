subroutine getNPatches(nPatches)

  ! Return the number of wall patches we have on this proc

  use BCTypes
  use blockPointers
 
  implicit none

  ! Output Variables
  integer(kind=intType),intent(out) :: nPatches

  ! Working Variables
  integer(kind=intType) :: nn, mm

  nPatches = 0_intType

  domains: do nn=1,nDom
     call setPointers(nn,1_intType,1_intType)

     bocos: do mm=1,nBocos
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or.&
             BCType(mm) == NSWallIsothermal) then
           nPatches = nPatches + 1
        end if
     end do bocos
  end do domains
end subroutine getNPatches

subroutine getPatchName(iPatch, patchName) 
  ! Get names one at a time since f2py can't do arrays of strings (nicely)
  use BCTypes
  use blockPointers
  use cgnsGrid

  implicit none

  ! Input Variables
  integer(kind=intType),intent(in) :: iPatch

  ! Output Variables 
  character*(*), intent(inout) :: patchName

  ! Working
  integer(kind=intType) :: nn, mm, patchCount, cgb

  patchCount = 0
  domains: do nn=1,nDom
     call setPointers(nn,1_intType,1_intType)

     bocos: do mm=1,nBocos
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or.&
             BCType(mm) == NSWallIsothermal) then
           cgb = flowDoms(nn,1,1)%cgnsBlockID
           if (patchCount + 1 == iPatch) then
              patchName = cgnsDoms(cgb)%bocoInfo(cgnsSubface(mm))%wallBCName
           end if
           patchCount = patchCount + 1
        end if
     end do bocos
  end do domains
end subroutine getPatchName

subroutine getPatchSize(iPatch, patchSize)
  ! Get size of one patach
  use BCTypes
  use blockPointers
  use cgnsGrid

  implicit none

  ! Input Variables
  integer(kind=intType),intent(in) :: iPatch

  ! Output Variables 
  integer(kind=intType), intent(out) :: patchSize(2)

  ! Working
  integer(kind=intType) :: nn, mm, patchCount, iBeg, iEnd, jBeg, jEnd

  patchCount = 0
  domains: do nn=1,nDom
     call setPointers(nn,1_intType,1_intType)

     bocos: do mm=1,nBocos
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or.&
             BCType(mm) == NSWallIsothermal) then
           if (patchCount + 1 == iPatch) then
              jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
              iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd
              patchSize(1) = (iEnd - iBeg + 1)
              patchSize(2) = (jEnd - jBeg + 1)
           end if
           patchCount = patchCount + 1
        end if
     end do bocos
  end do domains
end subroutine getPatchSize

subroutine getForceSize(size, sizeCell)
  ! Compute the number of points that will be returned from getForces
  ! or getForcePoints
  use BCTypes
  use blockPointers
  use inputTimeSpectral
  implicit none

  integer(kind=intType),intent(out) :: size, sizeCell
  integer(kind=intType) :: nn,mm
  integer(kind=intType) :: iBeg,iEnd,jBeg,jEnd

  size = 0_intType
  sizeCell = 0_intType
 
  domains: do nn=1,nDom
     call setPointers(nn,1_intType,1_intType)
     bocos: do mm=1,nBocos
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or.&
             BCType(mm) == NSWallIsothermal) then

           jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd
           size = size + (iEnd - iBeg + 1)*(jEnd - jBeg + 1)
           sizeCell = sizeCell + (iEnd - iBeg)*(jEnd - jBeg)
        end if
     end do bocos
  end do domains
end subroutine getForceSize

subroutine getForceConnectivity(conn, ncell)
  ! Return the connectivity list for the each of the patches

  use BCTypes
  use blockPointers
  use inputPhysics
  implicit none

  ! Input/Output
  integer(kind=intType), intent(in) :: ncell
  integer(kind=intType), intent(inout) :: conn(4*ncell)

  ! Working
  integer(kind=intType) :: nn, mm, cellCount, nodeCount, ni, nj, i, j
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
  logical regularOrdering
  cellCount = 0
  nodeCount = 0

  domains: do nn=1,nDom
     call setPointers(nn, 1_intType, 1_intType)
     bocos: do mm=1,nBocos
        if(BCType(mm) == EulerWall .or. &
           BCType(mm) == NSWallAdiabatic .or. &
           BCType(mm) == NSWallIsothermal) then

           jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd
           
           ni = iEnd - iBeg + 1
           nj = jEnd - jBeg + 1
           
           ! We want to ensure that all the normals of the faces are
           ! consistent. To ensure this, we enforce that all normals
           ! are "into" the domain. Therefore we must treat difference
           ! faces of a block differently. For example for an iLow
           ! face, when looping over j-k in the regular way, results
           ! in in a domain inward pointing normal for iLow but
           ! outward pointing normal for iHigh. The same is true for
           ! kMin and kMax. However, it is reverse for the J-faces:
           ! This is becuase the way the pointers are extracted i then
           ! k is the reverse of what "should" be for consistency. The
           ! other two, the pointers are cyclic consistent: i,j->k,
           ! j,k (wrap) ->i, but for the j-direction is is i,k->j when
           ! to be consistent with the others it should be
           ! k,i->j. Hope that made sense. 

           select case(BCFaceID(mm))
           case(iMin, jMax, kMin)
              regularOrdering = .True.
           case default
              regularOrdering = .False.
           end select

           ! Now this can be reversed *again* if we have a block that
           ! is left handed. 
           if (.not. rightHanded) then 
              regularOrdering = .not. (regularOrdering)
           end if

           if (regularOrdering) then 
                 ! Do regular ordering.
                 
                 ! Loop over generic face size...Note we are doing zero
                 ! based ordering!

                 ! This cartoon of a generic cell might help:
                 !
                 ! i, j+1 +-----+ i+1, j+1
                 !   n4   |     | n3
                 !        +-----+
                 !       i,j    i+1, j
                 !       n1     n2
                 !

                 do j=0,nj-2
                    do i=0,ni-2
                       conn(4*cellCount+1) = nodeCount + (j  )*ni + i     ! n1
                       conn(4*cellCount+2) = nodeCount + (j  )*ni + i + 1 ! n2
                       conn(4*cellCount+3) = nodeCount + (j+1)*ni + i + 1 ! n3
                       conn(4*cellCount+4) = nodeCount + (j+1)*ni + i     ! n4
                       cellCount = cellCount + 1
                    end do
                 end do
              else
                 ! Do reverse ordering:
                 do j=0,nj-2
                    do i=0,ni-2
                       conn(4*cellCount+1) = nodeCount + (j  )*ni + i     ! n1
                       conn(4*cellCount+2) = nodeCount + (j+1)*ni + i     ! n4
                       conn(4*cellCount+3) = nodeCount + (j+1)*ni + i + 1 ! n3
                       conn(4*cellCount+4) = nodeCount + (j  )*ni + i + 1 ! n2
                       cellCount = cellCount + 1
                    end do
                 end do
              end if
           nodeCount = nodeCount + ni*nj
        end if
     end do bocos
  end do domains
end subroutine getForceConnectivity

subroutine getForcePoints(points, npts, sps_in)

  use BCTypes
  use blockPointers
  use inputTimeSpectral
  implicit none
  !
  !      Local variables.
  !
  integer(kind=intType), intent(in) :: npts,sps_in
  real(kind=realType), intent(inout) :: points(3,npts)

  integer(kind=intType) :: mm, nn, i, j, ii,sps
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd

  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************

  sps = sps_in

  ii = 0 
  domains: do nn=1,nDom
     call setPointers(nn, 1_intType, sps)
     
     ! Loop over the number of boundary subfaces of this block.
     bocos: do mm=1,nBocos
        
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or.&
             BCType(mm) == NSWallIsothermal) then
           
           ! NODE Based
           jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd
           
           do j=jBeg, jEnd ! This is a node loop
              do i=iBeg, iEnd ! This is a node loop
                 ii = ii +1
                 select case(BCFaceID(mm))
                    
                 case(imin)
                    points(:,ii) = x(1,i,j,:)
                 case(imax)
                    points(:,ii) = x(il,i,j,:)
                 case(jmin) 
                    points(:,ii) = x(i,1,j,:)
                 case(jmax) 
                    points(:,ii) = x(i,jl,j,:)
                 case(kmin) 
                    points(:,ii) = x(i,j,1,:)
                 case(kmax) 
                    points(:,ii) = x(i,j,kl,:)
                 end select    
              end do
           end do
        end if
     end do bocos
  end do domains
  
end subroutine getForcePoints

subroutine getForces(forces, npts, sps_in)

  use BCTypes
  use blockPointers
  use flowVarRefState
  use inputTimeSpectral
  use communication
  use inputPhysics
  implicit none
  !
  !      Local variables.
  !
  integer(kind=intType), intent(in) :: npts, sps_in
  real(kind=realType), intent(out) :: forces(3,npts)

    real(kind=realType) :: area(npts) ! Dual area's
  integer(kind=intType) :: mm, nn, i, j, ii, jj,sps
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
  real(kind=realType) :: sss(3),v2(3),v1(3), qa, sepSensor, Cavitation
  integer(kind=intType) :: lower_left,lower_right,upper_left,upper_right
  real(kind=realType) :: cFp(3), cFv(3), cMp(3), cMv(3), yplusmax, qf(3)

  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************

  sps = sps_in

  ii = 0 
  domains: do nn=1,nDom
     call setPointers(nn,1_intType,sps)
     call forcesAndMoments(cFp, cFv, cMp, cMv, yplusMax, sepSensor, Cavitation)
     
     ! Loop over the number of boundary subfaces of this block.
     bocos: do mm=1,nBocos
        
        if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or. &
             BCType(mm) == NSWallIsothermal) then

           ! This is easy, just copy out F in continuous ordering. 
           do j=BCData(mm)%jnBeg,BCData(mm)%jnEnd
              do i=BCData(mm)%inBeg,BCData(mm)%inEnd
                 ii = ii + 1
                 Forces(:, ii) = bcData(mm)%F(i, j, :)
              end do
           end do
        end if
     end do bocos
  end do domains
end subroutine getForces

subroutine setFullMask

  ! Shortcut routine to set all the boundary mask values to 1
  use block
  use inputTimeSpectral
  implicit none

  integer(kind=intType) :: sps, mm, nn

  do sps=1,nTimeIntervalsSpectral
     do nn=1,nDom
        do mm=1,flowDoms(nn, 1, sps)%nBocos
           flowDoms(nn, 1, sps)%bcData(mm)%mask = 1
        end do
     end do
  end do
end subroutine setFullMask

subroutine setNMaskFams(n)

  use costFunctions
  implicit none
  integer(kind=intType), intent(in) ::n

  if (allocated(maskFams)) then 
     deallocate(maskFams)
  end if
  allocate(maskFams(n))
end subroutine setNMaskFams

subroutine setMaskFam(i, fam)

  use costFunctions
  implicit none
  integer(kind=intType), intent(in) :: i
  character(len=maxCGNSNameLen), intent(in) :: fam
  maskFams(i) = fam
end subroutine setMaskFam

subroutine setMask

  ! Check if the family of each boundary condition is in the supplied list.

  use block
  use BCTypes
  use inputTimeSpectral
  use costFunctions
  use cgnsGrid 
  implicit none

  integer(kind=intType) :: sps, mm, nn, kk, cgb, bc
  character(len=maxCGNSNameLen) :: patchName
  character(len=maxCGNSNameLen) :: StrUpCase
  ! Note that we are not setting pointers here for the sake of
  ! efficiency. 
  do sps=1,nTimeIntervalsSpectral
     do nn=1,nDom
        do mm=1,flowDoms(nn, 1, sps)%nBocos
           flowDoms(nn, 1, sps)%bcData(mm)%mask = 0
           
           ! This is an effiicent linear loop...it should be ok as
           ! long as the number of familes isn't excessive.
           bc = flowDoms(nn, 1, sps)%bcType(mm)
           if(bc == EulerWall .or. &
              bc == NSWallAdiabatic .or. &
              bc == NSWallIsothermal) then

              ! Indx of the original CGNS block
              cgb = flowDoms(nn,1,1)%cgnsBlockID

              ! CGNS Family name of this patch
              patchName = cgnsDoms(cgb)%bocoInfo(flowDoms(nn,1,1)%cgnsSubface(mm))%wallBCName

              ! ! Loop over the family names we 
              do kk=1, size(maskFams)
                 if (trim(StrUpCase(maskFams(kk)))  == &
                     trim(StrUpCase(patchName))) then 
                    flowDoms(nn, 1, sps)%bcData(mm)%mask = 1
                 end if
              end do
           end if
        end do
     end do
  end do
end subroutine setMask

FUNCTION StrUpCase ( Input_String ) RESULT ( Output_String )

  ! Borrowed from: http://www.ssec.wisc.edu/~paulv/Fortran90/Utility/
  use constants
  implicit none
  
  ! -- Argument and result
  CHARACTER(len=maxCGNSNameLen), INTENT( IN )     :: Input_String
  CHARACTER(len=maxCGNSNameLen) :: Output_String
  
  ! -- Local variables
  INTEGER(kind=intType) :: i, n

  ! -- Copy input string
  Output_String = Input_String
  
  ! -- Loop over string elements
  DO i = 1, LEN( Output_String )

     ! -- Find location of letter in lower case constant string
     n = INDEX( LOWER_CASE, Output_String( i:i ) )
     
     ! -- If current substring is a lower case letter, make it upper case
     IF ( n /= 0 ) Output_String( i:i ) = UPPER_CASE( n:n )
     
  END DO
  
END FUNCTION StrUpCase
