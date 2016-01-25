
subroutine initializeOBlock(oBlock, nn)

  ! This routine allocates the data for the supplied oBlock using the
  !  data currently in blockPointers
  use constants
  use overset
  use blockPointers
  use adtAPI
  use BCTypes
  use cgnsGrid
  use communication
  implicit none 

  ! Input Params
  type(oversetBlock), intent(inout) :: oBlock
  integer(kind=intType) :: nn, kk

  ! Working paramters
  integer(kind=intType) :: i, j, k, mm, nADT, nHexa, planeOffset
  integer(kind=intType) :: iStart, iEnd, jStart, jEnd, kStart, kEnd

  ! Set all the sizes for this block.
  oBlock%il = il
  oBlock%jl = jl
  oBlock%kl = kl

  oBlock%proc = myID
  oBlock%block = nn

  allocate( &
       oBlock%qualDonor(1, ie*je*ke), &
       oBlock%globalCell(0:ib, 0:jb, 0:kb), &
       oBlock%nearWall(1:ie, 1:je, 1:ke), &
       oBlock%invalidDonor(1:ie, 1:je, 1:ke))

  oBlock%nearWall = 0
  oBlock%invalidDonor = 0

  kk = 10
  do mm=1,nBocos
     select case (BCFaceID(mm))
     case (iMin)
        iStart=1; iEnd=2+kk;
        jStart=BCData(mm)%icBeg; jEnd=BCData(mm)%icEnd
        kStart=BCData(mm)%jcBeg; kEnd=BCData(mm)%jcEnd
     case (iMax)
        iStart=ie-kk; iEnd=ie;
        jStart=BCData(mm)%icBeg; jEnd=BCData(mm)%icEnd
        kStart=BCData(mm)%jcBeg; kEnd=BCData(mm)%jcEnd
     case (jMin)
        iStart=BCData(mm)%icBeg; iEnd=BCData(mm)%icEnd
        jStart=1; jEnd=1+kk;
        kStart=BCData(mm)%jcBeg; kEnd=BCData(mm)%jcEnd
     case (jMax)
        iStart=BCData(mm)%icBeg; iEnd=BCData(mm)%icEnd
        jStart=je-kk; jEnd=je;
        kStart=BCData(mm)%jcBeg; kEnd=BCData(mm)%jcEnd
     case (kMin)
        iStart=BCData(mm)%icBeg; iEnd=BCData(mm)%icEnd
        jStart=BCData(mm)%jcBeg; jEnd=BCData(mm)%jcEnd
        kStart=1; kEnd=1+kk;
     case (kMax)
        iStart=BCData(mm)%icBeg; iEnd=BCData(mm)%icEnd
        jStart=BCData(mm)%jcBeg; jEnd=BCData(mm)%jcEnd
        kStart=ke-kk; kEnd=ke;
     end select

     if (BCType(mm) == NSWallAdiabatic .or. &
          BCType(mm) == NSWallIsoThermal .or. &
          BCType(mm) == EulerWall) then 

        ! Clip Bounds:
        iStart = max(1, iStart)
        iEnd   = min(ie, iEnd)
        jStart = max(1, jStart)
        jEnd   = min(je, jEnd)
        kStart = max(1, kStart)
        kEnd   = min(ke, kEnd)

        do k=kStart, kEnd
           do j=jStart, jEnd
              do i=iStart, iEnd
                 oBlock%nearWall(i, j, k) = 1
              end do
           end do
        end do
     end if
  end do ! BocoLoop

  do mm=1,nBocos
     ! Just record the ranges necessary and we'll add in a generic
     ! loop. Why is it the first three? Well, the first level of halos
     ! off of an overset outer bound is completely
     ! meaningless. Essentially we ignore those. So the outer two
     ! layers of cells are indices 2 and 3. Therefore the first 3 on
     ! either side need to be flagged as invalid.

     select case (BCFaceID(mm))
     case (iMin)
        iStart=1; iEnd=3;
        jStart=BCData(mm)%inBeg+1; jEnd=BCData(mm)%inEnd
        kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
     case (iMax)
        iStart=nx; iEnd=ie;
        jStart=BCData(mm)%inBeg+1; jEnd=BCData(mm)%inEnd
        kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
     case (jMin)
        iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
        jStart=1; jEnd=3;
        kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
     case (jMax)
        iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
        jStart=ny; jEnd=je;
        kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
     case (kMin)
        iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
        jStart=BCData(mm)%jnBeg+1; jEnd=BCData(mm)%jnEnd
        kStart=1; kEnd=3;
     case (kMax)
        iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
        jStart=BCData(mm)%jnBeg+1; jEnd=BCData(mm)%jnEnd
        kStart=nz; kEnd=ke;
     end select


     ! select case (BCFaceID(mm))
     ! case (iMin)
     !    iStart=1; iEnd=2;
     !    jStart=BCData(mm)%inBeg+1; jEnd=BCData(mm)%inEnd
     !    kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
     ! case (iMax)
     !    iStart=nx; iEnd=il;
     !    jStart=BCData(mm)%inBeg+1; jEnd=BCData(mm)%inEnd
     !    kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
     ! case (jMin)
     !    iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
     !    jStart=1; jEnd=2;
     !    kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
     ! case (jMax)
     !    iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
     !    jStart=ny; jEnd=jl;
     !    kStart=BCData(mm)%jnBeg+1; kEnd=BCData(mm)%jnEnd
     ! case (kMin)
     !    iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
     !    jStart=BCData(mm)%jnBeg+1; jEnd=BCData(mm)%jnEnd
     !    kStart=1; kEnd=2;
     ! case (kMax)
     !    iStart=BCData(mm)%inBeg+1; iEnd=BCData(mm)%inEnd
     !    jStart=BCData(mm)%jnBeg+1; jEnd=BCData(mm)%jnEnd
     !    kStart=nz; kEnd=kl;
     ! end select

     if (BCType(mm) == OversetOuterBound) then
        do k=kStart, kEnd
           do j=jStart, jEnd
              do i=iStart, iEnd
                 ! Compute the index
                 oBlock%invalidDonor(i, j, k) = 1
              end do
           end do
        end do
     end if
  end do

  ! Copy Volume to qualDonor and do minVol while we're at it
  oBlock%minVol = Large
  mm = 0
  do k=1,ke
     do j=1,je
        do i=1,ie
           mm = mm + 1
           oBlock%qualDonor(1, mm) = vol(i, j, k)
           oBlock%minVol = min(oBlock%minVol, vol(i, j, k))
        end do
     end do
  end do

  do k=0,kb
     do j=0,jb
        do i=0,ib
           oBlock%globalCell(i, j, k) = globalCell(i, j, k)
        end do
     end do
  end do

  ! Now setup the data for the ADT
  nHexa = il * jl * kl
  nADT = ie * je * ke

  allocate(oBlock%xADT(3, nADT), oBlock%hexaConn(8, nHexa))

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
        end do
     end do
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

  ! Call the custom build routine -- Serial only, only Hexa volumes,
  ! we supply our own ADT Type

  call buildSerialHex(nHexa, nADT, oBlock%xADT, oBlock%hexaConn, oBlock%ADT)

  ! Flag this block as being allocated
  oBlock%allocated = .True.

end subroutine initializeOBlock
