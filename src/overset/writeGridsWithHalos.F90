subroutine writeGridsWithHalos()

  !
  !      ******************************************************************
  !      *                                                                *
  !      * Special utiliy function to write all unsplit meshes to an      *
  !      * unformatted plot3d file including the halos. We temporarily    *
  !      * use this to generate the overset connectivity externally       *
  !      *                                                                *
  !      ******************************************************************
  !
  use blockPointers
  implicit none
  
  ! Working variables
  integer(kind=intType) :: nn, i, j, k, logg
  character*40 :: tmp, integer_string

  ! Simply loop over blocks and dump each to its own file.
  do nn=1,nDom
     
     ! Get the filename
     write(integer_string,*) nn
     integer_string = adjustl(integer_string)
     tmp = 'g.'//trim(integer_string)
     
     ! Open up the file
     open(unit=1, file=tmp, form='unformatted', status='unknown')
     
     ! Set pointers for this block...only first level 
     call setPointers(nn, 1, 1)

     logg = 1
     write(logg) ie+1, je+1, ke+1
     write(logg) ((( x(i,j,k,1), i=0,ie), j=0,je), k=0,ke)
     write(logg) ((( x(i,j,k,2), i=0,ie), j=0,je), k=0,ke)
     write(logg) ((( x(i,j,k,3), i=0,ie), j=0,je), k=0,ke)
     close(logg)

  end do
end subroutine writeGridsWithHalos

