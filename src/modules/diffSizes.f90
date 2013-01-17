  module diffSizes
    use precision
    implicit none
    save
    
    integer(kind=intType), parameter :: ISIZE3ofviscsubface = 3
    !integer(kind=intType) :: ISIZE1ofviscsubface
    integer(kind=intType) :: ISIZE1OFDrfbcdata
    integer(kind=intType) :: ISIZE1OFDrfviscsubface
    integer(kind=intType) :: ISIZE1OFDrfflowdoms
    integer(kind=intType) :: ISIZE2OFDrfflowdoms
    integer(kind=intType) :: ISIZE3OFDrfflowdoms

  end module diffSizes
