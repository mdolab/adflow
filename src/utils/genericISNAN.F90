module genericISNAN
    implicit none
    interface myIsNAN
     module procedure myIsNAN_r
     module procedure myIsNAN_c
    end interface

    contains

    logical function myIsNAN_r(val)
    !
    !       myIsNAN_r determines whether or not the given real value is a NAN or INF and
    !       returns the according logical.
    !
    use precision
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
    implicit none
    !
    !      Function arguments.
    !
    real(kind=alwaysRealType), intent(in) :: val

    ! Check if NAN or INF
    myIsNAN_r = ieee_is_nan(val) .or. .not. ieee_is_finite(val)

    end function myIsNAN_r

    logical function myIsNAN_c(val)
    !
    !       myIsNAN_c determines whether or not the given complex value contains NAN of INF and
    !       returns the according logical.
    !
    use precision
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
    implicit none
    !
    !      Function arguments.
    !
    complex(kind=realType), intent(in) :: val

    ! Check if either real or imag part is NAN or INF
    myIsNAN_c = ieee_is_nan(real(val)) .or. .not. ieee_is_finite(real(val))
    myIsNAN_c = myIsNAN_c .or. ieee_is_nan(aimag(val)) .or. .not. ieee_is_finite(aimag(val))

    end function myIsNAN_c

end module genericISNAN