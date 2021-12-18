module vectors
    implicit none
    integer, parameter :: rk = selected_real_kind(10, 40), ik = selected_int_kind(5)
    real(rk), parameter :: pi = 4.0 * atan(1.0_rk)
    !real(rk), parameter :: alpha_g = 192.85948 * pi / 180.0, delta_g = 27.12825 * pi / 180.0, l_ncp = 122.93192 * pi / 180.0
    ! https://en.wikipedia.org/wiki/Astronomical_coordinate_systems J2000.0

    type :: vec
        real(rk) :: x, y, z
    end type vec

    interface operator(+)
        module procedure vadd
    end interface

    interface operator(-)
        module procedure vminus
    end interface

    interface operator(*)
        module procedure scalarprod
        module procedure scalar_mult_a
        module procedure scalar_mult_b
    end interface

    interface assignment(=)
        module procedure vec_to_vec
        module procedure scalar_to_vec
        module procedure arr_to_vec
    end interface

contains
    function vadd(a, b) result(c)
        type(vec), intent(in) :: a, b
        type(vec) :: c

        c%x = a%x + b%x
        c%y = a%y + b%y
        c%z = a%z + b%z
    end function vadd

    function vminus(a, b) result(c)
        type(vec), intent(in) :: a, b
        type(vec) :: c

        c%x = a%x - b%x
        c%y = a%y - b%y
        c%z = a%z - b%z
    end function vminus

    function scalarprod(a, b) result(c)
        type(vec), intent(in) :: a, b
        real(rk) :: c

        c = a%x * b%x + a%y * b%y + a%z * b%z
    end function scalarprod

    function scalar_mult_a(a, b) result(c)
        type(vec), intent(in) :: a
        real(rk), intent(in) :: b
        type(vec) :: c

        c%x = a%x * b
        c%y = a%y * b
        c%z = a%z * b
    end function scalar_mult_a

    function scalar_mult_b(a, b) result(c)
        real(rk), intent(in) :: a
        type(vec), intent(in) :: b
        type(vec) :: c

        c%x = b%x * a
        c%y = b%y * a
        c%z = b%z * a
    end function scalar_mult_b

    function norm(a) result(b)
        type(vec), intent(in) :: a
        real(rk) :: b

        b = sqrt(a * a)
    end function norm

    function cross_product(a, b) result(c)
        type(vec), intent(in) :: a, b
        type(vec) :: c

        c%x = a%y * b%z - a%z * b%y
        c%y = a%z * b%x - a%x * b%z
        c%z = a%x * b%y - a%y * b%x
    end function cross_product

    subroutine vec_to_vec(a, b)
        type(vec), intent(inout) :: a
        type(vec), intent(in) :: b

        a%x = b%x
        a%y = b%y
        a%z = b%z
    end subroutine vec_to_vec

    subroutine scalar_to_vec(a, b)
        type(vec), intent(inout) :: a
        real(rk), intent(in) :: b

        a%x = b
        a%y = b
        a%z = b
    end subroutine scalar_to_vec

    subroutine arr_to_vec(a, b)
        type(vec), intent(inout) :: a
        real(rk), dimension(3), intent(in) :: b

        a%x = b(1)
        a%y = b(2)
        a%z = b(3)
    end subroutine arr_to_vec
end module vectors