module vectors
    implicit none
    integer, parameter :: rk = selected_real_kind(10, 40), ik = selected_int_kind(5)
    real(rk), parameter :: pi = 4.0 * atan(1.0_rk)

    ! type definition
    type :: vec
        real(rk) :: x, y, z
    end type vec

    ! vectorial addition
    interface operator(+)
        module procedure vadd
    end interface

    ! vectorial subtraction
    interface operator(-)
        module procedure vminus
    end interface

    ! scalar product of two vectors, vector * scalar and scalar * vector
    interface operator(*)
        module procedure scalarprod
        module procedure scalar_mult_a
        module procedure scalar_mult_b
    end interface

    ! vector = vector, vector = scalar -> setting all compononents to one value, vector = real, dim(3) -> copying array entries to vector
    interface assignment(=)
        module procedure vec_to_vec
        module procedure scalar_to_vec
        module procedure arr_to_vec
    end interface

contains
    ! vectorial addition
    function vadd(a, b) result(c)
        type(vec), intent(in) :: a, b
        type(vec) :: c

        c%x = a%x + b%x
        c%y = a%y + b%y
        c%z = a%z + b%z
    end function vadd

    ! vectorial subtraction
    function vminus(a, b) result(c)
        type(vec), intent(in) :: a, b
        type(vec) :: c

        c%x = a%x - b%x
        c%y = a%y - b%y
        c%z = a%z - b%z
    end function vminus

    ! vectorial scalar product a*b = a_i*b_i
    function scalarprod(a, b) result(c)
        type(vec), intent(in) :: a, b
        real(rk) :: c

        c = a%x * b%x + a%y * b%y + a%z * b%z
    end function scalarprod

    ! vector * scalar
    function scalar_mult_a(a, b) result(c)
        type(vec), intent(in) :: a
        real(rk), intent(in) :: b
        type(vec) :: c

        c%x = a%x * b
        c%y = a%y * b
        c%z = a%z * b
    end function scalar_mult_a

    ! scalar * vector
    function scalar_mult_b(a, b) result(c)
        real(rk), intent(in) :: a
        type(vec), intent(in) :: b
        type(vec) :: c

        c%x = b%x * a
        c%y = b%y * a
        c%z = b%z * a
    end function scalar_mult_b

    ! standard 2-norm of euclidean space
    function norm(a) result(b)
        type(vec), intent(in) :: a
        real(rk) :: b

        b = sqrt(a * a)
    end function norm

    ! cross product of two vectors (a x b)_i = e_ijk a_j * b_k where e is the Levi-Civita symbol
    function cross_product(a, b) result(c)
        type(vec), intent(in) :: a, b
        type(vec) :: c

        c%x = a%y * b%z - a%z * b%y
        c%y = a%z * b%x - a%x * b%z
        c%z = a%x * b%y - a%y * b%x
    end function cross_product

    ! assigns vector = vector
    subroutine vec_to_vec(a, b)
        type(vec), intent(inout) :: a
        type(vec), intent(in) :: b

        a%x = b%x
        a%y = b%y
        a%z = b%z
    end subroutine vec_to_vec

    ! assigns vector a = scalar b -> a_x = a_y = a_z = b
    subroutine scalar_to_vec(a, b)
        type(vec), intent(inout) :: a
        real(rk), intent(in) :: b

        a%x = b
        a%y = b
        a%z = b
    end subroutine scalar_to_vec

    ! vector a = b(3) -> a_x = b(1), a_y = b(2), a_z = b(3)
    subroutine arr_to_vec(a, b)
        type(vec), intent(inout) :: a
        real(rk), dimension(3), intent(in) :: b

        a%x = b(1)
        a%y = b(2)
        a%z = b(3)
    end subroutine arr_to_vec
end module vectors