module io_manager
    use vectors

contains
    subroutine prettyprint_real(mat, n)
        real(rk), dimension(n, n), intent(in) :: mat
        integer(ik), intent(in) :: n
        integer(ik) :: i

        do i = 1, n
            write(6, *) mat(i, :)
        end do
        write(6, *) ""
    end subroutine prettyprint_real

    subroutine prettyprint_vec_mat(vecs, n)
        type(vec), dimension(n, n), intent(in) :: vecs
        integer(ik), intent(in) :: n
        integer(ik) :: i, j

        do i = 1, n
            do j = 1, n
                write(6, *) i, j, ': (', vecs(i, j)%x, vecs(i, j)%y, vecs(i, j)%z, ')'
            end do
        end do
        write(6, *) ""
    end subroutine prettyprint_vec_mat

    subroutine prettyprint_vec_arr(vecs, n)
        type(vec), dimension(n), intent(in) :: vecs
        integer(ik), intent(in) :: n
        integer(ik) :: i

        do i = 1, n
            write(6, *) i, ': (', vecs(i)%x, vecs(i)%y, vecs(i)%z, ')'
        end do
        write(6, *) ""
    end subroutine prettyprint_vec_arr

    subroutine prettyprint_vec_single(a)
        type(vec), intent(in) :: a

        write(6, *) '(', a%x, a%y, a%z, ')'
        write(6, *) ""
    end subroutine prettyprint_vec_single
end module io_manager