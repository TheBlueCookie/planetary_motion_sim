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
    end subroutine prettyprint_real

end module io_manager