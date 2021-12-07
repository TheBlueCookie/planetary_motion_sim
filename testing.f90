program tester
    use plan_sim
    implicit none

    type(vec), dimension(4) :: pos
    integer(ik) :: i
    real(rk), dimension(4, 4) :: dists

    pos(1) = vec(1, 2, 3)
    pos(2) = vec(4, 5, 6)
    pos(3) = vec(7, 8, 9)
    pos(4) = vec(10, 11, 12)

    do i = 1, 10000000
        dists = calc_distances(pos, 4)
    end do


end program tester