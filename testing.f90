program tester
    use plan_sim
    use io_manager
    implicit none

    integer(ik), parameter :: n = 2

    type(vec), dimension(n) :: pos, res_forces
    type(vec), dimension(n, n) :: dist_vecs, force_vecs
    integer(ik) :: i
    real(rk), dimension(n, n) :: dists, mass_mat
    real(rk), dimension(n) :: mass

    pos(1) = vec(-1, 0, 0)
    pos(2) = vec(1, 0, 0)
!    pos(3) = vec(7, 8, 9)
!    pos(4) = vec(10, 11, 12)

    mass(1) = 3
    mass(2) = 2
!    mass(3) = 3
!    mass(4) = 4

    dists = calc_distances(pos, n)
    mass_mat = get_mass_matrix(mass, n)
    res_forces = calc_forces(dists, pos, mass_mat, n)
    dist_vecs = calc_distance_vecs(pos, n)

    call prettyprint_real(dists, n)
    call prettyprint_real(mass_mat, n)
    print *, ""
    print *, dist_vecs
    print *, ""
    print *, res_forces
    print *, pot_energy(dists, mass_mat, n)

end program tester