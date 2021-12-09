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
    type(vec) :: a, b

    a = vec(1, 2, 0)
    b = vec(2, 1, 0)

    pos(1) = vec(-1., 0., 0.)
    pos(2) = vec(1., 0., 0.)
!    pos(3) = vec(7, 8, 9)
!    pos(4) = vec(10, 11, 12)

    mass(1) = 2.
    mass(2) = 5.
!    mass(3) = 3
!    mass(4) = 4

    mass_mat = get_mass_matrix(mass, n)
    dist_vecs = calc_distance_vecs(pos, n)
    dists = calc_scalar_distances_inv(dist_vecs, n)
    res_forces = calc_forces(dist_vecs, dists, mass_mat, n)

    call prettyprint_real(dists, n)
    call prettyprint_real(mass_mat, n)
    call prettyprint_vec_mat(dist_vecs, n)
    call prettyprint_vec_arr(res_forces, n)
    print *, pot_energy(dists, mass_mat, n)

end program tester