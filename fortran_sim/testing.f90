program tester
    use plan_sim
    use io_manager
    implicit none

    integer(ik), parameter :: n = 3, m = 100000, k = 10, fstep = 100000
    real(rk), parameter :: step = 300

    type(vec), dimension(n) :: sim_pos, sim_forces, sim_vel
    type(vec), dimension(n, n) :: dist_vecs, force_vecs
    integer(ik) :: i
    real(rk), dimension(n, n) :: dists, mass_mat
    real(rk), dimension(n) :: mass

    call initialize_sim(n, m, k)

    mass(1) = 1.989e30 ! sun
    mass(2) = 5.972e24 ! earth
    mass(3) = 2.972e24

    sim_pos(1) = vec(0., 0., 0.)
    sim_pos(2) = vec(150e9, 0., 0.)
    sim_pos(3) = vec(90e9, 20e9, 0.)

    sim_vel(1) = vec(0., 0., 0.)
    sim_vel(2) = vec(0., 20e3, 0.)
    sim_vel(3) = vec(5e3, 0., 0e3)

    call vel_verlet(sim_pos, sim_vel, mass, step, n, m, k, fstep)

    ! call prettyprint_vec_arr(pos(2, :), m)

    call export_all_sim_results(pos, vel, acc, forces, n, m, k)
end program tester