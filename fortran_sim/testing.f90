program tester
    use plan_sim
    use io_manager
    implicit none

    integer(ik), parameter :: n = 2, m = 100000
    real(rk), parameter :: step = 100

    type(vec), dimension(n) :: sim_pos, sim_forces, sim_vel
    type(vec), dimension(n, n) :: dist_vecs, force_vecs
    integer(ik) :: i
    real(rk), dimension(n, n) :: dists, mass_mat
    real(rk), dimension(n) :: mass

    call initialize(n, m)

    mass(1) = 1.989e30 ! sun
    mass(2) = 5.972e24 ! earth

    sim_pos(1) = vec(0., 0., 0.)
    sim_pos(2) = vec(150e9, 0., 0.)

    sim_vel(1) = vec(0., 0., 0.)
    sim_vel(2) = vec(0., 30e3, 0.)

    call vel_verlet(sim_pos, sim_vel, mass, step, m, n)

    call prettyprint_vec_arr(pos(2, :), m)
end program tester