program tester
    use plan_sim
    use io_manager
    implicit none

    character(len = 120) :: parfile = 'sim_params.txt'

    call initialize_sim(parfile, n_body, n_step, f_step, s_step, step, d_dur)
    !
    !    mass(1) = 1.989e30 ! sun
    !    mass(2) = 5.972e24 ! earth
    !    mass(3) = 2.972e24
    !
    !    sim_pos(1) = vec(0., 0., 0.)
    !    sim_pos(2) = vec(150e9, 0., 0.)
    !    sim_pos(3) = vec(90e9, 20e9, 0.)
    !
    !    sim_vel(1) = vec(0., 0., 0.)
    !    sim_vel(2) = vec(0., 20e3, 0.)
    !    sim_vel(3) = vec(5e3, 0., 0e3)

    call vel_verlet(pos, vel, masses, step, n_body, n_step, s_step, f_step)

    ! call prettyprint_vec_arr(pos(2, :), m)

    ! call export_all_sim_results(pos, vel, acc, forces, n_body, n_step, s_step)
end program tester