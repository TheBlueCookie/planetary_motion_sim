program simulation
    use plan_sim
    use io_manager
    implicit none

    character(len = 120) :: parfile

    call get_command_argument(1, parfile)

    call initialize_sim(parfile, n_body, n_step, f_step, s_step, step, d_dur)

    call vel_verlet(pos, vel, masses, step, n_body, n_step, s_step, f_step)

end program simulation