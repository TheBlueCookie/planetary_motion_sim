module plan_sim
    use vectors
    use io_manager

    ! global variables storing simulation results and parameters
    type(vec), dimension(:, :), allocatable :: pos, vel, acc, forces
    real(rk), dimension(:), allocatable :: masses, pot_en, kin_en
    real(rk), dimension(:, :), allocatable :: dists_temp
    real(rk) :: step, d_dur
    integer(ik) :: n_step, n_body, f_step, s_step

    ! initialization marker
    logical :: init_status = .false.

    ! gravitational constant and seconds in a day
    real(rk), parameter :: g = 6.67430e-11 ! https://physics.nist.gov/cgi-bin/cuu/Value?bg
    real(rk), parameter :: sec_per_day = 24 * 3600

contains

    ! reads simulation params from file and allocates arrays accordingly
    subroutine initialize_sim(parfile, nbody, nstep, fstep, sstep, step, ddur)
        character(len = 120), intent(in) :: parfile
        integer(ik), intent(inout) :: nbody, nstep, fstep, sstep
        real(rk), intent(inout) :: ddur, step
        real(rk) :: sdur, rkstep

        call get_arr_dims(parfile, nbody, nstep)

        allocate(pos(nbody, nstep), vel(nbody, nstep), acc(nbody, nstep), forces(nbody, nstep), &
                masses(nbody), pot_en(nstep), kin_en(nstep), dists_temp(nbody, nbody))

        call read_sim_params(parfile, nbody, masses, pos(:, 1), vel(:, 1), ddur, fstep, sstep, nstep)

        sdur = ddur * sec_per_day
        rkstep = nstep
        step = sdur / rkstep

        call info('Initialized simulation.')

        init_status = .true.
    end subroutine initialize_sim

    ! computes a matrix of masses where m_ij = m_i * m_j
    function get_mass_matrix(masses, nbody) result(mass_mat)
        real(rk), dimension(nbody), intent(in) :: masses
        integer(ik), intent(in) :: nbody
        real(rk), dimension(nbody, nbody) :: mass_mat
        integer(ik) :: i, j
        real(rk) :: m

        do i = 1, nbody
            m = masses(i)
            do j = 1, nbody
                if (i == j) then
                    mass_mat(i, j) = 0
                else
                    mass_mat(i, j) = m * masses(j)
                end if
            end do
        end do
    end function get_mass_matrix


    ! full Velocity Verlet algorithm
    subroutine vel_verlet(ipos, ivel, masses, step, nbody, nstep, savestep, fstep)
        type(vec), dimension(nbody), intent(in) :: ipos, ivel
        real(rk), dimension(nbody), intent(in) :: masses
        real(rk), intent(in) :: step
        integer(ik), intent(in) :: nbody, nstep, savestep, fstep
        type(vec), dimension(nbody) :: force_vecs
        type(vec), dimension(nbody, nbody) :: dist_vecs
        real(rk), dimension(nbody, nbody) :: scal_dists_inv, mass_matrix
        integer(ik) :: i, m
        character(len = 100) :: dirpath

        dirpath = prepare_subdir()

        mass_matrix = get_mass_matrix(masses, nbody)
        dist_vecs = calc_distance_vecs(ipos, nbody)
        scal_dists_inv = calc_scalar_distances_inv(dist_vecs, nbody)
        force_vecs = calc_forces(dist_vecs, scal_dists_inv, mass_matrix, nbody)

        kin_en(1) = kin_energy(ivel, masses, nbody)
        pot_en(1) = pot_energy(scal_dists_inv, mass_matrix, nbody)

        do i = 1, nbody
            pos(i, 1) = ipos(i)
            vel(i, 1) = ivel(i)
            acc(i, 1) = force_vecs(i) * (1 / masses(i))
            forces(i, 1) = force_vecs(i)
        end do

        do m = 2, nstep
            call verlet_step(m)
            kin_en(m) = kin_energy(vel(:, m), masses, nbody)
            pot_en(m) = pot_energy(dists_temp, mass_matrix, nbody)

            if (modulo(m, savestep) == 0) then
                call save_sim_step(pos(:, m), vel(:, m), acc(:, m), forces(:, m), kin_en(m), pot_en(m), step * m, &
                        nbody, dirpath)
            end if

            if (modulo(m, fstep) == 0) then
                call print_sim_info(pos(:, m), nbody, nstep, m, step, savestep)
            end if
        end do

    contains

        ! single verlet step
        subroutine verlet_step(j)
            integer(ik), intent(in) :: j
            integer(ik) :: i

            do i = 1, nbody
                pos(i, j) = pos(i, j - 1) + vel(i, j - 1) * step + &
                        forces(i, j - 1) * (1 / masses(i) * 0.5 * step ** 2)
            end do

            dist_vecs = calc_distance_vecs(pos(:, j), nbody)
            scal_dists_inv = calc_scalar_distances_inv(dist_vecs, nbody)
            dists_temp = scal_dists_inv
            forces(:, j) = calc_forces(dist_vecs, scal_dists_inv, mass_matrix, nbody)

            do i = 1, nbody
                acc(i, j) = forces(i, j) * (1 / masses(i))
                vel(i, j) = vel(i, j - 1) + (acc(i, j - 1) + acc(i, j)) * (step * 0.5)
            end do
        end subroutine verlet_step
    end subroutine vel_verlet

    ! calculates the resulting forces acting on all bodies
    function calc_forces(vec_dists, scal_dists, mass_mat, nbody) result(res_forces)
        type(vec), dimension(nbody, nbody), intent(in) :: vec_dists
        real(rk), dimension(nbody, nbody), intent(in) :: scal_dists, mass_mat
        integer(ik), intent(in) :: nbody
        type(vec), dimension(nbody) :: res_forces
        integer(ik) :: i, j
        real(rk) :: temp_fac

        do i = 1, nbody
            res_forces(i) = 0.0_rk
            do j = 1, nbody
                res_forces(i) = res_forces(i) + vec_dists(i, j) * mass_mat(i, j) * scal_dists(i, j) ** 3
            end do
            res_forces(i) = res_forces(i) * g
        end do
    end function calc_forces

    ! calculates total kinetic energy
    function kin_energy(ivel, masses, nbody) result(e_kin)
        type(vec), dimension(nbody), intent(in) :: ivel
        real(rk), dimension(nbody), intent(in) :: masses
        integer(ik), intent(in) :: nbody
        real(rk) :: e_kin
        integer(ik) :: i

        e_kin = 0

        do i = 1, nbody
            e_kin = e_kin + 0.5 * masses(i) * norm(ivel(i)) ** 2
        end do
    end function kin_energy

    ! calculates total potential energy
    function pot_energy(dists, mass_mat, nbody) result(e_pot)
        real(rk), dimension(nbody, nbody), intent(in) :: dists, mass_mat
        integer(ik), intent(in) :: nbody
        real(rk) :: e_pot

        e_pot = 0.5 * g * sum(matmul(mass_mat, dists ** 2))
    end function pot_energy

    ! calculates a n x n matrix of all inverse scalar distances between all bodies from their vectorial distances
    ! r_ij = r_ji = 1/norm(r_i - r_j) and r_ii = 0
    function calc_scalar_distances_inv(vec_dists, nbody) result(distances)
        type(vec), dimension(nbody, nbody), intent(in) :: vec_dists
        integer(ik), intent(in) :: nbody
        real(rk), dimension(nbody, nbody) :: distances
        type(vec) :: tvec
        integer(ik) :: i, j

        do i = 1, nbody
            do j = 1, i
                if (i == j) then
                    distances(i, j) = 0
                else
                    distances(i, j) = 1 / norm(vec_dists(i, j))
                    distances(j, i) = distances(i, j)
                end if
            end do
        end do
    end function calc_scalar_distances_inv

    ! calculates a n x n matrix of all vectorial distances between all bodies from their positions
    ! r_ij = - r_ji = r_i - r_j
    function calc_distance_vecs(ipos, nbody) result(dist_vecs)
        type(vec), dimension(nbody), intent(in) :: ipos
        integer(ik), intent(in) :: nbody
        type(vec), dimension(nbody, nbody) :: dist_vecs
        integer(ik) :: i, j
        type(vec) :: tvec

        do i = 1, nbody
            tvec = ipos(i)
            do j = 1, i
                if (i == j) then
                    dist_vecs(i, j) = 0.0_rk
                else
                    dist_vecs(i, j) = ipos(j) - tvec
                    dist_vecs(j, i) = dist_vecs(i, j) * (-1.0_rk)
                end if
            end do
        end do

    end function calc_distance_vecs
end module plan_sim