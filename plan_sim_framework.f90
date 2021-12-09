module plan_sim
    use vectors
    use io_manager

    type :: verlet_solution(nbody, nstep)
    integer(ik), len :: nbody, nstep
        type(vec) :: pos(nbody, nstep), vel(nbody, nstep), acc(nbody, nstep)
    end type verlet_solution

    real(rk), parameter :: g = 1 !6.67430e-11 ! https://physics.nist.gov/cgi-bin/cuu/Value?bg

contains

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
                    mass_mat(i, j) = 0 !m * masses(j)
                else
                    mass_mat(i, j) = m * masses(j)
                end if
            end do
        end do
    end function get_mass_matrix


    function vel_verlet(ipos, ivel, iacc, step, nstep, nbody) result(sol)
        type(vec), dimension(nbody), intent(in) :: ipos, ivel, iacc
        real(rk), intent(in) :: step
        integer(ik), intent(in) :: nstep, nbody
        type(verlet_solution(nbody, nstep)) :: sol

    end function vel_verlet

    !    function verlet_step(ipos, ivel, iacc, step, nbody) result(sol)
    !        type(grouped_vec(nbody)), intent(in) :: ipos, ivel, iacc
    !        real(rk), intent(in) :: step
    !        integer(ik), intent(in) :: nbody
    !        type(verlet_solution(nbody, 1)) :: sol
    !
    !        sol%pos = ipos + ivel * step + 0.5 * iacc * step **2
    !        sol%acc =
    !        sol$vel = ivel + 0.5 * ( )
    !
    !    end function verlet_step

    function calc_forces(dists, ipos, mass_mat, nbody) result(res_forces)
        real(rk), dimension(nbody, nbody), intent(in) :: dists, mass_mat
        type(vec), dimension(nbody) :: ipos
        integer(ik), intent(in) :: nbody
        type(vec), dimension(nbody) :: res_forces
        real(rk), dimension(nbody, nbody) :: t_grav_mat
        integer(ik) :: i

        t_grav_mat = g * matmul(mass_mat, dists ** 3)

        do i = 1, nbody
            res_forces(i)%x = sum(t_grav_mat(i, :)) * ipos(i)%x
            res_forces(i)%y = sum(t_grav_mat(i, :)) * ipos(i)%y
            res_forces(i)%z = sum(t_grav_mat(i, :)) * ipos(i)%z
        end do

        write(6, *) res_forces

    end function calc_forces

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

    function pot_energy(dists, mass_mat, nbody) result(e_pot)
        real(rk), dimension(nbody, nbody), intent(in) :: dists, mass_mat
        integer(ik), intent(in) :: nbody
        real(rk) :: e_pot, m
        integer(ik) :: i, j

        e_pot = 0.5 * g * sum(matmul(mass_mat, dists))
    end function pot_energy

    function calc_distances(ipos, nbody) result(distances)
        type(vec), dimension(nbody), intent(in) :: ipos
        integer(ik), intent(in) :: nbody
        real(rk), dimension(nbody, nbody) :: distances
        type(vec) :: tvec
        integer(ik) :: i, j

        do i = 1, nbody
            tvec = ipos(i)
            do j = 1, i
                if (j == i) then
                    distances(i, j) = 0
                else
                    distances(i, j) = 1 / norm(tvec - ipos(j))
                    distances(j, i) = distances(i, j)
                end if
            end do
        end do
    end function calc_distances

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
                end if
            end do
        end do

    end function calc_distance_vecs
end module plan_sim