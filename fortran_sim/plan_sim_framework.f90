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


    function vel_verlet(ipos, ivel, masses, step, nstep, nbody) result(sol)
        type(vec), dimension(nbody), intent(in) :: ipos, ivel
        real(rk), dimension(nbody), intent(in) :: masses
        real(rk), intent(in) :: step
        integer(ik), intent(in) :: nstep, nbody
        type(verlet_solution(nbody, nstep)) :: sol
        type(vec), dimension(nbody) :: force_vecs
        type(vec), dimension(nbody, nbody) :: dist_vecs
        real(rk), dimension(nbody, nbody) :: scal_dists_inv, mass_matrix
        integer(ik) :: i, j

        mass_matrix = get_mass_matrix(masses, nbody)
        dist_vecs = calc_distance_vecs(ipos, nbody)
        scal_dists_inv = calc_scalar_distances_inv(dist_vecs, nbody)
        force_vecs = calc_forces(dist_vecs, scal_dists_inv, mass_matrix, nbody)

        do i = 1, nbody
            sol%pos(i, 1) = ipos(i)
            sol%vel(i, 1) = ivel(i)
            sol%acc(i, 1) = norm(force_vecs(i)) / masses(i)
        end do

        do i = 2, nstep
            do j = 1, nbody
                sol%pos(j, i) = sol%pos(i-1, j) + sol%vel(i-1, j) + 0.5 * sol%acc(i-1, j) * step ** 2
            end do
        end do

    end function vel_verlet

!        function verlet_step(ipos, ivel, force_vecs, step, nbody) result(sol)
!            type(vec), dimension(nbody), intent(in) :: ipos, ivel, force_vecs
!            real(rk), intent(in) :: step
!            integer(ik), intent(in) :: nbody
!            type(verlet_solution(nbody, 1)) :: sol
!            integer(ik) :: i
!
!            do i = 1, nbody
!                sol%pos(i, 1) = ipos(i) + ivel(i) * norm(force_vecs(i)) * step + 0.5 * iacc * step **2
!                sol%acc =
!                sol$vel = ivel + 0.5 * ()
!            end do
!
!        end function verlet_step

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
        real(rk) :: e_pot

        e_pot = 0.5 * g * sum(matmul(mass_mat, dists ** 2))
    end function pot_energy

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