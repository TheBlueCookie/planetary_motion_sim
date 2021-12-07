module plan_sim
    use vectors

    type :: verlet_solution(nbody, nstep)
    integer(ik), len :: nbody, nstep
        type(vec) :: pos(nbody, nstep), vel(nbody, nstep), acc(nbody, nstep)
    end type verlet_solution

contains

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

!    function forces(ipos, masses, nbody) result(res_forces)
!        type(grouped_vec(nbody)), intent(in) :: ipos
!        real(rk), dimension(nbody), intent(in) :: masses
!        integer(ik), intent(in) :: nbody
!
!    end function forces

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

    function calc_distances(ipos, nbody) result(distances)
        type(vec), dimension(nbody), intent(in) :: ipos
        integer(ik), intent(in) :: nbody
        real(rk), dimension(nbody, nbody) :: distances
        type(vec) :: tvec
        integer(ik) :: i, j

        do i = 1, nbody
            tvec = ipos(i)
            do j = 1, i
                distances(i, j) = norm(tvec - ipos(j))
                distances(j, i) = distances(i, j)
            end do
        end do
    end function calc_distances
end module plan_sim