module plan_sim
    implicit none
    integer, parameter :: rk = selected_real_kind(10, 40), ik = selected_int_kind(5)

    type :: vec
        real(rk) :: x, y, z
    end type vec

    type :: grouped_vec(nbody)
        integer(ik), len :: nbody
        type(vec) :: vecs(nbody)
    end type grouped_vec

    type :: verlet_solution(nbody, nstep)
        integer(ik), len :: nbody, nstep
        type(vec) :: pos(nbody, nstep), vel(nbody, nstep), acc(nbody, nstep)
    end type verlet_solution

contains

    function vel_verlet(ipos, ivel, iacc, step, nstep, nbody) result(sol)
        type(grouped_vec(nbody)), intent(in) :: ipos, ivel, iacc
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

    function forces(ipos, masses, nbody) result(forces)
        type(grouped_vec(nbody)), intent(in) :: ipos
        real(rk), dimension(nbody), intent(in) :: masses
        integer(ik), intent(in) :: nbody



    end function forces


end module plan_sim