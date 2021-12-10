program tester
    use plan_sim
    use io_manager
    implicit none

    integer(ik), parameter :: n = 2, m = 100
    real(rk), parameter :: step = 1

    type(vec), dimension(n) :: pos, res_forces, vel
    type(vec), dimension(n, n) :: dist_vecs, force_vecs
    integer(ik) :: i
    real(rk), dimension(n, n) :: dists, mass_mat
    real(rk), dimension(n) :: mass
    type(vec) :: a, b
    type(verlet_solution(n, m)) :: solution

    mass(1) = 1.989e30 ! sun
    mass(2) = 5.972e24 ! earth

    pos(1) = vec(0., 0., 0.)
    pos(2) = vec(150e9, 0., 0.)

    vel(1) = vec(0., 0., 0.)
    vel(2) = vec(0., 30e3, 0.)

    solution = vel_verlet(pos, vel, mass, step, m, n)

end program tester