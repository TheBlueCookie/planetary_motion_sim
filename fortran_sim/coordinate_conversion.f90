program sim_preperation
    use plan_sim
    implicit none

    integer(ik), parameter :: body_count = 9 - 1
    real(rk), dimension(body_count) :: periapsis, apoapsis, eccentricity, a, vel_at_periapsis
    real(rk) :: gm
    integer(ik) :: i

    ! https://nssdc.gsfc.nasa.gov/planetary/factsheet/
    periapsis = [46.0, 107.5, 147.1, 206.6, 740.5, 1352.6, 2741.3, 4444.5]

    ! https://www.princeton.edu/~willman/planetary_systems/Sol/
    a = [0.3870993, 0.723336, 1.000003, 1.52371, 5.2029, 9.537, 19.189, 30.0699]

    periapsis = periapsis * 1e9
    apoapsis = apoapsis * 1e9
    a = a * 1.495978707e11

    gm = 1.327124400e20 ! https://ssd.jpl.nasa.gov/astro_par.html

    vel_at_periapsis = sqrt(gm * (2.0 / periapsis - 1 / a))

    print *, a
    print *, gm
    print *, periapsis
    print *, vel_at_periapsis

end program sim_preperation