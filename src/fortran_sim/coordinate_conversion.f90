program sim_preperation
    use vectors
    implicit none

    integer(ik), parameter :: body_count = 10 - 1
    real(rk), dimension(body_count) :: periapsis, a, vel_at_periapsis
    real(rk) :: gm, gme
    integer(ik) :: i

    ! https://nssdc.gsfc.nasa.gov/planetary/factsheet/
    periapsis = [46.0, 107.5, 147.1, 206.6, 740.5, 1352.6, 2741.3, 4444.5, 0.363]

    ! semi major axis
    ! https://www.princeton.edu/~willman/planetary_systems/Sol/ (planets) and https://en.wikipedia.org/wiki/Orbit_of_the_Moon (moon)
    a = [0.3870993, 0.723336, 1.000003, 1.52371, 5.2029, 9.537, 19.189, 30.0699, 0.0025718815]

    periapsis = periapsis * 1e9
    a = a * 1.495978707e11

    ! general graviational constant of the sun
    gm = 1.327124400e20 ! https://ssd.jpl.nasa.gov/astro_par.html
    ! general graviational constant of the earth
    gme = 398600.435507e9

    ! periapsis velocity
    vel_at_periapsis(:8) = sqrt(gm * (2.0 / periapsis(:8) - 1 / a(:8)))
    vel_at_periapsis(9) = sqrt(gme * (2.0 / periapsis(9) - 1 / a(9)))

    print *, a
    print *, gm, gme
    print *, periapsis
    print *, vel_at_periapsis

end program sim_preperation