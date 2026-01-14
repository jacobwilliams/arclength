    program main

    !! Example program to compute the arclength of a satellite orbit
    !! using numerical quadrature.

    use gooding_module
    use quadrature_module
    use iso_fortran_env, only: wp => real64

    implicit none

    real(wp),parameter :: mu_earth = 398600.4418_wp ! km^3/s^2

    type(integration_class_1d) :: integrator
    real(wp),dimension(6) :: rv0 !! initial state vector
    real(wp) :: arclength, err, t0, tf
    integer :: ierr, n_funcs

    t0 = 0.0_wp ! initial time
    tf = 3600.0_wp  ! final time (1 hour)
    rv0 = [8000.0_wp,1000.0_wp,0.0_wp,0.0_wp,7.546_wp,1.0_wp]  ! some orbit
    n_funcs = 0 ! number of function evaluations

    call integrator%initialize(func, t0, tf, tolx = 1.0e-13_wp, methodx = quad_gauss_8%n_points)
    call integrator%integrate(arclength, ierr, err)

    print*, 'Arclength            = ', arclength
    print*, 'Estimated error      = ', err
    print*, 'Error flag           = ', ierr
    print*, 'Number of func evals = ', n_funcs

    contains

        function func(me,x) result(f)
            !! 1d user function f(x)
            implicit none
            class(integration_class_1d),intent(inout) :: me
            real(wp), intent(in)                      :: x !! this is dt
            real(wp)                                  :: f
            real(wp),dimension(6) :: rvf !! final state vector
            call propagate(mu_earth, rv0, x, rvf)
            f = norm2(rvf(4:6))  ! integrate vmag get arc length
            n_funcs = n_funcs + 1
        end function func

    end program main
