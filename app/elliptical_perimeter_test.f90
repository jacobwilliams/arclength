program elliptical_perimeter_test

    use carlson_elliptic_module, wp => carlson_elliptic_module_wp
    use fortran_astrodynamics_toolkit
    use quadpack

    integer  :: ier !! status code
    real(wp) :: p !! semiparameter [km]
    real(wp) :: period !! orbital period [sec]
    real(wp) :: dt !! time from initial state to periapsis [sec]
    real(wp),dimension(3) :: r, v !! position and velocity vectors [km] and [km/s]
    real(wp),dimension(6) :: rv0 !! initial state vector [km, km/s]
    integer :: n_evals !! number of function evaluations
    real(wp) :: e !! eccentricity of the ellipse
    real(wp) :: result !! quadrature result
    real(wp) :: error

    real(wp),parameter :: a      = 8000.0_wp
    real(wp),parameter :: b      = 7000.0_wp
    real(wp),parameter :: mu     = 398600.4418_wp !! km^3/s^2
    real(wp),parameter :: epsrel = 10.0_wp**(int(log10(epsilon(1.0_wp) * 100))) !! tolerance (1e-13 for real64)

    write(*,*) ''
    write(*,*) '-----------------------------'
    write(*,*) 'epsrel=',epsrel
    write(*,*) '-----------------------------'

    write(*,*) ''
    write(*,*) 'elliptical_perimeter_test'
    write(*,*) 'perimeter:        ', perimeter_of_ellipse(a,b)
    write(*,*) 'perimeter:        ', perimeter_of_ellipse(a,a)

    write(*,*) 'test with circle: ', perimeter_of_ellipse(b,b), ', error: ', perimeter_of_ellipse(b,b) - (2.0_wp * pi * b)
    write(*,*) ''

    ! compare to quadrature integration of the ellipse
    write(*,*) 'compare circle to quadrature integration'
    call compare(a,b)
    call compare(a,a)
    call compare(b,b)

    ! this one can't get requested accuracy from the quadrature routines:
    call compare(a*3,b)

    contains

    subroutine compare(a,b)
        ! compare to quadrature integration of the ellipse
        real(wp),intent(in) :: a, b
        write(*,*) ''
        write(*,*) 'a,b,e = ', a, b, e
        e = sqrt(1.0_wp - (b*b)/(a*a))
        period = orbit_period(mu,a)
        p = a * (1.0_wp - e**2)
        call orbital_elements_to_rv(mu, p, e, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, r, v)
        rv0 = [r, v]
        ! integrate one orbit period to get the perimeter of the ellipse:
        call dquad(f, 0.0_wp, period, result, epsrel, n_evals, ier)  ! qquad for real128
        !call dgauss8(f, 0.0_wp, period, epsrel, result, ier, error) ! qgauss8 for real128

                  write(*,*) 'elliptic function result:       ', perimeter_of_ellipse(a,b)
        if (a==b) write(*,*) 'analytical result:              ', (2.0_wp * pi * a)
                  write(*,*) 'dquad result:                   ', result
                  write(*,*) 'dquad ier:                      ', ier
                  write(*,*) 'rel err from elliptic function: ', (result - perimeter_of_ellipse(a,b)) / result
        if (a==b) write(*,*) 'rel err from analytical result: ', (result - (2.0_wp * pi * a)) / result

    end subroutine compare

    function perimeter_of_ellipse(a,b) result(p)

        real(wp),intent(in) :: a, b
        real(wp) :: p !! perimeter

        associate(a2 => a*a, b2 => b*b)
            p = 4.0_wp * a2 * drf(0.0_wp, b2, a2) &
                - (4.0_wp / 3.0_wp) * a2 * (a2 - b2) * drd(0.0_wp, b2, a2)
        end associate

    end function perimeter_of_ellipse

    real(wp) function f(x)
        !! 1d user function f(x)
        real(wp), intent(in) :: x !! this is dt
        real(wp),dimension(6) :: rvf !! final state vector
        call propagate(mu, rv0, x, rvf)
        f = norm2(rvf(4:6))  ! integrate vmag get arc length
    end function f

end program elliptical_perimeter_test