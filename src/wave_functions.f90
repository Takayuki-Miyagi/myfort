module wave_functions
  use functions_from_c
  use physics_constants
  implicit none
contains
  function ho_radial_wf(n,l,anu,r) result(s)
    ! R = sqrt( 2 * nu**3 * Gamma(n+1) / Gamma(n+l+1.5) ) x^{l} e^{ -x^2/2 } L^{n}_{l+0.5}( x^2 )
    ! x = r * nu or p * nu
    ! nu = sqrt(m w / h) or sqrt( h / mw )
    ! Note that anu = nu^2
    integer(int32),intent(in) :: n,l
    real(real64),intent(in) :: anu,r
    real(real64) :: s
    real(real64) :: prefact, exp_component, nu, x, a
    nu = sqrt(anu)
    x = nu * r
    a = dble(l)+0.5d0
    prefact = sqrt( 2.d0 * nu**3 )
    exp_component = 0.5d0*ln_gamma(dble(n+1)) - 0.5d0*ln_gamma(dble(n+l)+1.5d0) - 0.5d0*x**2 &
        & + dble(l) * log(x)
    s = prefact * exp(exp_component) * laguerre(n,a,x**2)
  end function ho_radial_wf

  function ho_radial_wf_norm(n,l,anu,r) result(s)
    ! R = sqrt( 2 * nu * Gamma(n+1) / Gamma(n+l+1.5) ) x^{l+1} e^{ -x^2/2 } L^{n}_{l+0.5}( x^2 )
    ! x = r * nu or p * nu
    ! nu = sqrt(m w / h) or sqrt( h / mw )
    ! Note that anu = nu^2
    integer(int32),intent(in) :: n,l
    real(real64),intent(in) :: anu,r
    real(real64) :: s
    real(real64) :: prefact, exp_component, nu, x, a
    nu = sqrt(anu)
    x = nu * r
    a = dble(l)+0.5d0
    prefact = sqrt( 2.d0 * nu )
    exp_component = 0.5d0*ln_gamma(dble(n+1)) - 0.5d0*ln_gamma(dble(n+l)+1.5d0) - 0.5d0*x**2 &
        & + dble(l+1) * log(x)
    s = prefact * exp(exp_component) * laguerre(n,a,x**2)
  end function ho_radial_wf_norm

  function hydrogen_radial_wf(n,l,a_star,r) result(s)
    ! R_{nl} = sqrt( (2/n a* )**3 (n-l-1)!/2n(n+l)! ) x^{l} e^{ -x/2 } L^{2l+1}_{n-l-1}( x )
    ! x = 2 * r / (n*a*)
    ! a* is reduced Bohr radius
    integer(int32),intent(in) :: n,l
    real(real64),intent(in) :: a_star,r
    real(real64) :: s
    real(real64) :: prefact, exp_component,  x
    x = 2.d0 * r / ( a_star * dble(n))
    prefact = sqrt( 4.d0/ a_star**3 ) / dble(n**2)
    exp_component = 0.5d0*ln_gamma(dble(n-l)) - 0.5d0*ln_gamma(dble(n+l+1)) - 0.5d0*x + dble(l)*log(x)
    s = prefact * exp(exp_component) * laguerre(n-l-1,dble(2*l+1),x)
  end function hydrogen_radial_wf

  function hydrogen_radial_wf_norm(n,l,a_star,r) result(s)
    ! R_{nl} = sqrt( (2/n a* ) (n-l-1)!/2n(n+l)! ) x^{l+1} e^{ -x/2 } L^{2l+1}_{n-l-1}( x )
    ! x = 2 * r / (n*a*)
    ! a* is reduced Bohr radius
    integer(int32),intent(in) :: n,l
    real(real64),intent(in) :: a_star,r
    real(real64) :: s
    real(real64) :: prefact, exp_component,  x
    x = 2.d0 * r / ( a_star * dble(n))
    prefact = sqrt( 1.d0/a_star ) / dble(n)
    exp_component = 0.5d0*ln_gamma(dble(n-l)) - 0.5d0*ln_gamma(dble(n+l+1)) - 0.5d0*x + dble(l+1)*log(x)
    s = prefact * exp(exp_component) * laguerre(n-l-1,dble(2*l+1),x)
  end function hydrogen_radial_wf_norm

  function Laguerre_radial_wf(n,l,b,r) result(s)
    ! From Eq.(13) in A. E. McCoy and M. A. Caprio, J. Math. Phys. 57, (2016).
    ! b is length scale [L]
    integer,intent(in) :: n,l
    real(8),intent(in) :: b,r
    real(8) :: s
    real(8) :: prefact, exp_component, x, x2
    x = r / b
    x2 = 2.d0 * x
    prefact = sqrt( 2.d0 / b ) * 2.d0 / b
    exp_component = 0.5d0*ln_gamma(dble(n+1)) - &
        & 0.5d0*ln_gamma(dble(n+2*l+3)) - x
    if(abs(x2) > 1.d-16) then
      exp_component = exp_component  + dble(l)*log(x2)
    end if
    s = prefact * exp(exp_component) * laguerre(n,dble(2*l+2),x2)
  end function Laguerre_radial_wf

  function Laguerre_radial_wf_norm(n,l,b,r) result(s)
    ! From Eq.(13) in A. E. McCoy and M. A. Caprio, J. Math. Phys. 57, (2016).
    ! b is length scale [L]
    integer(int32),intent(in) :: n,l
    real(real64),intent(in) :: b,r
    real(real64) :: s
    real(real64) :: prefact, exp_component, x, x2
    x = r / b
    x2 = 2.d0 * x
    prefact = sqrt( 2.d0 / b )
    exp_component = 0.5d0*ln_gamma(dble(n+1)) - &
        & 0.5d0*ln_gamma(dble(n+2*l+3)) - x  + dble(l+1) * log(x2)
    s = prefact * exp(exp_component) * laguerre(n,dble(2*l+2),x2)
  end function Laguerre_radial_wf_norm

  function Mom_Laguerre_radial_wf_norm(n,l,b,p) result(s)
    ! b is length scale [L]
    integer,intent(in) :: n,l
    real(8),intent(in) :: b,p
    real(8) :: s
    real(8) :: exp_prefact, exp_component, x, xi
    integer :: k
    xi = b * p * 0.5d0
    x = 0.5d0 / sqrt(0.25d0 + xi*xi)
    s = 0.d0
    exp_prefact = 0.5d0*ln_gamma(dble(n+1)) + 0.5d0*ln_gamma(dble(n+2*l+3)) + &
        & dble(2*l+3)*log(2.d0*x) + dble(l)*log(2.d0*xi) + ln_gamma(dble(l+1))
    do k = 0, n
      exp_component = dble(k)*log(2.d0*x) - ln_gamma(dble(n+1-k)) - ln_gamma(dble(2*l+k+3)) + exp_prefact
      s = s + (-1.d0)**k * dble(k+1) * Gegenbauer_polynomial(k+1,dble(l+1),x) * exp(exp_component)
    end do
    s = s * sqrt(b / pi) * xi
  end function Mom_Laguerre_radial_wf_norm

  function hydrogen_radial_wf_mom_norm(n,l,a_star,p) result(s)
    ! from wikipedia https://en.wikipedia.org/wiki/Hydrogen_atom
    ! Do not forget phase (-i)^{l} when you calclate the matrix elements!!
    ! x = p * (n*a*)
    ! a* is reduced Bohr radius
    integer(int32),intent(in) :: n, l
    real(real64),intent(in) :: a_star, p
    real(real64) :: s
    real(real64) :: prefact, exp_component, x, z
    x = p * dble(n) * a_star
    z = (x**2-1.d0) / (x**2+1.d0)
    prefact = sqrt( 2.d0/pi * a_star) * dble(n)
    exp_component = 0.5d0*ln_gamma(dble(n-l)) - 0.5d0*ln_gamma(dble(n+l+1)) + dble(2*l+2)*log(2.d0) + &
        & ln_gamma(dble(l+1)) + dble(l+1)*log(x) - dble(l+2)*log(x**2+1)
    s = prefact * exp(exp_component) * Gegenbauer_polynomial(n-l-1,dble(l+1),z)
  end function hydrogen_radial_wf_mom_norm

end module wave_functions
