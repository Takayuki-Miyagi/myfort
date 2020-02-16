module physics_constants
  use general
  real(real64), parameter, public :: pi = 3.141592741012573d0 ! \pi
  real(real64), parameter, public :: hc = 197.32705d0         ! \hbar c [MeV fm] or [eV nm]
  real(real64), public :: amp = 938.27231d0        ! proton mass [MeV] can be changed in LQCD calc.
  real(real64), public :: amn = 939.56563d0        ! neutron mass [MeV] can be changed in LQCD calc.
  real(real64), parameter, public :: m_e = 510.9989461 ! electron mass [keV]
  real(real64), parameter, public :: alpha = 137.035999d0     ! electric fine structure constant
  real(real64) :: g_a = 1.29d0        ! axial vector coupling
  real(real64) :: f_pi = 92.4d0       ! pion decay constatnt [MeV]
  real(real64) :: m_pi = 138.d0       ! pion mass [MeV]
  real(real64) :: gs = 0.880          ! nucleon's magnetic moemnt g-factor isoscalar
  real(real64) :: gv = 4.706          ! nucleon's magnetic moment g-factor isovector  mu_{p/n} = (gs +/- gv)/2
end module physics_constants
