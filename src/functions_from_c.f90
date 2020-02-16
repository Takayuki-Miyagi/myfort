module functions_from_c
  use, intrinsic :: iso_c_binding
  implicit none
  !
  ! C interfaces
  !
  interface
    ! 3-j symbol
    function coupling_3j(j1,j2,j3,m1,m2,m3) bind(c,name='gsl_sf_coupling_3j')
      import c_int, c_double
      real(c_double) :: coupling_3j
      integer(c_int), value, intent(in) :: j1,j2,j3,m1,m2,m3
    end function coupling_3j

    ! 6-j symbol
    function coupling_6j(j1,j2,j3,j4,j5,j6) bind(c,name='gsl_sf_coupling_6j')
      import c_int, c_double
      real(c_double) :: coupling_6j
      integer(c_int), value, intent(in) :: j1,j2,j3,j4,j5,j6
    end function coupling_6j

    ! 9-j symbol
    function coupling_9j(j1,j2,j3,j4,j5,j6,j7,j8,j9) bind(c,name='gsl_sf_coupling_9j')
      import c_int, c_double
      real(c_double) :: coupling_9j
      integer(c_int), value, intent(in) :: j1,j2,j3,j4,j5,j6,j7,j8,j9
    end function coupling_9j

    ! factorial n!
    function factorial(n) bind(c,name='gsl_sf_fact')
      import c_double, c_int
      integer(c_int), value, intent(in) :: n
      real(c_double) :: factorial
    end function factorial

    ! double factorial n!!
    function double_factorial(n) bind(c,name='gsl_sf_doublefact')
      import c_double, c_int
      integer(c_int), value, intent(in) :: n
      real(c_double) :: double_factorial
    end function double_factorial

    ! Gamma function \Gamma(x)
    function gamma_function(x) bind(c,name='gsl_sf_gamma')
      import c_double
      real(c_double), value, intent(in) :: x
      real(c_double) :: gamma_function
    end function gamma_function

    ! log of Gamma function ln \Gamma(x)
    function ln_gamma(x) bind(c,name='gsl_sf_lngamma')
      import c_double
      real(c_double), value, intent(in) :: x
      real(c_double) :: ln_gamma
    end function ln_gamma

    ! spherical bessel function j_l(x)
    function spherical_bessel(l,x) bind(c,name='gsl_sf_bessel_jl')
      import c_int, c_double
      integer(c_int), value, intent(in) :: l
      real(c_double), value, intent(in) :: x
      real(c_double) :: spherical_bessel
    end function spherical_bessel

    ! Legendre polynomial P_l(x)
    function legendre_polynomial(l,x) bind(c,name='gsl_sf_legendre_Pl')
      import c_int, c_double
      integer(c_int), value, intent(in) :: l
      real(c_double), value, intent(in) :: x
      real(c_double) :: legendre_polynomial
    end function legendre_polynomial

    ! Associated Legendre polynomial sutable for pherical harmonics
    ! \sqrt{\frac{2l+1}{4\pi}} \sqrt{\frac{(l-m)!}{(l+m)!}} P_{l}^{m}(x)
    function assoc_legendre_spharm(l,m,x) bind(c,name='gsl_sf_legendre_sphPlm')
      import c_int, c_double
      integer(c_int), value, intent(in) :: l,m
      real(c_double), value, intent(in) :: x
      real(c_double) :: assoc_legendre_spharm
    end function assoc_legendre_spharm

    ! associated Laguerre polynomial L^{(a)}_{n}(x)
    function laguerre(n,a,x) bind(c,name='gsl_sf_laguerre_n')
      import c_int, c_double
      integer(c_int), value, intent(in) :: n
      real(c_double), value, intent(in) :: a, x
      real(c_double) :: laguerre
    end function laguerre

    ! Gegenbauer polynomial C^{(lambda)}_{n}(x)
    function Gegenbauer_polynomial(n,lambda,x) bind(c,name='gsl_sf_gegenpoly_n')
      import c_int, c_double
      integer(c_int), value, intent(in) :: n
      real(c_double), value, intent(in) :: lambda, x
      real(c_double) :: gegenbauer_polynomial
    end function Gegenbauer_polynomial

    ! Gauss-Legendre quadrature
    function gauss_legendre_allocate(n) &
          & bind(c,name='gsl_integration_glfixed_table_alloc')
      import c_int, c_ptr
      integer(c_int), value, intent(in) :: n
      type(c_ptr) :: gauss_legendre_allocate
    end function gauss_legendre_allocate
    function gauss_legendre_ith_point_weight(a,b,i,xi,wi,t) &
          & bind(c,name='gsl_integration_glfixed_point')
      import c_int, c_double, c_ptr
      real(c_double), value, intent(in) :: a, b
      integer(c_int), value, intent(in) :: i
      real(c_double) :: xi, wi
      type(c_ptr), value, intent(in) :: t
      integer(c_int) :: gauss_legendre_ith_point_weight
    end function gauss_legendre_ith_point_weight
    subroutine gauss_legendre_release(t) &
          & bind(c,name='gsl_integration_glfixed_table_free')
      import c_ptr
      type(c_ptr), value :: t
    end subroutine gauss_legendre_release

    ! fixed-point quadrature
    function integration_fixed_allocate(T, n, a, b, alpha, beta) &
          & bind(c,name='gsl_integration_fixed_alloc')
      import c_int, c_double, c_ptr
      integer(c_int), value, intent(in) :: n
      real(c_double), value, intent(in) :: a, b, alpha, beta
      type(c_ptr), value :: T
      type(c_ptr) :: integration_fixed_allocate
    end function integration_fixed_allocate
    subroutine integration_fixed_free( workspace ) bind(c, name="gsl_integration_fixed_free")
      import c_ptr
      type(c_ptr), value :: workspace
    end subroutine integration_fixed_free
    function integration_fixed_n( workspace ) bind(c, name="gsl_integration_fixed_n")
      import c_ptr, c_int
      type(c_ptr), value :: workspace
      integer(c_int) :: integration_fixed_n
    end function integration_fixed_n
    function integration_fixed_nodes( workspace ) bind(c, name="gsl_integration_fixed_nodes")
      import c_ptr
      type(c_ptr), value :: workspace
      type(c_ptr) :: integration_fixed_nodes
    end function integration_fixed_nodes
    function integration_fixed_weights( workspace ) bind(c, name="gsl_integration_fixed_weights")
      import c_ptr
      type(c_ptr), value :: workspace
      type(c_ptr) :: integration_fixed_weights
    end function integration_fixed_weights

    ! open, read, write, and close gzip file (additional interface gzip_open below)
    ! When you pass string to C, you need add NULL (achar(0)) in the end of strings.
    ! It is done in gzip_open.
    function gz_open(filename, mode) bind(c, name='gzopen')
      import c_char, c_ptr
      character(c_char) :: filename(*), mode(*)
      type(c_ptr) :: gz_open
    end function gz_open
    function gzip_read(f, buf, len) bind(c, name='gzgets')
      import c_int, c_char, c_ptr
      type(c_ptr), value :: f
      character(c_char) :: buf(*)
      integer(c_int), value, intent(in) :: len
      type(c_ptr) :: gzip_read
    end function gzip_read
    function gzip_write(f, buf, len) bind(c, name='gzwrite')
      import c_int, c_char, c_ptr
      type(c_ptr), value :: f
      character(c_char) :: buf(*)
      integer(c_int), value, intent(in) :: len
      type(c_ptr) :: gzip_write
    end function gzip_write
    function gzip_close( f ) bind(c, name='gzclose')
      import c_ptr
      type(c_ptr), value :: f
      type(c_ptr) :: gzip_close
    end function gzip_close
  end interface
  type(c_ptr) :: integration_fixed_legendre; bind(c, name="gsl_integration_fixed_legendre") :: integration_fixed_legendre
  type(c_ptr) :: integration_fixed_chebyshev; bind(c, name="gsl_integration_fixed_chebyshev") :: integration_fixed_chebyshev
  type(c_ptr) :: integration_fixed_gegenbauer; bind(c, name="gsl_integration_fixed_gegenbauer") :: integration_fixed_gegenbauer
  type(c_ptr) :: integration_fixed_jacobi; bind(c, name="gsl_integration_fixed_jacobi") :: integration_fixed_jacobi
  type(c_ptr) :: integration_fixed_laguerre; bind(c, name="gsl_integration_fixed_laguerre") :: integration_fixed_laguerre
  type(c_ptr) :: integration_fixed_hermite; bind(c, name="gsl_integration_fixed_hermite") :: integration_fixed_hermite
  type(c_ptr) :: integration_fixed_exponential; bind(c, name="gsl_integration_fixed_exponential") :: integration_fixed_exponential
  type(c_ptr) :: integration_fixed_rational; bind(c, name="gsl_integration_fixed_rational") :: integration_fixed_rational
  type(c_ptr) :: integration_fixed_chebyshev2; bind(c, name="gsl_integration_fixed_chebyshev2") :: integration_fixed_chebyshev2
  !
  ! end C interfaces
  !
contains

  function dcg(j1, m1, j2, m2, j3, m3) result(s)
    !
    !  Clebsch-Gordan coefficient
    !  dcg(j1, m1, j2, m2, j3, m3)
    !  = ((j1)/2, (m1)/2, (j2)/2, (m2)/2 | (j3)/2, (m3)/2)
    integer, intent(in) :: j1, j2, j3, m1, m2, m3
    real(8) :: s
    s = coupling_3j(j1,j2,j3,m1,m2,-m3) * sqrt(dble(j3+1)) * (-1.d0) ** ((j1-j2+m3)/2)
  end function dcg

  function tjs(j1, j2, j3, m1, m2, m3) result(r)
    real(8) :: r
    integer, intent(in) :: j1, j2, j3, m1, m2, m3
    r = coupling_3j(j1,j2,j3,m1,m2,m3)
  end function tjs

  function sjs(j1, j2, j3, l1, l2, l3) result(s)
    !
    !  6j coefficient
    !  d6j(j1, j2, j3, l1, l2, l3) = {(j1)/2 (j2)/2 (j3)/2}
    !                                {(l1)/2 (l2)/3 (l3)/2}
    integer, intent(in) :: j1, j2, j3, l1, l2, l3
    real(8) :: s
    s = coupling_6j(j1,j2,j3,l1,l2,l3)
  end function sjs

  function snj(j11, j12, j13, j21, j22, j23, j31, j32, j33) result(s)
    !
    !  9j coefficient
    !  d9j(j11, j12, j13, j21, j22, j23, j31, j32, j33)
    !    {(j11)/2 (j12)/2 (j13)/2}
    !  = {(j21)/2 (j22)/2 (j23)/2}
    !    {(j31)/2 (j32)/2 (j33)/2}
    integer, intent(in) :: j11, j12, j13, j21, j22, j23, j31, j32, j33
    real(8) :: s
    s = coupling_9j(j11,j12,j13,j21,j22,j23,j31,j32,j33)
  end function snj

  subroutine gauss_legendre(x1,x2,x,w,n)
    ! input:
    ! x1   : lower limit of the integration interval
    ! x2   : upper limit ---------- "" -------------
    ! n    : the desired number of mesh points
    ! output :
    ! x     : gauss-legendre mesh points on the interval (x1,x2)
    ! w     : the corresponding weights
    integer, intent(in) :: n
    real(8), intent(in) :: x1, x2
    real(8), intent(out), allocatable :: x(:), w(:)
    real(8) :: xi, wi
    integer :: info, i
    type(c_ptr) :: t

    if(allocated(x)) deallocate(x)
    if(allocated(w)) deallocate(w)
    allocate(x(n))
    allocate(w(n))
    t = gauss_legendre_allocate(n)
    do i = 1, n
      info = gauss_legendre_ith_point_weight(x1,x2,i-1,xi,wi,t)
      x(i) = xi
      w(i) = wi
    end do
    call gauss_legendre_release(t)
  end subroutine gauss_legendre

  subroutine gauss_legendre_(x1,x2,x,w,n)
    ! input:
    ! x1   : lower limit of the integration interval
    ! x2   : upper limit ---------- "" -------------
    ! n    : the desired number of mesh points
    ! output :
    ! x     : gauss-legendre mesh points on the interval (x1,x2)
    ! w     : the corresponding weights
    integer, intent(in) :: n
    real(8), intent(in) :: x1, x2
    real(8), intent(out), allocatable :: x(:), w(:)
    call fixed_point_quadrature("legendre", n, x, w, a_in=x1, b_in=x2)
  end subroutine gauss_legendre_

  subroutine fixed_point_quadrature(quad_name, n, x, w, a_in, b_in, alpha_in, beta_in, weight_renorm)
    character(*), intent(in) :: quad_name
    integer, intent(in) :: n
    real(8), intent(in), optional :: a_in, b_in, alpha_in, beta_in
    logical, intent(in), optional :: weight_renorm
    real(8), allocatable :: x(:), w(:)
    real(8) :: a=0.d0, b=0.d0, alpha=0.d0, beta=0.d0
    integer :: i
    type(c_ptr) :: workspace, nodes, weights
    real(c_double), pointer :: x_(:), w_(:)

    if(allocated(x)) deallocate(x)
    if(allocated(w)) deallocate(w)
    allocate(x(n))
    allocate(w(n))
    if(present(a_in)) a = a_in
    if(present(b_in)) b = b_in
    if(present(alpha_in)) alpha = alpha_in
    if(present(beta_in)) beta = beta_in

    select case(quad_name)
    case("legendre","Legendre")
      workspace = integration_fixed_allocate(integration_fixed_legendre, n, a, b, alpha, beta)
    case("chebyshev1","Chebyshev1","Chebyshev type 1", "Chebyshev Type 1")
      workspace = integration_fixed_allocate(integration_fixed_chebyshev, n, a, b, alpha, beta)
    case("gegenbauer","Gegenbauer")
      workspace = integration_fixed_allocate(integration_fixed_gegenbauer, n, a, b, alpha, beta)
    case("jacobi","Jacobi")
      workspace = integration_fixed_allocate(integration_fixed_jacobi, n, a, b, alpha, beta)
    case("laguerre","Laguerre")
      workspace = integration_fixed_allocate(integration_fixed_laguerre, n, a, b, alpha, beta)
    case("hermite","Hermite")
      workspace = integration_fixed_allocate(integration_fixed_hermite, n, a, b, alpha, beta)
    case("exponential","Exponential")
      workspace = integration_fixed_allocate(integration_fixed_exponential, n, a, b, alpha, beta)
    case("rational","Rational")
      workspace = integration_fixed_allocate(integration_fixed_rational, n, a, b, alpha, beta)
    case("chebyshev2","Chebyshev2","Chebyshev type 2", "Chebyshev Type 2")
      workspace = integration_fixed_allocate(integration_fixed_chebyshev2, n, a, b, alpha, beta)
    case default
      write(*,"(a)") "Unknown quadrature name"
      stop
    end select

    nodes = integration_fixed_nodes( workspace )
    weights = integration_fixed_weights( workspace )
    call c_f_pointer(nodes, x_, [n])
    call c_f_pointer(weights, w_, [n])
    x(:) = x_(:)
    w(:) = w_(:)
    call integration_fixed_free(workspace)
    if(.not. present(weight_renorm)) return
    if(.not. weight_renorm) return
    do i = 1, n
      select case(quad_name)
      case("legendre","Legendre")
        w(i) = w(i) * 1.d0
      case("chebyshev1","Chebyshev1","Chebyshev type 1", "Chebyshev Type 1")
        w(i) = w(i) * sqrt((b-x(i)) * (x(i)-a))
      case("gegenbauer","Gegenbauer")
        w(i) = w(i) / ( (b-x(i)) * (x(i)-a) )**alpha
      case("jacobi","Jacobi")
        w(i) = w(i) / ( (b-x(i))**alpha * (x(i)-a)**beta )
      case("laguerre","Laguerre")
        w(i) = w(i) * exp( b * (x(i)-a) ) / (x(i)-a)**alpha
      case("hermite","Hermite")
        w(i) = w(i) * exp( b * (x(i)-a)**2 ) / abs(x(i)-a)**alpha
      case("exponential","Exponential")
        w(i) = w(i) / abs( x(i)-(a+b)*0.5d0 )**alpha
      case("rational","Rational")
        w(i) = w(i) / ( (x(i)-a)**alpha * (x(i)+b)*beta )
      case("chebyshev2","Chebyshev2","Chebyshev type 2", "Chebyshev Type 2")
        w(i) = w(i) / sqrt( (b-x(i)) * (x(i)-a) )
      case default
        write(*,"(a)") "Unknown quadrature name"
        stop
      end select
    end do
  end subroutine fixed_point_quadrature

  function gzip_open( filename, mode ) result(p)
    character(*), intent(in) :: filename, mode
    type(c_ptr) :: p
    p = gz_open(trim(filename)//achar(0), trim(mode)//achar(0))
  end function gzip_open

  function gzip_writeline( f, buf, len) result(p)
    type(c_ptr) :: f, p
    character(*), intent(in) :: buf
    integer, intent(in) :: len
    p = gzip_write( f, trim(buf)//achar(10), len+1)
  end function gzip_writeline

  function gzip_readline( f, buf, len) result(p)
    ! note
    ! len_trim returns length removed space (32 in ascii code)
    ! -2 means removing null (0) and line feed (10)
    type(c_ptr) :: f, p
    character(*), intent(inout) :: buf
    integer, intent(in) :: len
    p = gzip_read( f, buf, len)
    buf = buf(1:len_trim(buf) - 2)
  end function gzip_readline
end module functions_from_c

