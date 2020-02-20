module angular_momentum_couplings
  !
  ! Most of functions defined here are originally written by the members at the Nuclear Theory group in the University of Tokyo
  !
  use basic_types
  use functions_from_c
  ! cache for Talmi-Moshinsky bracket
  integer(int32), private, parameter :: n_trinomial = 100, n_dbinomial=1000, n_triag=500
  real(real64), private, allocatable  :: dtrinomial(:,:,:)
  real(real64), private, allocatable  :: dbinomial(:,:), triangle_c(:,:,:)

  private :: dtrinomial_func
  private :: triangle
  private :: dbinomial_func
contains
  ! 3-j symbol for Wigner-Eckert theorem
  function geometry_part(jbra, jop, jket, mbra, mop, mket) result(r)
    integer(int32), intent(in) :: jbra, mbra, jop, mop, jket, mket
    real(real64) :: r
    r = (-1.d0) ** ((jbra - mbra)/2) * &
        & tjs(jbra, jop, jket, -mbra, mop, mket)
  end function geometry_part

  !!!
  ! For HO transformation bracket
  !!!
  function gmosh(nl, ll, nr, lr, n1, l1, n2, l2, lm, d) result(r)
    !
    ! reference G.P.Kamuntavicius, R.K.Kalinauskas, B.R.Barrett, S.Mickevicius, and D. Germanas, Nucl. Phys. A 695, 191 (2001).
    !
    real(real64) :: r
    integer(int32), intent(in) :: nl, ll, nr, lr, n1, l1, n2, l2, lm
    real(real64), intent(in) :: d
    integer(int32) :: ee, er, e1, e2, m, ed, eb, ec, ea, ld, lb, lc, la
    real(real64) :: t, s

    r = 0.d0
    if(.not. allocated(dtrinomial) ) then
      write(*,'(a)') "you need to call init_dtrinomial first!"
      return
    end if
    ee = 2*nl + ll
    er = 2*nr + lr
    e1 = 2*n1 + l1
    e2 = 2*n2 + l2
    if(ee + er /= e1 + e2) return
    if(triag(ll, lr, lm)) return
    if(triag(l1, l2, lm)) return
    t = dsqrt((d ** (e1 - er)) / ((1.d0 + d) ** (e1 + e2)))
    m = min(er, e2)
    s = 1.d0
    do ed = 0, m
      eb = er - ed
      ec = e2 - ed
      ea = e1 - er + ed

      do ld = ed, 0, -2
        do lb = eb, 0, -2
          if(triag(ld,lb,lr)) cycle
          do lc = ec, 0, -2
            if(triag(ld,lc,l2)) cycle
            do la = ea, 0, -2
              if(triag(la,lb,l1)) cycle
              if(triag(la,ll,lc)) cycle

              r = r + s * t * &
                  & ninej(2*la, 2*lb, 2*l1, 2*lc, 2*ld, 2*l2, 2*ll, 2*lr, 2*lm) * &
                  & g(e1, l1, ea, la, eb, lb) * g(e2, l2, ec, lc, ed, ld) * &
                  & g(ee, ll, ea, la, ec, lc) * g(er, lr, eb, lb, ed, ld)

            end do
          end do
        end do
      end do
      s = s * (-d)
    end do
    r = r * (-1.d0) ** (n1 + n2 + nr + nl)
  contains
    function g(e1, l1, ea, la, eb, lb) result(r)
      real(real64) :: r
      integer(int32), intent(in) :: e1, l1, ea, la, eb, lb

      r = cg_coef(2*la, 0, 2*lb, 0, 2*l1, 0) * dsqrt((2*la + 1) * (2*lb + 1) * &
          & dtrinomial(e1 - l1, ea - la, eb - lb) * &
          & dtrinomial(e1 + l1 + 1, ea + la + 1, eb + lb + 1))

    end function g

    function cg_coef(j1, m1, j2, m2, j3, m3) result(s)
      !
      !  Clebsch-Gordan coefficient
      !
      !  cg_coef(j1, m1, j2, m2, j3, m3)
      !  = ((j1)/2, (m1)/2, (j2)/2, (m2)/2 | (j3)/2, (m3)/2)
      !
      !  using the formula by Racah (1942)
      !
      implicit none
      integer(int32), intent(in) :: j1, j2, j3, m1, m2, m3
      integer(int32) :: js, jm1, jm2, jm3, k1, k2, k3
      integer(int32) :: iz, izmin, izmax, isign
      double precision :: s
      double precision :: tmp, delta

      s = 0.0d0
      jm1 = j1 - m1
      jm2 = j2 - m2
      jm3 = j3 - m3

      ! error and trivial-value check
      if (abs(m1)>j1 .or. abs(m2)>j2 .or. abs(m3)>j3 .or. mod(jm1,2)/=0 .or. mod(jm2,2)/=0 .or. mod(jm3,2)/=0) then
        write(*,*) 'error [cg_coef]: invalid j or m'
        write(*,'(1a, 1i4, 2x, 1a, 1i4)') 'j1 =', j1, 'm1 =', m1
        write(*,'(1a, 1i4, 2x, 1a, 1i4)') 'j2 =', j2, 'm2 =', m2
        write(*,'(1a, 1i4, 2x, 1a, 1i4)') 'j3 =', j3, 'm3 =', m3
        stop
      end if

      ! call triangle(j1, j2, j3, delta, info)
      ! if (info == 1) return
      if (max(j1, j2, j3) > n_triag) stop 'increase n_triag'
      delta = triangle_c( j1, j2, j3 )
      if (delta == 0.d0) return

      if (m3 /= m1+m2) return

      jm1 = jm1 / 2
      jm2 = jm2 / 2
      jm3 = jm3 / 2
      js = (j1 + j2 + j3)/2
      k1 = (j2 + j3 - j1)/2
      k2 = (j3 + j1 - j2)/2
      k3 = (j1 + j2 - j3)/2

      if (max(j1, j2, j3) > n_dbinomial) stop 'increase n_dbinomial'

      tmp = sqrt(dbinomial(j1, k2)/dbinomial(j1, jm1)) &
          & * sqrt(dbinomial(j2, k3)/dbinomial(j2, jm2)) &
          & * sqrt(dbinomial(j3, k1)/dbinomial(j3, jm3)) &
          & * sqrt((j3+1.0d0)) * delta

      izmin = max(0, jm1-k2, k3-jm2)
      izmax = min(k3, jm1, j2-jm2)

      if (izmax > n_dbinomial) stop 'increase n_dbinomial'

      isign = (-1)**izmin
      do iz = izmin, izmax
        s = s + isign * dbinomial(k3,iz) * dbinomial(k2,jm1-iz) &
            & * dbinomial(k1,j2-jm2-iz)
        isign = isign * (-1)
      end do

      s = s * tmp

    end function cg_coef

    function sixj(j1, j2, j3, l1, l2, l3) result(s)
      !
      !  6j coefficient
      !
      !  sixj(j1, j2, j3, l1, l2, l3) = {(j1)/2 (j2)/2 (j3)/2}
      !                                {(l1)/2 (l2)/3 (l3)/2}
      !
      !  see I. Talmi, Simple Models of Complex Nuclei, p. 158
      !
      implicit none
      integer(int32), intent(in) :: j1, j2, j3, l1, l2, l3
      double precision :: s
      double precision :: d
      integer(int32) :: izmin, izmax, iz, isign
      integer(int32) :: js, k1, k2, k3, jl1, jl2, jl3

      s = 0.0d0

      if (max(j1, j2, j3, l1, l2, l3) > n_triag) stop 'increase n_triag'

      d = triangle_c(j1, j2, j3)
      if (d == 0.d0) return

      d =    triangle_c(j1, l2, l3) * triangle_c(l1, j2, l3) &
          * triangle_c(l1, l2, j3) / d
      if ( d == 0.d0 ) return

      js = (j1 + j2 + j3)/2
      k1 = (j2 + j3 - j1)/2
      k2 = (j3 + j1 - j2)/2
      k3 = (j1 + j2 - j3)/2
      jl1 = (j1 + l2 + l3)/2
      jl2 = (l1 + j2 + l3)/2
      jl3 = (l1 + l2 + j3)/2

      izmin = max(0, js, jl1, jl2, jl3)
      izmax = min(k1+jl1, k2+jl2, k3+jl3)

      if (izmax+1 > n_dbinomial) stop 'increase n_dbinomial'

      isign = (-1)**izmin
      do iz = izmin, izmax
        s = s + isign * dbinomial(iz+1, iz-js) &
            & * dbinomial(k1, iz-jl1) * dbinomial(k2, iz-jl2) &
            & * dbinomial(k3, iz-jl3)
        isign = isign * (-1)
      end do
      s = s * d

    end function sixj

    function ninej(j11, j12, j13, j21, j22, j23, j31, j32, j33) result(s)
      !
      !  9j coefficient
      !
      !  ninej(j11, j12, j13, j21, j22, j23, j31, j32, j33)
      !
      !    {(j11)/2 (j12)/2 (j13)/2}
      !  = {(j21)/2 (j22)/2 (j23)/2}
      !    {(j31)/2 (j32)/2 (j33)/2}
      !
      !  see I. Talmi, Simple Models of Complex Nuclei, p. 968
      !
      implicit none
      integer(int32), intent(in) :: j11, j12, j13, j21, j22, j23, j31, j32, j33
      double precision :: s
      integer(int32) :: k, kmin, kmax

      kmin = max(abs(j11-j33), abs(j12-j23), abs(j21-j32))
      kmax = min(j11+j33, j12+j23, j21+j32)

      s = 0.0d0
      do k = kmin, kmax, 2
        s = s + (k+1.0d0) &
            & * sixj(j11, j12, j13, j23, j33, k) &
            & * sixj(j21, j22, j23, j12, k, j32) &
            & * sixj(j31, j32, j33, k, j11, j21)
      end do
      s = s * (-1)**kmin

    end function ninej

  end function gmosh

  subroutine init_dtrinomial()
    ! cache for Talmi-Moshinksy bracket
    integer(int32) :: i, j, k, n, m, info
    real(real64) :: d
    allocate(dtrinomial(0:n_trinomial, 0:n_trinomial, 0:n_trinomial))
    allocate( dbinomial(0:n_dbinomial, 0:n_dbinomial) )
    allocate(triangle_c(0:n_triag, 0:n_triag, 0:n_triag))
    dbinomial(:,:) = 0.d0

    !$omp parallel do private(n, m) schedule(dynamic)
    do n = 0, n_dbinomial
      do m = 0, n
        dbinomial(n, m) = dbinomial_func(n, m)
      end do
    end do

    triangle_c(:,:,:) = 0.d0
    !$omp parallel do private( i, j, k, d, info ) schedule (dynamic)
    do k = 0, n_triag
      do j = 0, n_triag
        do i = 0, n_triag
          call triangle(i, j, k, d, info)
          if (info == 1) then
            triangle_c(i, j, k) = 0.d0
            cycle
          end if
          triangle_c(i, j, k) = d
        end do
      end do
    end do

    !$omp parallel do private( i, j, k ) schedule (dynamic)
    do k = 0, n_trinomial
      do j = 0, n_trinomial
        do i = 0, n_trinomial
          dtrinomial(i, j, k) = dtrinomial_func(i, j, k)
        end do
      end do
    end do
  end subroutine init_dtrinomial

  subroutine fin_dtrinomial()
    deallocate(dtrinomial)
    deallocate(dbinomial)
    deallocate(triangle_c)
  end subroutine fin_dtrinomial


  function dtrinomial_func(i, j, k) result(s)
    !
    !  trinomial coefficient: i!! / j!! / k!!
    !
    integer(int32), intent(in) :: i, j, k
    real(real64) :: s
    integer(int32) :: m
    s = 1.d0
    m = max(i, j, k)
    if(m == 0) return
    if(m > n_trinomial) then

      write(*,'(a)') 'in trinomial_func, index is too large'
      return

    end if
    s = double_factorial(i) / (double_factorial(j) * double_factorial(k))
  end function dtrinomial_func

  subroutine triangle(j1, j2, j3, delta, info)
    !
    !  triangle often used in calculation of 3j, 6j etc.
    !  delta
    !  = sqrt(((j1+j2-j3)/2)!((j1-j2+j3)/2)!((-j1+j2+j3)/2)!/(1+(j1+j2+j3)/2)!)
    !
    implicit none
    integer(int32), intent(in) :: j1, j2, j3
    double precision, intent(out) :: delta
    integer(int32), intent(out) :: info
    integer(int32) :: js, k1, k2, k3

    info = 0
    js = j1 + j2 + j3
    k1 = j2 + j3 - j1
    k2 = j3 + j1 - j2
    k3 = j1 + j2 - j3

    if (j1 < 0 .or. j2 < 0 .or. j3 < 0) then
      write(*,*) 'error [triangle]: invalid j'
      write(*,'(1a, 1i4, 2x, 1a, 1i4, 2x, 1a, 1i4)') &
          & 'j1 =', j1, 'j2 =', j2, 'j3 =', j3
      stop
    end if
    if (k1 < 0 .or. k2 < 0 .or. k3 <0 .or. mod(js, 2) /=0) then
      info = 1
      return
    endif

    ! exclude too large arguments to prevent from round-off error
    if (js > n_triag*3) then
      write(*,'(1a, 1i5, 1a)') '[triangle]: j1+j2+j3 =', js, ' is too large'
      stop
    end if

    js = js / 2
    k1 = k1 / 2

    if (js > n_dbinomial) stop 'increase n_dbinomial'

    delta = 1.0d0 / &
        & (sqrt(dbinomial(js, j3)) * sqrt(dbinomial(j3, k1)) * sqrt(js+1.0d0))

  end subroutine triangle

  function dbinomial_func(n, m) result(s)
    !
    !  binomial coefficient: n_C_m
    !  s: double precision
    !
    integer(int32), intent(in) :: n, m
    double precision :: s, s1, s2
    integer(int32) :: i, m1

    s = 1.0d0
    m1 = min(m, n-m)
    if (m1 == 0) return
    if (n > 1000) then
      write(*,'(1a, 1i6, 1a)') '[dbinomial]: n =', n, ' is too large'
      stop
    end if

    if (n < 250) then
      s1 = 1.0d0
      s2 = 1.0d0
      do i = 1, m1
        s1 = s1 * (n-i+1)
        s2 = s2 * (m1-i+1)
      end do
      s = s1 / s2
    else
      do i = 1, m1
        s = (s * (n-i+1)) / (m1-i+1)
      end do
    endif

  end function dbinomial_func
  !!!
  ! end  HO transformation bracket
  !!!
  function triag(i,j,k)
    integer(int32),intent(in)::i,j,k
    logical :: triag
    triag = ((i-(j+k))*(i-abs(j-k)) > 0)
  end function triag

  function spherical_harmonics(l,m,cth,phi) result(r)
    !
    ! Y_{l}^{m}(\theta,\phi) = (-1)^{(m+|m|}/2} \sqrt{\frac{2l+1}{4\pi}}
    !            \times \sqrt{\frac{(l-m)!}{(l+m)!}} P_{l}^{m}(\cos(\theta)) e^{i m \phi}
    ! cth = \cos(\theta)
    !
    integer(int32), intent(in) :: l, m
    real(real64), intent(in) :: cth, phi
    complex(real64) :: ex, ei=(0.0_real64, 1.0_real64)
    complex(real64) :: r
    ex = cos(dble(m)*phi) + ei * sin(dble(m)*phi)
    if(m < 0) then
      r = assoc_legendre_spharm(l,m,cth) * ex
      return
    end if
    r = (-1.d0)**m * assoc_legendre_spharm(l,m,cth) * ex
  end function spherical_harmonics

end module angular_momentum_couplings
