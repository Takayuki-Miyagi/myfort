module iteration_methods
  !
  ! This module is used when one wants to find x* such that
  !
  !   x* = f(x*)  (1)
  !
  ! In the iterative method, Eq.(1) can be rewritten as
  !
  !   x^{out}_{n+1} = f(x^{in}_{n}) (2)
  !
  ! The calculation is done repeatedly until || x^{out}_{n+1} - x^{out}_{n} || < epsilon is satisfied.
  ! The class "optimizer" can be used to get x^{in}_{n}
  ! User can chose methods:
  !   direct iteration x^{in}_{n} = x^{out}_{n}
  !   linear mixing    x^{in}_{n} = a * x^{out}_{n} + (1-a) * x^{in}_{n-1}
  !   Broyden method
  !   BFGS method
  !   modified Broyden method
  !   L-BFGS method
  ! See also:
  !   https://en.wikipedia.org/wiki/Quasi-Newton_method (Broyden method, BFGS method)
  !   A. Baran et al., Phys. Rev. C 78, 014318 (2008). (Modified Broyden method)
  !   https://en.wikipedia.org/wiki/Limited-memory_BFGS (L-BFGS)
  !
  use linear_algebra
  implicit none

  public :: Optimizer
  private :: InitOptimizer
  private :: FinOptimizer
  private :: SetInitialCondition
  private :: GetNextInput
  private :: DirectIteration
  private :: LinearMixing
  private :: BroydenMethod
  private :: ModifiedBroydenMethod
  private :: BFGSMethod
  private :: LBFGSMethod
  private :: get_norm_diff

  type :: Optimizer
    integer(int32) :: n = 0    ! dimension of vector (problem size)
    integer(int32) :: nite = 0 ! number of iteration
    real(real64) :: r = 1.d0 ! abs error
    real(real64) :: a = 1.d0 ! mixing ratio (0 < a <= 1)
    integer(int32) :: m = 10   ! number of previous step taken into account
    type(DMat) :: Bj, Bi
    type(DVec) :: xj, fj, xi, fi
    type(DVec), allocatable :: df(:), dx(:)
    character(40) :: Method = "direct"
  contains
    procedure :: InitOptimizer
    procedure :: FinOptimizer
    procedure :: SetInitialCondition
    procedure :: GetNextInput
    procedure :: DirectIteration
    procedure :: BroydenMethod
    procedure :: LinearMixing
    procedure :: ModifiedBroydenMethod
    procedure :: BFGSMethod
    procedure :: LBFGSMethod

    generic :: init => InitOptimizer
    generic :: fin => Finoptimizer
    generic :: get_next => GetNextInput
    generic :: set_init => SetInitialCondition
  end type Optimizer
contains

  subroutine InitOptimizer(this, n, a, method, m)
    class(Optimizer), intent(inout) :: this
    integer(int32), intent(in), optional :: n, m
    character(*), intent(in), optional :: method
    real(real64), intent(in), optional :: a
    integer(int32) :: i
    if(present(n)) this%n = n
    if(present(m)) this%m = m
    if(present(method)) this%method = method
    if(present(a)) this%a = a

    if(this%n == 0) then
      write(*,"(a,i8)") "Error: problem size n=", this%n
      stop
    end if
    call this%xi%zeros(this%n)
    call this%xj%zeros(this%n)
    call this%fi%zeros(this%n)
    call this%fj%zeros(this%n)

    select case(this%method)
    case("Broyden", "broyden", "BFGS", "bfgs")
      call this%Bi%zeros(this%n,this%n)
      call this%Bj%zeros(this%n,this%n)
      do i = 1, this%n
        this%Bi%m(i,i) = - this%a
        this%Bj%m(i,i) = - this%a
      end do
    case default
    end select

    select case(this%method)
    case("ModifiedBroyden", "modifiedbroyden", "mbroyden", &
          & "LBFGS", "lbfgs", "L-BFGS", "l-bfgs")
      allocate(this%df(this%m))
      allocate(this%dx(this%m))
      do i = 1, this%m
        call this%df(i)%zeros(this%n)
        call this%dx(i)%zeros(this%n)
      end do
    case default
    end select
!    write(*,*)
!    write(*,"(3a)") "## Iteration method is ", trim(this%method), " method"
!    write(*,*)
  end subroutine InitOptimizer

  subroutine FinOptimizer(this)
    class(Optimizer), intent(inout) :: this
    integer(int32) :: i
    call this%xi%fin()
    call this%xj%fin()
    call this%fi%fin()
    call this%fj%fin()
    select case(this%method)
    case("Broyden", "broyden", "BFGS", "bfgs")
      call this%Bi%fin()
      call this%Bj%fin()
    case default
    end select
    select case(this%method)
    case("ModifiedBroyden", "modifiedbroyden", "mbroyden", &
          & "LBFGS", "lbfgs", "L-BFGS", "l-bfgs")
      do i = 1, this%m
        call this%df(i)%fin()
        call this%dx(i)%fin()
      end do
      deallocate(this%df)
      deallocate(this%dx)
    case default
    end select
    this%r = 1.d0
  end subroutine FinOptimizer

  subroutine SetInitialCondition(this, x)
    class(Optimizer), intent(inout) :: this
    real(real64), intent(in) :: x(:)
    this%nite = 1
    this%xi%v = x
  end subroutine SetInitialCondition

  subroutine GetNextInput(this, v, ite)
    class(Optimizer), intent(inout) :: this
    real(real64), intent(inout) :: v(:)
    integer(int32), intent(in) :: ite

    this%nite = ite
    this%fj%v = v

    select case(this%method)
    case("direct")
      call this%DirectIteration()
    case("linearmixing", "linear", "lm")
      call this%LinearMixing()
    case("Broyden", "broyden")
      call this%BroydenMethod()
    case("BFGS", "bfgs")
      call this%BFGSMethod()
    case("ModifiedBroyden", "modifiedbroyden", "mbroyden")
      call this%ModifiedBroydenMethod()
    case("LBFGS", "L-BFGS", "lbfgs", "l-bfgs")
      call this%LBFGSMethod()
    case default
      write(*,"(a)") "Unkown method, direct iteration is assumed"
      call this%DirectIteration()
    end select
!    write(*,"(es18.8)") this%r
    v = this%xj%v
  end subroutine GetNextInput

  subroutine DirectIteration(this)
    class(Optimizer), intent(inout) :: this
    this%r = get_norm_diff(this%fj, this%xj)
    this%xj = this%fj
  end subroutine DirectIteration

  subroutine LinearMixing(this)
    class(Optimizer), intent(inout) :: this
    this%r = get_norm_diff(this%fj, this%xj)
    this%xj = this%a * this%fj + (1.d0 - this%a) * this%xj
  end subroutine LinearMixing

  subroutine BroydenMethod(this)
    ! find x* so that x* = f(x*)
    ! x_{j+1} = x_{j} - alpha B_{j} f_{j}
    ! B_{j+1} = B_{j} + (|dx> - Bj|df>) <dx|Bj / <dx|Bj|df>
    ! xj : input of j-th iteration
    ! fj : f(xj)
    ! Bj : Jacobian of j-th iteration
    class(Optimizer), intent(inout) :: this
    type(DVec) :: df, dx
    real(real64) :: a
    this%r = get_norm_diff(this%fj, this%xj)
    this%fj = (this%fj - this%xj)
    if(this%nite == 1) then
      this%xj = this%xi - ((this%Bi * this%fj) * this%a)
      this%fi = this%fj
      return
    else
      df = this%fj - this%fi
      dx = this%xj - this%xi
      a = 1.d0 / (dx * (this%Bi * df))
      this%Bj = this%Bi + ((dx - (this%Bi * df)) .x. (dx * this%Bi)) * a
      this%xi = this%xj
      this%fi = this%fj
      this%xj = this%xi - (this%Bj * this%fj)
      this%Bi = this%Bj
      call df%Fin()
      call dx%Fin()
      return
    end if
  end subroutine BroydenMethod

  subroutine BFGSMethod(this)
    ! find x* so that x* = f(x*)
    ! x_{j+1} = x_{j} - alpha B_{j} f_{j}
    ! B_{j+1} = (I - |df><dx| / <df|dx>)^T Bj (I - |df><dx| / <df|dx>) + |dx> <dx| / <df|dx>
    ! xj : input of j-th iteration
    ! fj : f(xj)
    ! Bj : Jacobian of j-th iteration
    class(Optimizer), intent(inout) :: this
    type(DVec) :: df, dx
    type(DMat) :: I, m1, m2, m3
    real(real64) :: a
    this%r = get_norm_diff(this%fj, this%xj)
    this%fj = this%fj - this%xj
    if(this%nite == 1) then
      this%xj = this%xi - (this%Bi * this%fj)
      this%fi = this%fj
      return
    else
      df = this%fj - this%fi
      dx = this%xj - this%xi
      a = 1.d0 / (df * dx)
      call I%eye(this%n)
      m1 = I - ((df .x. dx) * a)
      m2 = m1%T() * this%Bi
      m3 = (dx .x. dx) * a
      this%Bj = m2 + m3
      this%xi = this%xj
      this%fi = this%fj
      this%xj = this%xi - (this%Bj * this%fj)
      this%Bi = this%Bj
      call df%Fin()
      call dx%Fin()
      call m1%Fin()
      call m2%Fin()
      call m3%Fin()
      call I%Fin()
      return
    end if
  end subroutine BFGSMethod

  subroutine ModifiedBroydenMethod(this)
    ! find x* so that x* = f(x*)
    ! x_{j+1} = x_{j} - alpha B_{j} f_{j}
    ! xj : input of j-th iteration
    ! fj : f(xj)
    ! Bj : Jacobian of j-th iteration
    class(Optimizer), intent(inout) :: this
    type(DMat) :: beta, mat, u
    type(DVec) :: vec, gamm, c
    real(real64) :: a, curv, w0
    integer(int32) :: iused, inext, icurr, i, j

    this%r = get_norm_diff(this%fj, this%xj)
    this%fj = (this%fj - this%xj)
    if(this%m == 0 .or. this%nite == 0) then
      this%xj = this%xi + this%fj * this%a * 1.d-1
      return
    end if
    iused = min(this%nite, this%m)
    icurr = mod(this%nite-1, this%m) + 1
    inext = mod(icurr, this%m) + 1
    w0 = 1.d-2
    if(this%nite > 1) then
      this%df(icurr) = this%fj - this%df(icurr)
      this%dx(icurr) = this%xj - this%dx(icurr)
      a = 1.d0 / dsqrt((this%df(icurr) * this%df(icurr)))
      this%df(icurr) = this%df(icurr) * a
      this%dx(icurr) = this%dx(icurr) * a
    end if
    call mat%zeros(iused, iused)
    call beta%zeros(iused, iused)
    do i = 1, iused
      do j = 1, iused
        mat%m(i,j) = (this%df(i) * this%df(j))
      end do
      mat%m(i,i) = 1.d0 + w0*w0
    end do
    if(iused > 0) beta = mat%inv()
    call c%ini(iused)
    call u%ini(this%n, iused)
    do i = 1, iused
      c%v(i) = this%df(i) * this%fj
      u%m(:,i) = (this%df(i)%v(:) * this%a) + this%dx(i)%v(:)
    end do
    gamm = c * beta
    vec = (this%fj * this%a) - (u * gamm)
    curv = this%fj * vec
    this%df(inext) = this%fj
    this%dx(inext) = this%xj
    if(curv > -1.d0) then
      this%xj = this%xj + vec
    else
      this%xj = this%xj + this%fj * (this%a * 0.5)
    end if
    call mat%Fin()
    call vec%Fin()
    call beta%Fin()
    call gamm%Fin()
  end subroutine ModifiedBroydenMethod

  subroutine LBFGSMethod(this)
    ! find x* so that x* = f(x*)
    class(Optimizer), intent(inout) :: this
    type(DMat) :: beta, mat, u
    type(DVec) :: vec, gamm, c
    real(real64) :: a, curv, w0
    integer(int32) :: iused, inext, icurr, i, j

    this%r = get_norm_diff(this%fj, this%xj)
    this%fj = (this%fj - this%xj)
    if(this%m == 0 .or. this%nite == 0) then
      this%xj = this%xi + this%fj * this%a
      return
    end if

    iused = min(this%nite, this%m)
    icurr = mod(this%nite-1, this%m) + 1
    inext = mod(icurr, this%m) + 1
    w0 = 1.d-2
    if(this%nite > 1) then
      this%df(icurr) = this%fj - this%df(icurr)
      this%dx(icurr) = this%xj - this%dx(icurr)
      a = 1.d0 / dsqrt((this%df(icurr) * this%df(icurr)))
      this%df(icurr) = this%df(icurr) * a
      this%dx(icurr) = this%dx(icurr) * a
    end if
    call mat%ini(iused, iused)
    call beta%ini(iused, iused)
    do i = 1, iused
      do j = 1, iused
        mat%m(i,j) = (this%df(i) * this%df(j))
      end do
      mat%m(i,i) = 1.d0 + w0*w0
    end do
    if(iused > 0) beta = mat%inv()
    call c%ini(iused)
    call u%ini(this%n, iused)
    do i = 1, iused
      c%v(i) = this%df(i) * this%fj
      u%m(:,i) = (this%df(i)%v(:) * this%a) + this%dx(i)%v(:)
    end do
    gamm = c * beta
    vec = (this%fj * this%a) - (u * gamm)
    curv = this%fj * vec
    this%df(inext) = this%fj
    this%dx(inext) = this%xj
    if(curv > -1.d0) then
      this%xj = this%xj + vec
    else
      this%xj = this%xj + this%fj * (this%a * 0.5)
    end if
    call mat%Fin()
    call vec%Fin()
    call beta%Fin()
    call gamm%Fin()
  end subroutine LBFGSMethod

  function get_norm_diff(a, b) result(c)
    type(DVec), intent(in) :: a, b
    real(real64) :: c
    integer(int32) :: i
    c = 0.d0
    do i = 1, a%n_size
      c = c + sqrt( (a%v(i) - b%v(i))**2 )
    end do
  end function get_norm_diff

end module iteration_methods
