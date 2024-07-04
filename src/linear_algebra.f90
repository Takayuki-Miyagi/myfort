module linear_algebra
  use vector_definitions
  use matrix_definitions

  implicit none
  public :: EigenSolSymD
  public :: EigenSolHermite
  public :: exp


  ! diagonalization methods
  private :: InitEigenSolSymD
  private :: FinEigenSolSymD
  private :: DiagSymD
  private :: EigenvalD
  private :: InitEigenSolHermite
  private :: FinEigenSolHermite
  private :: DiagHermite
  !private :: EigenvalHermite

  interface exp
    module procedure :: ExpD, ExpC
  end interface exp

  type :: EigenSolSymD
    type(DVec) :: eig
    type(DMat) :: vec
  contains
    procedure :: init => InitEigenSolSymD
    procedure :: fin => FinEigenSolSymD
    procedure :: DiagSym => DiagSymD   ! eigen values and eigen vectors
    procedure :: Eigenval => EigenvalD ! only eigen values
  end type EigenSolSymD

  type :: GenEigenSolSymD
    type(DVec) :: eig
    type(DMat) :: vec
    integer :: itype = 1
  contains
    procedure :: init => InitGenEigenSolSymD
    procedure :: fin => FinGenEigenSolSymD
    procedure :: DiagSym => DiagGenSymD   ! eigen values and eigen vectors
  end type GenEigenSolSymD

  type :: EigenSolHermite
    type(DVec) :: eig
    type(CMat) :: vec
  contains
    procedure :: init => InitEigenSolHermite
    procedure :: fin => FinEigenSolHermite
    procedure :: DiagSym => DiagHermite      ! eigen values and eigen vectors
    !procedure :: Eigenval => EigenvalHermite ! only eigen values
  end type EigenSolHermite
contains

  subroutine InitEigenSolSymD(this, A)
    class(EigenSolSymD) :: this
    type(DMat), intent(in) :: A
    integer(int32) :: n
    n = size(A%m, 1)
    call this%eig%ini(n)
    call this%vec%ini(n,n)
  end subroutine InitEigenSolSymD

  subroutine FinEigenSolSymD(this)
    class(EigenSolSymD) :: this
    call this%eig%fin()
    call this%vec%fin()
  end subroutine FinEigenSolSymD

  subroutine DiagSymD(this, A, qmin, qmax, m, error)
    class(EigenSolSymD) :: this
    type(DMat), intent(in) :: A
    real(real64), intent(in), optional :: qmin, qmax
    integer(int32), intent(in), optional :: m
    integer(int32), intent(in), optional :: error
    real(real64), allocatable :: work(:), rcondz(:), zerrbd(:), mat(:,:)
    integer(int32), allocatable :: iwork(:), ifailv(:)
    integer(int32) :: info, lwork, n, i, num
    real(real64) :: lw, dlamch, e, eerbd
    n = size(A%M, 1)
    this%vec = A

    if(.not. present(m) .and. .not. present(qmin) .and. .not. present(qmax)) then
      !
      ! solve all eigen values and eigen vectors
      !

      allocate(work(1))
      call dsyev('v', 'u', n, this%vec%m, n, this%eig%v, work, -1, info)
      lwork = int(work(1))
      deallocate(work)

      allocate(work(lwork))
      call dsyev('v', 'u', n, this%vec%m, n, this%eig%v, work, lwork, info)
      if(info /= 0) then
        write(*,'(a, i6)') 'Error in DiagSym: info = ', info
        stop
      end if
      deallocate(work)
      do i = 1, n
        if(this%vec%m(1,i) < 0.d0) this%vec%m(:,i) = this%vec%m(:,i) * (-1.d0)
      end do
      if(present(error)) then
        allocate(rcondz(n), zerrbd(n))
        e = epsilon(1.d0)
        eerbd = e * max(abs(this%eig%v(1)), abs(this%eig%v(n)))
        call ddisna('Eigenvectors', n, n, this%eig%v, rcondz, info)
        do i = 1, n
          zerrbd(i) = eerbd / rcondz(i)
        end do
        write(*,'(a)') 'Error estimate for eigen values'
        write(*, '(1es12.4)') eerbd
        write(*,'(a)') 'Error estimate for eigen vectors'
        write(*, '(10es12.4)') zerrbd
        deallocate(rcondz, zerrbd)
      end if

    elseif(present(m)) then
      !
      ! solve lowest m eigen values and eigen vectors
      !
      this%eig%v(:) = 0.d0
      allocate(mat(n,n))
      allocate(iwork(5*n), ifailv(n))
      mat = A%m
      allocate(work(1))
      call dsyevx('v', 'i', 'u', n, mat, n, -1.d100, 1.d100, 1, m, dlamch('S'), &
        &  num, this%eig%v, this%vec%m, n, work, -1, iwork, ifailv, info)
      lwork = int(work(1))
      deallocate(work)

      allocate(work(1:lwork))
      call dsyevx('v', 'i', 'u', n, mat, n, -1.d100, 1.d100, 1, m, dlamch('S'), &
          &  num, this%eig%v, this%vec%m, n, work, lwork, iwork, ifailv, info)
      this%vec%m(:,num+1:n) = 0.d0
      deallocate( iwork, ifailv, work, mat)
      do i = 1, num
        if(this%vec%m(1,i) < 0.d0) this%vec%m(:,i) = this%vec%m(:,i) * (-1.d0)
      end do

    elseif(present(qmin) .and. present(qmax)) then
      !
      ! solve eigen values in (qmin, qmax) and
      ! corrsponding eigen vectors
      !

      allocate(mat(n,n))
      allocate(iwork(5*n), ifailv(n))
      mat = A%m
      this%eig%v(:) = 0.d0
      allocate(work(1))
      call dsyevx('v', 'v', 'u', n, mat, n, qmin, qmax, 1, n, dlamch('S'), &
        &  num, this%eig%v, this%vec%m, n, work, -1, iwork, ifailv, info)
      lwork = int(work(1))
      deallocate(work)
      allocate(work(1:lwork))
      call dsyevx('v', 'v', 'u', n, mat, n, qmin, qmax, 1, n, dlamch('S'), &
          &  num, this%eig%v, this%vec%m, n, work, lwork, iwork, ifailv, info)
      this%vec%m(:,num+1:n) = 0.d0
      deallocate( iwork, ifailv, work, mat )
      do i = 1, num
        if(this%vec%m(1,i) < 0.d0) this%vec%m(:,i) = this%vec%m(:,i) * (-1.d0)
      end do
    end if
  end subroutine DiagSymD

  subroutine EigenvalD(this, A, m)
    class(EigenSolSymD) :: this
    type(DMat), intent(in) :: A
    integer(int32), intent(in) :: m
    integer(int32), allocatable :: iwork(:), iblock(:), isplit(:)
    real(real64), allocatable :: work(:), d(:), e(:), tau(:), w(:)
    real(real64) :: dlamch, lw
    integer(int32) :: n, info, lwork, nsplit, i

    n = size(A%M, 1)
    allocate(d(n), e(max(1, n-1)), tau(max(1, n-1)), w(n))
    allocate(iblock(n), isplit(n))
    allocate(work(1))
    call dsytrd('u',n,A%m,n,d,e,tau,work,-1,info)
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))
    call dsytrd('u',n,A%m,n,d,e,tau,work,lwork,info)
    call dstebz('i','e',n,0.d0,0.d0,1,min(n,m),dlamch('S'), &
        & d, e, m, nsplit, w, iblock, isplit, work, iwork, info)
    do i = 1, min(n,m)
      this%eig%v(i) = w(i)
    end do
    deallocate(work)
    deallocate(d, e, tau,iblock, isplit)
  end subroutine EigenvalD

  subroutine InitGenEigenSolSymD(this, A, B, itype)
    class(GenEigenSolSymD) :: this
    type(DMat), intent(in) :: A, B
    integer(int32), intent(in), optional :: itype
    integer(int32) :: n
    if(present(itype)) this%itype = itype
    n = size(A%m, 1)
    if(this%itype == 2) n = size(A%m, 1)
    if(this%itype == 3) n = size(B%m, 1)
    call this%eig%ini(n)
    call this%vec%ini(n,n)
  end subroutine InitGenEigenSolSymD

  subroutine FinGenEigenSolSymD(this)
    class(GenEigenSolSymD) :: this
    call this%eig%fin()
    call this%vec%fin()
  end subroutine FinGenEigenSolSymD

  subroutine DiagGenSymD(this, A, B)
    class(GenEigenSolSymD) :: this
    type(DMat), intent(in) :: A, B
    integer(int32) :: n, lda, ldb, lwork, liwork, info, idummy
    real(real64) :: dummy
    integer(int32), allocatable :: iwork(:)
    real(real64), allocatable :: work(:)
    this%vec = A
    n = size(A%M, 1)
    lda = size(A%m,1)
    ldb = size(B%m,1)
    allocate(work(1), iwork(1))
    call dsygvd(this%itype, 'V', 'U', n, A%m, lda, B%m, ldb, this%eig%v, work, -1, iwork, -1, info)
    lwork = int(work(1))
    liwork = iwork(1)
    deallocate(work,iwork)
    allocate(work(lwork), iwork(liwork))
    call dsygvd(this%itype, 'V', 'U', n, A%m, lda, B%m, ldb, this%eig%v, work, lwork, iwork, liwork, info)
    deallocate(work, iwork)
    this%vec = A
  end subroutine DiagGenSymD

  subroutine InitEigenSolHermite(this, A)
    class(EigenSolHermite) :: this
    type(CMat), intent(in) :: A
    integer(int32) :: n
    n = size(A%m, 1)
    call this%eig%ini(n)
    call this%vec%ini(n,n)
  end subroutine InitEigenSolHermite

  subroutine FinEigenSolHermite(this)
    class(EigenSolHermite) :: this
    call this%eig%fin()
    call this%vec%fin()
  end subroutine FinEigenSolHermite

  !subroutine DiagHermite(this, A, qmin, qmax, m, error)
  subroutine DiagHermite(this, A, qmin, qmax, m)
    class(EigenSolHermite) :: this
    type(CMat), intent(in) :: A
    real(real64), intent(in), optional :: qmin, qmax
    integer(int32), intent(in), optional :: m
    !integer(int32), intent(in), optional :: error
    complex(real64), allocatable :: work(:)
    real(real64), allocatable :: rwork(:)
    integer(int32) :: info, lwork, n
    real(real64) :: lw
    n = size(A%M, 1)
    this%vec = A

    if(.not. present(m) .and. .not. present(qmin) .and. .not. present(qmax)) then

      !
      ! solve all eigen values and eigen vectors
      !

      allocate(rwork(3*n - 2))
      allocate(work(1))
      call zheev('v', 'u', n, this%vec%m, n, this%eig%v, work, -1, rwork, info)
      lwork = int(work(1))
      deallocate(work)

      allocate(work(lwork))
      call zheev('v', 'u', n, this%vec%m, n, this%eig%v, work, lwork, rwork, info)
      if(info /= 0) then
        write(*,'(a, i6)') 'Error in DiagSym: info = ', info
        stop
      end if
      deallocate(work, rwork)
      !if(present(error)) then
      !  allocate(rcondz(n), zerrbd(n))
      !  e = epsilon(1.d0)
      !  eerbd = e * max(abs(this%eig%v(1)), abs(this%eig%v(n)))
      !  call ddisna('Eigenvectors', n, n, this%eig%v, rcondz, info)
      !  do i = 1, n
      !    zerrbd(i) = eerbd / rcondz(i)
      !  end do
      !  write(*,'(a)') 'Error estimate for eigen values'
      !  write(*, '(1es12.4)') eerbd
      !  write(*,'(a)') 'Error estimate for eigen vectors'
      !  write(*, '(10es12.4)') zerrbd
      !  deallocate(rcondz, zerrbd)
      !end if

      !elseif(present(m)) then
      !  !
      !  ! solve lowest m eigen values and eigen vectors
      !  !
      !  this%eig%v(:) = 0.d0
      !  allocate(mat(n,n))
      !  allocate(iwork(5*n), ifailv(n))
      !  mat = A%m
      !  call dsyevx('v', 'i', 'u', n, mat, n, -1.d100, 1.d100, 1, m, dlamch('S'), &
      !      &  num, this%eig%v, this%vec%m, n, lw, -1, iwork, ifailv, info)
      !  lwork = int(lw)
      !  allocate(work(1:lwork))
      !  call dsyevx('v', 'i', 'u', n, mat, n, -1.d100, 1.d100, 1, m, dlamch('S'), &
      !      &  num, this%eig%v, this%vec%m, n, work, lwork, iwork, ifailv, info)
      !  this%vec%m(:,num+1:n) = 0.d0
      !  deallocate( iwork, ifailv, work, mat)

      !elseif(present(qmin) .and. present(qmax)) then
      !  !
      !  ! solve eigen values in (qmin, qmax) and
      !  ! corrsponding eigen vectors
      !  !

      !  allocate(mat(n,n))
      !  allocate(iwork(5*n), ifailv(n))
      !  mat = A%m
      !  this%eig%v(:) = 0.d0
      !  call dsyevx('v', 'v', 'u', n, mat, n, qmin, qmax, 1, n, dlamch('S'), &
      !      &  num, this%eig%v, this%vec%m, n, lw, -1, iwork, ifailv, info)
      !  lwork = int(lw)
      !  allocate(work(1:lwork))
      !  call dsyevx('v', 'v', 'u', n, mat, n, qmin, qmax, 1, n, dlamch('S'), &
      !      &  num, this%eig%v, this%vec%m, n, work, lwork, iwork, ifailv, info)
      !  this%vec%m(:,num+1:n) = 0.d0
      !  deallocate( iwork, ifailv, work, mat )
    end if
  end subroutine DiagHermite

  !subroutine EigenvalHermite(this, A, m)
  !  class(EigenSolHermite) :: this
  !  type(DMat), intent(in) :: A
  !  integer(int32), intent(in) :: m
  !  integer(int32), allocatable :: iwork(:), iblock(:), isplit(:)
  !  real(real64), allocatable :: work(:), d(:), e(:), tau(:), w(:)
  !  real(real64) :: dlamch, lw
  !  integer(int32) :: n, info, lwork, nsplit, i

  !  n = size(A%M, 1)
  !  allocate(d(n), e(max(1, n-1)), tau(max(1, n-1)), w(n))
  !  allocate(iblock(n), isplit(n))
  !  call dsytrd('u',n,A%m,n,d,e,tau,lw,-1,info)
  !  lwork = int(lw)
  !  allocate(work(lwork))
  !  call dsytrd('u',n,A%m,n,d,e,tau,work,lwork,info)
  !  call dstebz('i','e',n,0.d0,0.d0,1,min(n,m),dlamch('S'), &
  !      & d, e, m, nsplit, w, iblock, isplit, work, iwork, info)
  !  do i = 1, min(n,m)
  !    this%eig%v(i) = w(i)
  !  end do
  !  deallocate(work)
  !  deallocate(d, e, tau,iblock, isplit)
  !end subroutine EigenvalHermite

  type(DMat) function ExpD(a, ord, tol_in, show_process) result(r)
    type(DMat), intent(in) :: a
    type(DMat) :: b
    integer(int32), intent(in), optional :: ord
    real(real64), intent(in), optional :: tol_in
    logical, intent(in), optional :: show_process
    real(real64) :: tol = 1.d-8
    logical :: show = .false.
    integer(int32) :: i
    integer(int32) :: iord = 12
    if(present(ord)) iord = ord
    if(present(tol_in)) tol = tol_in
    if(present(show_process)) show = show_process
    call r%eye(size(a%m, 1))
    b = r
    do i = 1, iord
      b = (b * a) / dble(i)
      r = r + b
      if(show) then
        write(*,"(a,i4,a,es18.8)") "max value of matrix: order=",i,&
            & " maxval(A)",maxval(b%m)
      end if
      if((maxval(b%m)) < tol) return
    end do
    write(*,"(a,i4)") "warning: increase ord, current value is ", iord
  end function ExpD

  type(CMat) function ExpC(a, ord) result(r)
    type(CMat), intent(in) :: a
    type(CMat) :: b
    integer(int32), intent(in), optional :: ord
    integer(int32) :: i
    integer(int32) :: iord = 12
    if(present(ord)) iord = ord
    call r%eye(size(a%m, 1))
    b = r
    do i = 1, iord
      b = b * a / (cmplx(i,kind=real64))
      r = r + b
    end do
  end function ExpC
end module linear_algebra
