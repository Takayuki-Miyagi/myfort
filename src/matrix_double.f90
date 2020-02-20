module matrix_double
  use basic_types
  implicit none
  public :: DMat
  public :: MatrixCopyD
  public :: MatrixProductD
  public :: MatrixSumD
  public :: MatrixSubtractD
  public :: MatrixScaleLD
  public :: MatrixScaleRD
  public :: MatrixScaleDivideD

  private :: IniM
  private :: FinM
  private :: zeros
  private :: eye
  private :: Trans
  private :: Inverse
  private :: Det
  private :: GetRandomMatrix
  private :: MatrixPrint
  private :: DiagMat
  private :: block_dmat


  type :: DMat
    real(real64), allocatable :: M(:,:)
    integer(int32) :: n_row = 0
    integer(int32) :: n_col = 0
  contains
    procedure :: Ini => IniM
    procedure :: Fin => FinM
    procedure :: zeros
    procedure :: eye
    procedure :: T => Trans
    procedure :: Inv => Inverse
    procedure :: Det
    procedure :: Random => GetRandomMatrix
    procedure :: blk => block_dmat
    procedure :: prt => MatrixPrint
    procedure :: DiagMat
  end type DMat
contains
  subroutine iniM(a, m, n)
    class(DMat), intent(inout) :: a
    integer(int32), intent(in) :: m,n
    if(allocated(a%m)) call a%fin()
    if(m < 1 .or. n < 1) return
    a%n_row=m
    a%n_col=n
    allocate(a%m(a%n_row,a%n_col))
  end subroutine iniM

  subroutine zeros(a, m, n)
    class(DMat), intent(inout) :: a
    integer(int32), intent(in) :: m,n
    if(allocated(a%m)) call a%fin()
    if(m < 1 .or. n < 1) return
    call a%ini(m,n)
    a%m = 0.d0
  end subroutine zeros

  subroutine eye(a, n)
    class(DMat), intent(inout) :: a
    integer(int32), intent(in) :: n
    integer(int32) :: i
    if(allocated(a%m)) call a%fin()
    if(n < 1) return
    call a%zeros(n,n)
    !$omp parallel
    !$omp do private(i)
    do i = 1, n
      a%m(i,i) = 1.d0
    end do
    !$omp end do
    !$omp end parallel
  end subroutine eye

  subroutine FinM(a)
    class(DMat), intent(inout) :: a
    if(allocated(a%m)) deallocate(a%m)
    a%n_row=0
    a%n_col=0
  end subroutine FinM

  subroutine MatrixCopyD(b, a)
    type(DMat), intent(inout) :: b
    type(DMat), intent(in) :: a
    integer(int32) :: m, n, i
    m = a%n_row
    n = a%n_col
    if(m < 1 .or. n < 1) return
    call b%Ini(m,n)
    !$omp parallel
    !$omp do private(i)
    do i = 1, n
      call dcopy(m, a%m(:,i), 1, b%m(:,i), 1)
    end do
    !$omp end do
    !$omp end parallel
  end subroutine MatrixCopyD

  type(DMat) function MatrixProductD(a, b) result(c)
    type(DMat), intent(in) :: a, b
    integer(int32) :: m, k, n
    m = a%n_row
    k = a%n_col
    if(a%n_col /= b%n_row) then
      write(*, '(a)') 'Error in MatrixProduct'
      stop
    end if
    n = b%n_col
    if(m < 1 .or. n < 1) return
    call c%Ini(m,n)
    call dgemm('n','n',m,n,k,1.d0,a%m,m,b%m,k,0.d0,c%m,m)
  end function MatrixProductD

  type(DMat) function MatrixSumD(a, b) result(c)
    type(DMat), intent(in) :: a, b
    integer(int32) :: m, n, i
    if(a%n_row /= b%n_row .or. a%n_col /= b%n_col) then
      write(*, '(a)') 'Error in MatrixSumD'
      stop
    end if
    m = a%n_row
    n = a%n_col
    if(m < 1 .or. n < 1) return
    call MatrixCopyD(c, a)
    !$omp parallel
    !$omp do private(i)
    do i = 1, n
      call daxpy(m, 1.d0, b%m(:,i), 1, c%m(:,i), 1)
    end do
    !$omp end do
    !$omp end parallel
  end function MatrixSumD

  type(DMat) function MatrixSubtractD(a, b) result(c)
    type(DMat), intent(in) :: a, b
    integer(int32) :: m, n, i
    if(a%n_row /= b%n_row .or. a%n_col /= b%n_col) then
      write(*, '(a)') 'Error in MatrixSubtractD'
      stop
    end if
    m = a%n_row
    n = a%n_col
    if(m < 1 .or. n < 1) return
    call MatrixCopyD(c, a)
    !$omp parallel
    !$omp do private(i)
    do i = 1, n
      call daxpy(m, -1.d0, b%m(:,i), 1, c%m(:,i), 1)
    end do
    !$omp end do
    !$omp end parallel
  end function MatrixSubtractD

  type(DMat) function MatrixScaleLD(b, a) result(c)
    type(DMat), intent(in) :: b
    real(real64), intent(in) :: a
    integer(int32) :: m, n, i
    m = b%n_row
    n = b%n_col
    if(m < 1 .or. n < 1) return
    call MatrixCopyD(c, b)
    !$omp parallel
    !$omp do private(i)
    do i = 1, n
      call dscal(m, a, c%m(:,i), 1)
    end do
    !$omp end do
    !$omp end parallel
  end function MatrixScaleLD

  type(DMat) function MatrixScaleRD(a, b) result(c)
    type(DMat), intent(in) :: b
    real(real64), intent(in) :: a
    integer(int32) :: m, n, i
    m = b%n_row
    n = b%n_col
    if(m < 1 .or. n < 1) return
    call MatrixCopyD(c, b)
    !$omp parallel
    !$omp do private(i)
    do i = 1, n
      call dscal(m, a, c%m(:,i), 1)
    end do
    !$omp end do
    !$omp end parallel
  end function MatrixScaleRD

  type(DMat) function MatrixScaleDivideD(b, a) result(c)
    type(DMat), intent(in) :: b
    real(real64), intent(in) :: a
    integer(int32) :: m, n, i
    m = b%n_row
    n = b%n_col
    if(m < 1 .or. n < 1) return
    call MatrixCopyD(c, b)
    !$omp parallel
    !$omp do private(i)
    do i = 1, n
      call dscal(m, 1.d0 / a, c%m(:,i), 1)
    end do
    !$omp end do
    !$omp end parallel
  end function MatrixScaleDivideD

  type(DMat) function Trans(a) result(b)
    class(DMat), intent(in) :: a
    integer(int32) :: n, m
    m = a%n_row
    n = a%n_col
    if(m < 1 .or. n < 1) return
    call b%Ini(n,m)
    b%M = transpose(a%M)
  end function Trans

  type(DMat) function inverse(r) result(s)
    class(DMat), intent(in) :: r
    real(real64), allocatable :: a(:,:)
    real(real64), allocatable :: work(:)
    integer(int32), allocatable :: ipvt(:)
    integer(int32) :: info, n
    n = r%n_row
    if(n < 1) return
    call s%Ini(n,n)
    allocate(work(n*n),ipvt(n))
    allocate(a(n,n))
    a = r%m
    call dgetrf(n,n,a,n,ipvt,info)
    call dgetri(n,a,n,ipvt,work,n**2,info)
    s%m = a
    deallocate(a,work,ipvt)
  end function inverse

  real(real64) function Det(r) result(d)
    class(DMat), intent(in) :: r
    integer(int32) :: n, i, info
    real(real64), allocatable :: a(:,:)
    integer(int32), allocatable :: ipiv(:)

    d = 0.d0
    n = r%n_row
    if(n < 1) return
    allocate(ipiv(n), a(n,n))
    a = r%m
    call dgetrf(n, n, a, n, ipiv, info)
    if(info /= 0) then
      write(*,'(a, i3)') "error in det: info = ", info
      stop
    end if
    d = 1.d0
    do i = 1, n
      if(ipiv(i) .ne. i) then
        d = -d * a(i, i)
      else
        d = d * a(i, i)
      end if
    end do
    deallocate(ipiv, a)
  end function Det

  subroutine MatrixPrint(this, msg, iunit, binary)
    class(DMat), intent(in) :: this
    integer, intent(in), optional :: iunit
    character(*), intent(in), optional :: msg
    logical, intent(in), optional :: binary
    logical :: bin
    character(12) :: cfmt
    integer(int32) :: i, j, n, m
    integer :: unt

    if(this%n_row<1 .or. this%n_col<1) return

    if(present(iunit)) then; unt = iunit
    else; unt = 6; end if

    if(present(binary)) then; bin = binary
    else; bin = .false.; end if

    if(bin) then
      write(unt) this%m
    else
      if(present(msg)) write(unt,*) msg
      n = this%n_row
      m = this%n_col
      if(unt == 6) then
        cfmt = '( xf10.4)'
        write(cfmt(2:3), '(I2)') m
        do i=1,n
          write(unt,cfmt) this%m(i,:)
        end do
      else
        do i=1,n
          do j=1,m
            write(unt,'(2i8,f14.6)') i,j,this%m(i,j)
          end do
        end do
      end if
    end if
  end subroutine MatrixPrint

  function block_dmat(this, m1, m2, n1, n2) result(r)
    class(DMat), intent(in) :: this
    type(DMat) :: r
    integer(int32), intent(in) :: m1, m2, n1, n2
    integer(int32) :: m, n
    m = m2 - m1 + 1
    n = n2 - n1 + 1
    if(m < 1 .or. n < 1) return
    call r%ini(m,n)
    r%m(:,:) = this%m(m1:m2,n1:n2)
  end function block_dmat

  subroutine GetRandomMatrix(mat, m, n)
    use vector_double, only: DVec
    class(DMat), intent(inout) :: mat
    integer(int32), intent(in) :: m, n
    integer(int32) :: i
    type(DVec) :: v
    call mat%ini(m,n)
    if(m < 1 .or. n < 1) return
    do i = 1, n
      call v%Random(m)
      mat%m(:,i) = v%v(:)
      call v%fin()
    end do
  end subroutine GetRandomMatrix

  subroutine DiagMat(b, a)
    use vector_double, only: DVec
    class(DMat), intent(inout) :: b
    type(DVec), intent(in) :: a
    integer(int32) :: n, i
    n = a%n_size
    if(n < 1) return
    call b%zeros(n,n)
    !$omp parallel
    !$omp do private(i)
    do i = 1, n
      b%M(i,i) = a%V(i)
    end do
    !$omp end do
    !$omp end parallel
  end subroutine DiagMat
end module matrix_double
