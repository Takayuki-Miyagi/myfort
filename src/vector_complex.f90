module vector_complex
  use basic_types
  implicit none

  public :: CVec
  public :: VectorCopyC
  public :: VectorSumC
  public :: VectorSubtractC
  public :: VectorScaleRC
  public :: VectorScaleLC
  public :: VectorDivideC
  public :: InnerProductC

  private :: iniV
  private :: FinV
  private :: zeros
  private :: ComplexConjugate
  private :: VectorPrint
  private :: GetRandomVector
  private :: Nrm
  private :: Nrm2
  private :: block_cvec

  type :: CVec
    complex(real64), allocatable :: v(:)
    integer(int32) :: n_size = 0
  contains
    procedure :: Ini => iniV
    procedure :: zeros
    procedure :: Fin => FinV
    procedure :: CC => ComplexConjugate
    procedure :: prt => VectorPrint
    procedure :: Random => GetRandomVector
    procedure :: blk => block_cvec
    procedure :: Nrm
    procedure :: Nrm2
  end type CVec
contains
  subroutine IniV(a, n)
    class(CVec), intent(inout) :: a
    integer(int32), intent(in) :: n
    if(n < 1) return
    if(allocated(a%V)) deallocate(a%V)
    a%n_size = n
    allocate(a%V(a%n_size))
  end subroutine IniV

  subroutine zeros(a, n)
    class(CVec), intent(inout) :: a
    integer(int32), intent(in) :: n
    if(allocated(a%V)) deallocate(a%V)
    if(n < 1) return
    a%n_size = n
    allocate(a%V(a%n_size))
    a%V(:) = (0.d0,0.d0)
  end subroutine zeros

  subroutine FinV(a)
    class(CVec), intent(inout) :: a
    if(allocated(a%V)) deallocate(a%V)
    a%n_size = 0
  end subroutine FinV

  subroutine VectorCopyC(b, a)
    type(CVec), intent(inout) :: b
    type(CVec), intent(in) :: a
    integer(int32) :: n
    n = a%n_size
    if(n < 1) return
    call b%Ini(n)
    call zcopy(n, a%v, 1, b%v, 1)
  end subroutine VectorCopyC

  type(CVec) function VectorSumC(a, b) result(c)
    type(CVec), intent(in) :: a, b
    integer(int32) :: n
    if(a%n_size /= b%n_size) then
      write(*,'(a)') 'Error in DVectorSum'
      stop
    end if
    n = a%n_size
    if(n < 1) return
    call VectorCopyC(c,a)
    call zaxpy(n, (1.d0,0.d0), b%v, 1, c%v, 1)
  end function VectorSumC

  type(CVec) function VectorSubtractC(a, b) result(c)
    type(CVec), intent(in) :: a, b
    integer(int32) :: n
    if(a%n_size /= b%n_size) then
      write(*,'(a)') 'Error in DVectorSubtract'
      stop
    end if
    n = a%n_size
    if(n < 1) return
    call VectorCopyC(c,a)
    call zaxpy(n, (-1.d0,0.d0), b%v, 1, c%v, 1)
  end function VectorSubtractC

  type(CVec) function VectorScaleRC(a, b) result(c)
    type(CVec), intent(in) :: a
    complex(real64), intent(in) :: b
    integer(int32) :: n
    n = a%n_size
    if(n < 1) return
    call VectorCopyC(c,a)
    call zscal(n, b, c%v, 1)
  end function VectorScaleRC


  type(CVec) function VectorScaleLC(b, a) result(c)
    type(CVec), intent(in) :: a
    complex(real64), intent(in) :: b
    integer(int32) :: n
    n = a%n_size
    if(n < 1) return
    call VectorCopyC(c,a)
    call zscal(n, b, c%v, 1)
  end function VectorScaleLC

  type(CVec) function VectorDivideC(a, b) result(c)
    type(CVec), intent(in) :: a
    real(real64), intent(in) :: b
    complex(real64) :: bb
    integer(int32) :: n
    n = a%n_size
    if(n < 1) return
    call VectorCopyC(c,a)
    bb = 1.d0/b
    call zscal(n, bb, c%v, 1)
  end function VectorDivideC

  type(CVec) function ComplexConjugate(a) result(b)
    class(CVec), intent(in) :: a
    integer(int32) :: n
    n = a%n_size
    if(n < 1) return
    call b%ini(n)
    b%v = conjg(a%v)
  end function ComplexConjugate

  real(real64) function InnerProductC(a, b) result(c)
    type(CVec), intent(in) :: a, b
    integer(int32) :: n
    real(real64) :: zdotu
    c = 0.d0
    if(a%n_size /= b%n_size) then
      write(*,'(a)') 'Error in InnerProduct'
      stop
    end if
    n = a%n_size
    if(n < 1) return
    c = zdotu(n, a%v, 1, b%v, 1)
  end function InnerProductC

  real(real64) function Nrm(a) result(b)
    class(CVec), intent(in) :: a
    integer(int32) :: n
    real(real64) :: dznrm2
    b = 0.d0
    n = a%n_size
    if(n < 1) return
    b = dznrm2(n, a%v, 1)
  end function Nrm

  real(real64) function Nrm2(a) result(b)
    class(CVec), intent(in) :: a
    integer(int32) :: n
    real(real64) :: dznrm2
    n = a%n_size
    b = 0.d0
    if(n < 1) return
    b = dznrm2(n, a%v, 1) ** 2
  end function Nrm2

  subroutine VectorPrint(this, msg, iunit, binary)
    class(CVec),intent(in)::this
    integer, intent(in), optional :: iunit
    character(*), intent(in), optional :: msg
    logical, intent(in), optional :: binary
    logical :: bin
    integer(int32) :: n, i, unt
    if(present(iunit)) then; unt = iunit
    else; unt = 6; end if

    if(present(binary)) then; bin = binary
    else; bin = .false.; end if

    if(bin) then
      write(unt) this%v
    else
      if(present(msg)) write(unt,*) msg
      if(unt == 6) then
        write(unt,'(a)',advance='no') 'Real:'
        write(unt,'(10f10.4)') dble(this%v(:))
        write(unt,'(a)',advance='no') 'Imag:'
        write(unt,'(10f10.4)') aimag(this%v(:))
      else
        n = size(this%v, 1)
        do i = 1, n
          write(unt,'(10f10.4)') this%v(i)
        end do
      end if
    end if

  end subroutine VectorPrint

  subroutine GetRandomVector(v, n, dist)
    ! idist = = 1:  real and imaginary parts each uniform (0,1)
    ! idist = = 2:  real and imaginary parts each uniform (-1,1)
    ! idist = = 3:  real and imaginary parts each normal (0,1)
    ! idist = = 4:  uniformly distributed on the disc abs(z) < 1
    ! idist = = 5:  uniformly distributed on the circle abs(z) = 1
    class(CVec), intent(inout) :: v
    integer(int32), intent(in) :: n
    integer(int32), intent(in), optional :: dist
    integer(int32) :: idist = 4
    if(n < 1) return
    if(present(dist)) idist = dist
    call v%ini(n)
    call zlarnv(idist, iseed, n, v%v)
  end subroutine GetRandomVector

  function block_cvec(this, n1, n2) result(r)
    class(CVec), intent(in) :: this
    type(CVec) :: r
    integer(int32), intent(in) :: n1, n2
    integer(int32) :: n
    n = n2 - n1 + 1
    if(n < 1) return
    call r%ini(n)
    r%v(:) = this%v(n1:n2)
  end function block_cvec
end module vector_complex
