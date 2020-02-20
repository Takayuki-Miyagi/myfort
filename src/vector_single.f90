module vector_single
  use basic_types
  implicit none

  public :: SVec
  public :: VectorCopyS
  public :: VectorSumS
  public :: VectorSubtractS
  public :: VectorScaleRS
  public :: VectorScaleLS
  public :: VectorDivideS
  public :: InnerProductS

  private :: IniV
  private :: zeros
  private :: FinV
  private :: VectorPrint
  private :: GetRandomVector
  private :: Nrm
  private :: Nrm2
  private :: block_svec

  type :: SVec
    real(real32), allocatable :: V(:)
    integer(int32) :: n_size = 0
  contains
    procedure :: Ini => iniV
    procedure :: zeros
    procedure :: Fin => FinV
    procedure :: prt => VectorPrint
    procedure :: Random => GetRandomVector
    procedure :: blk => block_svec
    procedure :: Nrm
    procedure :: Nrm2
  end type SVec
contains
  subroutine IniV(a, n)
    class(SVec), intent(inout) :: a
    integer(int32), intent(in) :: n
    if(allocated(a%V)) deallocate(a%V)
    a%n_size = n
    allocate(a%V(a%n_size))
  end subroutine IniV

  subroutine zeros(a, n)
    class(SVec), intent(inout) :: a
    integer(int32), intent(in) :: n
    call a%ini(n)
    a%v(:) = 0.0
  end subroutine zeros

  subroutine FinV(a)
    class(SVec), intent(inout) :: a
    if(allocated(a%V)) deallocate(a%V)
    a%n_size = 0
  end subroutine FinV

  subroutine VectorCopyS(b, a)
    type(SVec), intent(inout) :: b
    type(SVec), intent(in) :: a
    integer(int32) :: n
    n = a%n_size
    if(n < 1) return
    call b%Ini(n)
    call scopy(n, a%v, 1, b%v, 1)
  end subroutine VectorCopyS

  type(SVec) function VectorSumS(a, b) result(c)
    type(SVec), intent(in) :: a, b
    integer(int32) :: n
    if(a%n_size /= b%n_size) then
      write(*,'(a)') 'Error in DVectorSum'
      stop
    end if
    n = a%n_size
    if(n < 1) return
    call VectorCopyS(c, a)
    call saxpy(n, 1.0, b%v, 1, c%v, 1)
  end function VectorSumS

  type(SVec) function VectorSubtractS(a, b) result(c)
    type(SVec), intent(in) :: a, b
    integer(int32) :: n
    if(a%n_size /= b%n_size) then
      write(*,'(a)') 'Error in SVectorSubtract'
      stop
    end if
    n = a%n_size
    if(n < 1) return
    call VectorCopyS(c, a)
    call saxpy(n, -1.0, b%v, 1, c%v, 1)
  end function VectorSubtractS

  type(SVec) function VectorScaleRS(a, b) result(c)
    type(SVec), intent(in) :: a
    real(real32), intent(in) :: b
    integer(int32) :: n
    n = a%n_size
    if(n < 1) return
    call VectorCopyS(c, a)
    call sscal(n, b, c%v, 1)
  end function VectorScaleRS

  type(SVec) function VectorScaleLS(b, a) result(c)
    type(SVec), intent(in) :: a
    real(real32), intent(in) :: b
    integer(int32) :: n
    n = a%n_size
    if(n < 1) return
    call VectorCopyS(c, a)
    call sscal(n, b, c%v, 1)
  end function VectorScaleLS

  type(SVec) function VectorDivideS(a, b) result(c)
    type(SVec), intent(in) :: a
    real(real32), intent(in) :: b
    integer(int32) :: n
    n = a%n_size
    if(n < 1) return
    call VectorCopyS(c, a)
    call sscal(n, 1.0 / b, c%v, 1)
  end function VectorDivideS

  real(real32) function InnerProductS(a, b) result(c)
    type(SVec), intent(in) :: a, b
    integer(int32) :: n
    if(a%n_size /= b%n_size) then
      write(*,'(a)') 'Error in InnerProduct'
      stop
    end if
    c = 0.0
    n = a%n_size
    if(n < 1) return
    c = dot_product(a%v, b%v)
  end function InnerProductS

  real(real32) function Nrm(a) result(b)
    class(SVec), intent(in) :: a
    integer(int32) :: n
    b = 0.0
    n = a%n_size
    if(n < 1) return
    b = sqrt( dot_product(a%v, a%v) )
  end function Nrm

  real(real32) function Nrm2(a) result(b)
    class(SVec), intent(in) :: a
    integer(int32) :: n
    b = 0.0
    n = a%n_size
    if(n < 1) return
    b = dot_product(a%v, a%v)
  end function Nrm2

  subroutine VectorPrint(this, msg, iunit, binary)
    class(SVec),intent(in)::this
    integer(int32) :: i, n, unt
    integer(int32), intent(in), optional :: iunit
    character(*), intent(in), optional :: msg
    logical, intent(in), optional :: binary
    logical :: bin

    if(present(iunit)) then; unt = iunit
    else; unt = 6; end if

    if(present(binary)) then; bin = binary
    else; bin = .false.; end if

    if(bin) then
      write(unt) this%v
    else
      if(present(msg)) write(unt,*) msg
      if(unt == 6) then
        write(unt,'(10f10.4)') this%v(:)
      else
        n = size(this%v, 1)
        do i = 1, n
          write(unt,'(10f10.4)') this%v(i)
        end do
      end if
    end if
  end subroutine VectorPrint

  subroutine GetRandomVector(v, n, dist)
    ! idist = 1: uniform (0, 1)
    ! idist = 2: uniform (-1, 1)
    ! idist = 3: normal (-1, 1)
    class(SVec), intent(inout) :: v
    integer(int32), intent(in) :: n
    integer(int32), intent(in), optional :: dist
    integer(int32) :: idist = 3
    if(n < 1) return
    if(present(dist)) idist = dist
    call v%ini(n)
    call slarnv(idist, iseed, n, v%v)
  end subroutine GetRandomVector

  function block_SVec(this, n1, n2) result(r)
    class(SVec), intent(in) :: this
    type(SVec) :: r
    integer(int32), intent(in) :: n1, n2
    integer(int32) :: n
    n = n2 - n1 + 1
    if(n < 1) return
    call r%ini(n)
    r%v(:) = this%v(n1:n2)
  end function block_SVec
end module vector_single
