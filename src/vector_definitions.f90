module vector_definitions
  use basic_types
  implicit none

  private :: initialize_svec
  private :: zeros_svec
  private :: finalize_svec
  private :: print_svec
  private :: set_random_svec
  private :: get_size_svec
  private :: initialize_dvec
  private :: zeros_dvec
  private :: finalize_dvec
  private :: print_dvec
  private :: set_random_dvec
  private :: get_size_dvec
  private :: initialize_cvec
  private :: zeros_cvec
  private :: finalize_cvec
  private :: print_cvec
  private :: set_random_cvec
  private :: get_size_cvec
  private :: cvec_c_conjugate
  type :: SVec
    real(real32), allocatable :: v(:)
    integer(int32) :: n_size = 0
  contains
    procedure :: initialize_svec
    procedure :: zeros_svec
    procedure :: finalize_svec
    procedure :: print_svec
    procedure :: set_random_svec
    procedure :: get_size_svec
    generic :: ini => initialize_svec
    generic :: fin => finalize_svec
    generic :: zeros => zeros_svec
    generic :: prnt => print_svec
    generic :: rndm => set_random_svec
    generic :: GetSize => get_size_svec
  end type SVec

  type :: DVec
    real(real64), allocatable :: v(:)
    integer(int32) :: n_size = 0
  contains
    procedure :: initialize_dvec
    procedure :: zeros_dvec
    procedure :: finalize_dvec
    procedure :: print_dvec
    procedure :: set_random_dvec
    procedure :: get_size_dvec
    generic :: ini => initialize_dvec
    generic :: fin => finalize_dvec
    generic :: zeros => zeros_dvec
    generic :: prnt => print_dvec
    generic :: rndm => set_random_dvec
    generic :: GetSize => get_size_dvec
  end type DVec

  type :: CVec
    complex(real64), allocatable :: v(:)
    integer(int32) :: n_size = 0
  contains
    procedure :: initialize_cvec
    procedure :: zeros_cvec
    procedure :: finalize_cvec
    procedure :: print_cvec
    procedure :: set_random_cvec
    procedure :: get_size_cvec
    procedure :: cvec_c_conjugate
    generic :: ini => initialize_cvec
    generic :: fin => finalize_cvec
    generic :: zeros => zeros_cvec
    generic :: prnt => print_cvec
    generic :: rndm => set_random_cvec
    generic :: GetSize => get_size_cvec
    generic :: star => cvec_c_conjugate
  end type CVec
contains
  subroutine initialize_svec( this, n )
    class(svec), intent(inout) :: this
    integer(int32), intent(in) :: n
    if(allocated(this%v)) call this%fin()
    this%n_size = n
    allocate(this%v( this%n_size ))
  end subroutine initialize_svec

  subroutine finalize_svec( this )
    class(svec), intent(inout) :: this
    this%n_size = 0
    deallocate( this%v )
  end subroutine finalize_svec

  subroutine zeros_svec( this, n )
    class(svec), intent(inout) :: this
    integer, intent(in) :: n
    call this%ini( n )
    this%v(:) = 0.0
  end subroutine zeros_svec

  function get_size_svec( this ) result(n)
    class(SVec), intent(in) :: this
    integer(int32) :: n
    n = this%n_size
  end function get_size_svec

  subroutine print_svec(this, msg, iunit, binary)
    class(SVec),intent(in)::this
    integer(int32), intent(in), optional :: iunit
    character(*), intent(in), optional :: msg
    logical, intent(in), optional :: binary
    integer(int32) :: i, n, unt
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
  end subroutine print_svec

  subroutine set_random_svec(this, n, dist)
    ! idist = 1: uniform (0, 1)
    ! idist = 2: uniform (-1, 1)
    ! idist = 3: normal (-1, 1)
    class(SVec), intent(inout) :: this
    integer(int32), intent(in) :: n
    integer(int32), intent(in), optional :: dist
    integer(int32) :: idist = 3
    if(n < 1) return
    if(present(dist)) idist = dist
    call this%ini(n)
    call slarnv(idist, iseed, n, this%v)
  end subroutine set_random_svec

  subroutine set_svec_from_array( this, v )
    type(SVec), intent(inout) :: this
    real(real32), intent(in) :: v(:)
    call this%ini( size(v) )
    this%v(:) = v(:)
  end subroutine set_svec_from_array

  subroutine copy_svec( a, b )
    type(SVec), intent(inout) :: a
    class(*), intent(in) :: b
    select type( b )
    type is ( SVec )
      call a%ini( b%n_size )
      a%v(:) = b%v(:)
    type is ( DVec )
      call a%ini( b%n_size )
      a%v(:) = real(b%v(:), kind=real32)
    type is ( CVec )
      call a%ini( b%n_size )
      a%v(:) = real(b%v(:), kind=real32)
    end select
  end subroutine copy_svec

  function sum_svec( a, b ) result(c)
    type(SVec), intent(in) :: a, b
    type(SVec) :: c
    integer(int32) :: n
    if(a%n_size /= b%n_size) then
      write(*,'(a)') 'Error in sum_svec'
      stop
    end if
    n = a%n_size
    if(n < 1) return
    call copy_svec(c, a)
    call saxpy(n, 1.0, b%v, 1, c%v, 1)
  end function sum_svec

  function subtract_svec( a, b ) result(c)
    type(SVec), intent(in) :: a, b
    type(SVec) :: c
    integer(int32) :: n
    if(a%n_size /= b%n_size) then
      write(*,'(a)') 'Error in subtract_svec'
      stop
    end if
    n = a%n_size
    if(n < 1) return
    call copy_svec(c, a)
    call saxpy(n, -1.0, b%v, 1, c%v, 1)
  end function subtract_svec

  function scale_l_svec( b, a ) result(c)
    type(SVec), intent(in) :: a
    real(real32), intent(in) :: b
    type(SVec) :: c
    integer(int32) :: n
    n = a%n_size
    if(n < 1) return
    call copy_svec(c, a)
    call sscal(n, b, c%v, 1)
  end function scale_l_svec

  function scale_r_svec( a, b ) result(c)
    type(SVec), intent(in) :: a
    real(real32), intent(in) :: b
    type(SVec) :: c
    integer(int32) :: n
    n = a%n_size
    if(n < 1) return
    call copy_svec(c, a)
    call sscal(n, b, c%v, 1)
  end function scale_r_svec

  function scale_inv_svec( a, b ) result(c)
    type(SVec), intent(in) :: a
    real(real32), intent(in) :: b
    type(SVec) :: c
    integer(int32) :: n
    n = a%n_size
    if(n < 1) return
    call copy_svec(c, a)
    call sscal(n, 1.0/b, c%v, 1)
  end function scale_inv_svec

  function inner_product_svec(a, b) result(c)
    type(SVec), intent(in) :: a, b
    real(real32) :: c
    integer(int32) :: n
    c = 0.0
    if(a%n_size /= b%n_size) then
      write(*,'(a)') 'error in inner_product_svec'
      return
    end if
    n = a%n_size
    if(n < 1) return
    c = dot_product(a%v, b%v)
  end function inner_product_svec

  subroutine initialize_dvec( this, n )
    class(Dvec), intent(inout) :: this
    integer(int32), intent(in) :: n
    if(allocated(this%v)) call this%fin()
    this%n_size = n
    allocate(this%v( this%n_size ))
  end subroutine initialize_dvec

  subroutine finalize_dvec( this )
    class(Dvec), intent(inout) :: this
    this%n_size = 0
    deallocate( this%v )
  end subroutine finalize_dvec

  subroutine zeros_dvec( this, n )
    class(Dvec), intent(inout) :: this
    integer, intent(in) :: n
    call this%ini( n )
    this%v(:) = 0.d0
  end subroutine zeros_dvec

  function get_size_dvec( this ) result(n)
    class(Dvec), intent(in) :: this
    integer(int32) :: n
    n = this%n_size
  end function get_size_dvec

  subroutine print_dvec(this, msg, iunit, binary)
    class(Dvec),intent(in)::this
    integer(int32), intent(in), optional :: iunit
    character(*), intent(in), optional :: msg
    logical, intent(in), optional :: binary
    integer(int32) :: i, n, unt
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
  end subroutine print_dvec

  subroutine set_random_dvec(this, n, dist)
    ! idist = 1: uniform (0, 1)
    ! idist = 2: uniform (-1, 1)
    ! idist = 3: normal (-1, 1)
    class(Dvec), intent(inout) :: this
    integer(int32), intent(in) :: n
    integer(int32), intent(in), optional :: dist
    integer(int32) :: idist = 3
    if(n < 1) return
    if(present(dist)) idist = dist
    call this%ini(n)
    call dlarnv(idist, iseed, n, this%v)
  end subroutine set_random_dvec

  subroutine set_dvec_from_array( this, v )
    type(Dvec), intent(inout) :: this
    real(real64), intent(in) :: v(:)
    call this%ini( size(v) )
    this%v(:) = v(:)
  end subroutine set_dvec_from_array

  subroutine copy_dvec( a, b )
    type(Dvec), intent(inout) :: a
    class(*), intent(in) :: b
    select type( b )
    type is ( Svec )
      call a%ini( b%n_size )
      a%v(:) = dble(b%v(:))
    type is ( DVec )
      call a%ini( b%n_size )
      a%v(:) = b%v(:)
    type is ( CVec )
      call a%ini( b%n_size )
      a%v(:) = dble(b%v(:))
    end select
  end subroutine copy_dvec

  function sum_dvec( a, b ) result(c)
    type(DVec), intent(in) :: a, b
    type(DVec) :: c
    integer(int32) :: n
    if(a%n_size /= b%n_size) then
      write(*,'(a)') 'Error in sum_dvec'
      stop
    end if
    n = a%n_size
    if(n < 1) return
    call copy_dvec(c, a)
    call daxpy(n, 1.d0, b%v, 1, c%v, 1)
  end function sum_dvec

  function subtract_dvec( a, b ) result(c)
    type(DVec), intent(in) :: a, b
    type(DVec) :: c
    integer(int32) :: n
    if(a%n_size /= b%n_size) then
      write(*,'(a)') 'Error in subtract_dvec'
      stop
    end if
    n = a%n_size
    if(n < 1) return
    call copy_dvec(c, a)
    call daxpy(n, -1.d0, b%v, 1, c%v, 1)
  end function subtract_dvec

  function scale_l_dvec( b, a ) result(c)
    type(DVec), intent(in) :: a
    real(real64), intent(in) :: b
    type(DVec) :: c
    integer(int32) :: n
    n = a%n_size
    if(n < 1) return
    call copy_dvec(c, a)
    call dscal(n, b, c%v, 1)
  end function scale_l_dvec

  function scale_r_dvec( a, b ) result(c)
    type(DVec), intent(in) :: a
    real(real64), intent(in) :: b
    type(DVec) :: c
    integer(int32) :: n
    n = a%n_size
    if(n < 1) return
    call copy_dvec(c, a)
    call dscal(n, b, c%v, 1)
  end function scale_r_dvec

  function scale_inv_dvec( a, b ) result(c)
    type(DVec), intent(in) :: a
    real(real64), intent(in) :: b
    type(DVec) :: c
    integer(int32) :: n
    n = a%n_size
    if(n < 1) return
    call copy_dvec(c, a)
    call dscal(n, 1.d0/b, c%v, 1)
  end function scale_inv_dvec

  function inner_product_dvec(a, b) result(c)
    type(dvec), intent(in) :: a, b
    real(real64) :: c
    real(real64) :: ddot
    integer(int32) :: n
    c = 0.d0
    if(a%n_size /= b%n_size) then
      write(*,'(a)') 'error in inner_product_dvec'
      return
    end if
    n = a%n_size
    if(n < 1) return
    c = ddot(n, a%v, 1, b%v, 1)
  end function inner_product_dvec

  subroutine initialize_cvec( this, n )
    class(CVec), intent(inout) :: this
    integer(int32), intent(in) :: n
    if(allocated(this%v)) call this%fin()
    this%n_size = n
    allocate(this%v( this%n_size ))
  end subroutine initialize_cvec

  subroutine finalize_cvec( this )
    class(CVec), intent(inout) :: this
    this%n_size = 0
    deallocate( this%v )
  end subroutine finalize_cvec

  subroutine zeros_cvec( this, n )
    class(CVec), intent(inout) :: this
    integer, intent(in) :: n
    call this%ini( n )
    this%v(:) = (0.d0,0.d0)
  end subroutine zeros_cvec

  function get_size_cvec( this ) result(n)
    class(CVec), intent(in) :: this
    integer(int32) :: n
    n = this%n_size
  end function get_size_cvec

  subroutine print_cvec(this, msg, iunit, binary)
    class(CVec),intent(in)::this
    integer(int32), intent(in), optional :: iunit
    character(*), intent(in), optional :: msg
    logical, intent(in), optional :: binary
    integer(int32) :: i, n, unt
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
  end subroutine print_cvec

  subroutine set_random_cvec(this, n, dist)
    ! idist = = 1:  real and imaginary parts each uniform (0,1)
    ! idist = = 2:  real and imaginary parts each uniform (-1,1)
    ! idist = = 3:  real and imaginary parts each normal (0,1)
    ! idist = = 4:  uniformly distributed on the disc abs(z) < 1
    ! idist = = 5:  uniformly distributed on the circle abs(z) = 1
    class(CVec), intent(inout) :: this
    integer(int32), intent(in) :: n
    integer(int32), intent(in), optional :: dist
    integer(int32) :: idist = 4
    if(n < 1) return
    if(present(dist)) idist = dist
    call this%ini(n)
    call zlarnv(idist, iseed, n, this%v)
  end subroutine set_random_cvec

  subroutine set_cvec_from_array( this, v )
    type(cvec), intent(inout) :: this
    real(real64), intent(in) :: v(:)
    call this%ini( size(v) )
    this%v(:) = v(:)
  end subroutine set_cvec_from_array

  subroutine copy_cvec( a, b )
    type(CVec), intent(inout) :: a
    class(*), intent(in) :: b
    select type( b )
    type is ( SVec )
      call a%ini( b%n_size )
      a%v(:) = dble(b%v(:))
    type is ( DVec )
      call a%ini( b%n_size )
      a%v(:) = b%v(:)
    type is ( CVec )
      call a%ini( b%n_size )
      a%v(:) = b%v(:)
    end select
  end subroutine copy_cvec

  function sum_cvec( a, b ) result(c)
    type(CVec), intent(in) :: a, b
    type(CVec) :: c
    integer(int32) :: n
    if(a%n_size /= b%n_size) then
      write(*,'(a)') 'Error in sum_cvec'
      stop
    end if
    n = a%n_size
    if(n < 1) return
    call copy_cvec(c, a)
    call zaxpy(n, (1.d0,0.d0), b%v, 1, c%v, 1)
  end function sum_cvec

  function subtract_cvec( a, b ) result(c)
    type(CVec), intent(in) :: a, b
    type(CVec) :: c
    integer(int32) :: n
    if(a%n_size /= b%n_size) then
      write(*,'(a)') 'Error in subtract_cvec'
      stop
    end if
    n = a%n_size
    if(n < 1) return
    call copy_cvec(c, a)
    call zaxpy(n, (-1.d0,0.d0), b%v, 1, c%v, 1)
  end function subtract_cvec

  function scale_l_cvec( b, a ) result(c)
    type(CVec), intent(in) :: a
    complex(real64), intent(in) :: b
    type(cvec) :: c
    integer(int32) :: n
    n = a%n_size
    if(n < 1) return
    call copy_cvec(c, a)
    call zscal(n, b, c%v, 1)
  end function scale_l_cvec

  function scale_r_cvec( a, b ) result(c)
    type(CVec), intent(in) :: a
    complex(real64), intent(in) :: b
    type(CVec) :: c
    integer(int32) :: n
    n = a%n_size
    if(n < 1) return
    call copy_cvec(c, a)
    call zscal(n, b, c%v, 1)
  end function scale_r_cvec

  function scale_inv_cvec( a, b ) result(c)
    type(CVec), intent(in) :: a
    complex(real64), intent(in) :: b
    type(CVec) :: c
    integer(int32) :: n
    n = a%n_size
    if(n < 1) return
    call copy_cvec(c, a)
    call zscal(n, (1.d0,0.d0)/b, c%v, 1)
  end function scale_inv_cvec

  function cvec_c_conjugate(a) result(b)
    class(CVec), intent(in) :: a
    type(CVec) :: b
    integer(int32) :: n
    n = a%n_size
    if(n < 1) return
    call b%ini(n)
    b%v = conjg(a%v)
  end function cvec_c_conjugate

  function inner_product_cvec(a, b) result(c)
    type(cvec), intent(in) :: a, b
    complex(real64) :: c
    complex(real64) :: zdotu
    integer(int32) :: n
    c = (0.d0,0.d0)
    if(a%n_size /= b%n_size) then
      write(*,'(a)') 'error in inner_product_cvec'
      return
    end if
    n = a%n_size
    if(n < 1) return
    c = zdotu(n, a%v, 1, b%v, 1)
  end function inner_product_cvec
end module vector_definitions
