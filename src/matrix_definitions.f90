module matrix_definitions
  use basic_types
  use vector_definitions
  implicit none

  public :: SVec
  public :: DVec
  public :: CVec
  public :: SMat
  public :: DMat
  public :: CMat
  public :: assignment(=)
  public :: operator(+)
  public :: operator(-)
  public :: operator(*)
  public :: operator(/)
  public :: operator(.x.)

  private :: initialize_smat
  private :: zeros_smat
  private :: finalize_smat
  private :: eye_smat
  private :: print_smat
  private :: set_random_smat
  private :: get_size_smat
  private :: trans_smat
  private :: inverse_smat
  private :: determinant_smat
  private :: initialize_dmat
  private :: zeros_dmat
  private :: finalize_dmat
  private :: eye_dmat
  private :: print_dmat
  private :: set_random_dmat
  private :: get_size_dmat
  private :: trans_dmat
  private :: inverse_dmat
  private :: determinant_dmat
  private :: initialize_cmat
  private :: zeros_cmat
  private :: finalize_cmat
  private :: eye_cmat
  private :: print_cmat
  private :: set_random_cmat
  private :: get_size_cmat
  private :: trans_cmat
  private :: c_conjugate_cmat
  private :: h_conjugate_cmat
  private :: inverse_cmat
  private :: determinant_cmat
  private :: copy_svec
  private :: copy_dvec
  private :: copy_cvec
  private :: set_svec_from_array
  private :: set_dvec_from_array
  private :: set_cvec_from_array
  private :: copy_smat
  private :: copy_dmat
  private :: copy_cmat
  private :: set_smat_from_array
  private :: set_dmat_from_array
  private :: set_cmat_from_array
  private :: sum_svec
  private :: sum_dvec
  private :: sum_cvec
  private :: sum_smat
  private :: sum_dmat
  private :: sum_cmat
  private :: subtract_svec
  private :: subtract_dvec
  private :: subtract_cvec
  private :: subtract_smat
  private :: subtract_dmat
  private :: subtract_cmat
  private :: scale_l_svec
  private :: scale_r_svec
  private :: inner_product_svec
  private :: product_smat
  private :: scale_l_smat
  private :: scale_r_smat
  private :: smat_svec_product
  private :: svec_smat_product
  private :: scale_l_dvec
  private :: scale_r_dvec
  private :: inner_product_dvec
  private :: product_dmat
  private :: scale_l_dmat
  private :: scale_r_dmat
  private :: dmat_dvec_product
  private :: dvec_dmat_product
  private :: scale_l_cvec
  private :: scale_r_cvec
  private :: inner_product_cvec
  private :: product_cmat
  private :: scale_l_cmat
  private :: scale_r_cmat
  private :: cmat_cvec_product
  private :: cvec_cmat_product
  private :: outer_product_smat
  private :: outer_product_dmat
  private :: outer_product_cmat

  type :: SMat
    real(real32), allocatable :: m(:,:)
    integer(int32) :: n_row = 0
    integer(int32) :: n_col = 0
  contains
    procedure :: initialize_smat
    procedure :: zeros_smat
    procedure :: finalize_smat
    procedure :: eye_smat
    procedure :: print_smat
    procedure :: set_random_smat
    procedure :: get_size_smat
    procedure :: trans_smat
    procedure :: inverse_smat
    procedure :: determinant_smat
    procedure :: set_diag_smat
    generic :: ini => initialize_smat
    generic :: fin => finalize_smat
    generic :: zeros => zeros_smat
    generic :: eye => eye_smat
    generic :: prnt => print_smat
    generic :: rndm => set_random_smat
    generic :: t => trans_smat
    generic :: inv => inverse_smat
    generic :: det => determinant_smat
    generic :: GetSize => get_size_smat
    generic :: smat => set_diag_smat
  end type SMat

  type :: DMat
    real(real64), allocatable :: m(:,:)
    integer(int32) :: n_row = 0
    integer(int32) :: n_col = 0
  contains
    procedure :: initialize_dmat
    procedure :: zeros_dmat
    procedure :: finalize_dmat
    procedure :: eye_dmat
    procedure :: print_dmat
    procedure :: set_random_dmat
    procedure :: get_size_dmat
    procedure :: trans_dmat
    procedure :: inverse_dmat
    procedure :: determinant_dmat
    procedure :: set_diag_dmat
    generic :: ini => initialize_dmat
    generic :: fin => finalize_dmat
    generic :: zeros => zeros_dmat
    generic :: eye => eye_dmat
    generic :: prnt => print_dmat
    generic :: rndm => set_random_dmat
    generic :: GetSize => get_size_dmat
    generic :: t => trans_dmat
    generic :: inv => inverse_dmat
    generic :: det => determinant_dmat
    generic :: diag => set_diag_dmat
  end type DMat

  type :: CMat
    complex(real64), allocatable :: m(:,:)
    integer(int32) :: n_row = 0
    integer(int32) :: n_col = 0
  contains
    procedure :: initialize_cmat
    procedure :: zeros_cmat
    procedure :: finalize_cmat
    procedure :: eye_cmat
    procedure :: print_cmat
    procedure :: set_random_cmat
    procedure :: get_size_cmat
    procedure :: trans_cmat
    procedure :: c_conjugate_cmat
    procedure :: h_conjugate_cmat
    procedure :: inverse_cmat
    procedure :: determinant_cmat
    procedure :: set_diag_cmat
    generic :: ini => initialize_cmat
    generic :: fin => finalize_cmat
    generic :: zeros => zeros_cmat
    generic :: eye => eye_cmat
    generic :: prnt => print_cmat
    generic :: rndm => set_random_cmat
    generic :: GetSize => get_size_cmat
    generic :: t => trans_cmat
    generic :: star => c_conjugate_cmat
    generic :: h => h_conjugate_cmat
    generic :: inv => inverse_cmat
    generic :: det => determinant_cmat
    generic :: diag => set_diag_cmat
  end type CMat

  interface assignment(=)
    module procedure :: copy_svec
    module procedure :: copy_dvec
    module procedure :: copy_cvec
    module procedure :: set_svec_from_array
    module procedure :: set_dvec_from_array
    module procedure :: set_cvec_from_array
    module procedure :: copy_smat
    module procedure :: copy_dmat
    module procedure :: copy_cmat
    module procedure :: set_smat_from_array
    module procedure :: set_dmat_from_array
    module procedure :: set_cmat_from_array
  end interface assignment(=)

  interface operator(+)
    module procedure :: sum_svec
    module procedure :: sum_dvec
    module procedure :: sum_cvec
    module procedure :: sum_smat
    module procedure :: sum_dmat
    module procedure :: sum_cmat
  end interface operator(+)

  interface operator(-)
    module procedure :: subtract_svec
    module procedure :: subtract_dvec
    module procedure :: subtract_cvec
    module procedure :: subtract_smat
    module procedure :: subtract_dmat
    module procedure :: subtract_cmat
  end interface operator(-)

  interface operator(*)
    module procedure :: scale_l_svec
    module procedure :: scale_r_svec
    module procedure :: inner_product_svec
    module procedure :: product_smat
    module procedure :: scale_l_smat
    module procedure :: scale_r_smat
    module procedure :: smat_svec_product
    module procedure :: svec_smat_product
    module procedure :: scale_l_dvec
    module procedure :: scale_r_dvec
    module procedure :: inner_product_dvec
    module procedure :: product_dmat
    module procedure :: scale_l_dmat
    module procedure :: scale_r_dmat
    module procedure :: dmat_dvec_product
    module procedure :: dvec_dmat_product
    module procedure :: scale_l_cvec
    module procedure :: scale_r_cvec
    module procedure :: inner_product_cvec
    module procedure :: product_cmat
    module procedure :: scale_l_cmat
    module procedure :: scale_r_cmat
    module procedure :: cmat_cvec_product
    module procedure :: cvec_cmat_product
  end interface operator(*)

  interface operator(.x.)
    module procedure :: outer_product_smat
    module procedure :: outer_product_dmat
    module procedure :: outer_product_cmat
  end interface operator(.x.)

  interface operator(/)
    module procedure :: scale_inv_svec
    module procedure :: scale_inv_dvec
    module procedure :: scale_inv_cvec
    module procedure :: scale_inv_smat
    module procedure :: scale_inv_dmat
    module procedure :: scale_inv_cmat
  end interface operator(/)

contains
  subroutine initialize_smat( this, n_row, n_col )
    class(SMat), intent(inout) :: this
    integer(int32), intent(in) :: n_row, n_col
    if( allocated( this%m )) call this%fin()
    if( n_row < 1 .or. n_col < 1 ) return
    this%n_row = n_row
    this%n_col = n_col
    allocate( this%m( this%n_row, this%n_col ))
  end subroutine initialize_smat

  subroutine finalize_smat( this )
    class(SMat), intent(inout) :: this
    deallocate( this%m )
    this%n_row = 0
    this%n_col = 0
  end subroutine finalize_smat

  subroutine zeros_smat( this, n_row, n_col )
    class(SMat), intent(inout) :: this
    integer(int32), intent(in) :: n_row, n_col
    if(allocated( this%m )) call this%fin()
    if(n_row < 1 .or. n_col < 1) return
    call this%ini( n_row, n_col )
    this%m(:,:) = 0.0
  end subroutine zeros_smat

  subroutine eye_smat( this, n )
    class(SMat), intent(inout) :: this
    integer(int32), intent(in) :: n
    integer(int32) :: i
    if(allocated( this%m )) call this%fin()
    if(n < 1) return
    call this%zeros( n, n )
    do i = 1, n
      this%m(i,i) = 1.0
    end do
  end subroutine eye_smat

  function get_size_smat( this ) result( ns )
    class(SMat), intent(in) :: this
    integer(int32) :: ns(2)
    ns(1) = this%n_row
    ns(2) = this%n_col
  end function get_size_smat

  subroutine copy_smat( a, b )
    type(SMat), intent(inout) :: a
    class(*), intent(in) :: b
    select type( b )
    type is ( SMat )
      call a%ini( b%n_row, b%n_col )
      a%m(:,:) = b%m(:,:)
    type is ( DMat )
      call a%ini( b%n_row, b%n_col )
      a%m(:,:) = real(b%m(:,:), kind=real32)
    type is ( CMat )
      call a%ini( b%n_row, b%n_col )
      a%m(:,:) = real(b%m(:,:), kind=real32)
    end select
  end subroutine copy_smat

  function product_smat( a, b ) result(c)
    type(SMat), intent(in) :: a, b
    type(SMat) :: c
    integer(int32) :: m, k, n
    m = a%n_row
    k = a%n_col
    if(a%n_col /= b%n_row) then
      write(*, '(a)') 'Error in product_smat'
      stop
    end if
    n = b%n_col
    call c%ini(m,n)
    call sgemm('n','n',m,n,k,1.0,a%m,m,b%m,k,0.0,c%m,m)
  end function product_smat

  function sum_smat( a, b ) result(c)
    type(SMat), intent(in) :: a, b
    type(SMat) :: c
    integer(int32) :: m, n, i
    if(a%n_row /= b%n_row .or. a%n_col /= b%n_col) then
      write(*, '(a)') 'Error in sum_smat'
      stop
    end if
    m = a%n_row
    n = a%n_col
    if(m < 1 .or. n < 1) return
    call copy_smat(c, a)
    do i = 1, n
      call saxpy(m, 1.0, b%m(:,i), 1, c%m(:,i), 1)
    end do
  end function sum_smat

  function subtract_smat( a, b ) result(c)
    type(SMat), intent(in) :: a, b
    type(SMat) :: c
    integer(int32) :: m, n, i
    if(a%n_row /= b%n_row .or. a%n_col /= b%n_col) then
      write(*, '(a)') 'Error in subtract_smat'
      stop
    end if
    m = a%n_row
    n = a%n_col
    if(m < 1 .or. n < 1) return
    call copy_smat(c, a)
    do i = 1, n
      call saxpy(m, -1.0, b%m(:,i), 1, c%m(:,i), 1)
    end do
  end function subtract_smat

  function scale_l_smat( b, a ) result(c)
    type(SMat), intent(in) :: a
    real(real32), intent(in) :: b
    type(SMat) :: c
    integer(int32) :: m, n, i
    m = a%n_row
    n = a%n_col
    if(m < 1 .or. n < 1) return
    call copy_smat(c, a)
    do i = 1, n
      call sscal(m, b, c%m(:,i), 1)
    end do
  end function scale_l_smat

  function scale_r_smat( a, b ) result(c)
    type(SMat), intent(in) :: a
    real(real32), intent(in) :: b
    type(SMat) :: c
    integer(int32) :: m, n, i
    m = a%n_row
    n = a%n_col
    if(m < 1 .or. n < 1) return
    call copy_smat(c, a)
    do i = 1, n
      call sscal(m, b, c%m(:,i), 1)
    end do
  end function scale_r_smat

  function scale_inv_smat( a, b ) result(c)
    type(SMat), intent(in) :: a
    real(real32), intent(in) :: b
    type(SMat) :: c
    integer(int32) :: m, n, i
    m = a%n_row
    n = a%n_col
    if(m < 1 .or. n < 1) return
    call copy_smat(c, a)
    do i = 1, n
      call sscal(m, 1.0/b, c%m(:,i), 1)
    end do
  end function scale_inv_smat

  function trans_smat( a ) result(b)
    class(SMat), intent(in) :: a
    type(SMat) :: b
    integer(int32) :: n, m
    m = a%n_row
    n = a%n_col
    if(m < 1 .or. n < 1) return
    call b%ini(n,m)
    b%m = transpose(a%m)
  end function trans_smat

  function inverse_smat(r) result(s)
    class(SMat), intent(in) :: r
    type(SMat) :: s
    real(real32), allocatable :: a(:,:)
    real(real32), allocatable :: work(:)
    integer(int32), allocatable :: ipvt(:)
    integer(int32) :: info, n
    if( r%n_row /= r%n_col ) then
      write(*,*) "Error in inverse_smat"
      return
    end if
    n = r%n_row
    if(n < 1) return
    call s%ini(n,n)
    allocate(work(n*n),ipvt(n))
    allocate(a(n,n))
    a = r%m
    call sgetrf(n,n,a,n,ipvt,info)
    call sgetri(n,a,n,ipvt,work,n**2,info)
    s%m = a
    deallocate(a,work,ipvt)
  end function inverse_smat

  function determinant_smat(r) result(d)
    class(SMat), intent(in) :: r
    real(real32) :: d
    integer(int32) :: n, i, info
    real(real32), allocatable :: a(:,:)
    integer(int32), allocatable :: ipiv(:)
    d = 0.0
    if( r%n_row /= r%n_col ) then
      write(*,*) "Error in determinant_smat"
      return
    end if
    n = r%n_row
    if(n < 1) return
    allocate(ipiv(n), a(n,n))
    a = r%m
    call sgetrf(n, n, a, n, ipiv, info)
    if(info /= 0) then
      write(*,'(a, i3)') "error in det: info = ", info
      stop
    end if
    d = 1.0
    do i = 1, n
      if(ipiv(i) .ne. i) then
        d = -d * a(i, i)
      else
        d = d * a(i, i)
      end if
    end do
    deallocate(ipiv, a)
  end function determinant_smat

  subroutine set_random_smat(this, n_row, n_col)
    class(SMat), intent(inout) :: this
    integer(int32), intent(in) :: n_row, n_col
    integer(int32) :: i
    type(SVec) :: v
    call this%ini(n_row,n_col)
    if(n_row < 1 .or. n_col < 1) return
    do i = 1, n_col
      call v%rndm(n_row)
      this%m(:,i) = v%v(:)
      call v%fin()
    end do
  end subroutine set_random_smat

  subroutine set_smat_from_array( this, m )
    type(SMat), intent(inout) :: this
    real(real32), intent(in) :: m(:,:)
    call this%ini( size(m,1), size(m,2) )
    this%m(:,:) = m(:,:)
  end subroutine set_smat_from_array

  subroutine set_diag_smat(b, a)
    class(SMat), intent(inout) :: b
    type(SVec), intent(in) :: a
    integer(int32) :: n, i
    n = a%n_size
    if(n < 1) return
    call b%zeros(n,n)
    do i = 1, n
      b%M(i,i) = a%V(i)
    end do
  end subroutine set_diag_smat

  function outer_product_smat(a, b) result(c)
    type(SVec), intent(in) :: a, b
    type(SMat) :: c
    integer(int32) :: n, m
    n = a%n_size
    m = b%n_size
    call c%zeros(n,m)
    call sger(n, m, 1.0, a%v, 1, b%v, 1, c%m, n)
  end function outer_product_smat

  function smat_svec_product(a, b) result(c)
    type(SMat), intent(in) :: a
    type(SVec), intent(in) :: b
    type(SVec) :: c
    integer(int32) :: m, k
    m = a%n_row
    k = a%n_col
    if( k /= b%n_size ) then
      write(*, '(a)') 'Error in smat_svec_product'
      stop
    end if
    call c%ini(m)
    call sgemv('n',m,k,1.0,a%m,m,b%v,1,0.0,c%v,1)
  end function smat_svec_product

  function svec_smat_product(a, b) result(c)
    type(SMat), intent(in) :: b
    type(SVec), intent(in) :: a
    type(SVec) :: c
    integer(int32) :: m, k, n
    m = b%n_row
    k = b%n_col
    n = a%n_size
    if(m /= n) then
      write(*, '(a)') 'Error in svec_smat_product'
      stop
    end if
    call c%ini(k)
    call sgemv('t',m,k,1.0,b%m,m,a%v,1,0.0,c%v,1)
  end function svec_smat_product

  subroutine print_smat(this, msg, iunit, binary)
    class(SMat), intent(in) :: this
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
  end subroutine print_smat

  subroutine initialize_dmat( this, n_row, n_col )
    class(DMat), intent(inout) :: this
    integer(int32), intent(in) :: n_row, n_col
    if( allocated( this%m )) call this%fin()
    if( n_row < 1 .or. n_col < 1 ) return
    this%n_row = n_row
    this%n_col = n_col
    allocate( this%m( this%n_row, this%n_col ))
  end subroutine initialize_dmat

  subroutine finalize_dmat( this )
    class(DMat), intent(inout) :: this
    deallocate( this%m )
    this%n_row = 0
    this%n_col = 0
  end subroutine finalize_dmat

  subroutine zeros_dmat( this, n_row, n_col )
    class(DMat), intent(inout) :: this
    integer(int32), intent(in) :: n_row, n_col
    if(allocated( this%m )) call this%fin()
    if(n_row < 1 .or. n_col < 1) return
    call this%ini( n_row, n_col )
    this%m(:,:) = 0.d0
  end subroutine zeros_dmat

  subroutine eye_dmat( this, n )
    class(DMat), intent(inout) :: this
    integer(int32), intent(in) :: n
    integer(int32) :: i
    if(allocated( this%m )) call this%fin()
    if(n < 1) return
    call this%zeros( n, n )
    do i = 1, n
      this%m(i,i) = 1.d0
    end do
  end subroutine eye_dmat

  function get_size_dmat( this ) result( ns )
    class(DMat), intent(in) :: this
    integer(int32) :: ns(2)
    ns(1) = this%n_row
    ns(2) = this%n_col
  end function get_size_dmat

  subroutine copy_dmat( a, b )
    type(DMat), intent(inout) :: a
    class(*), intent(in) :: b
    select type( b )
    type is ( SMat )
      call a%ini( b%n_row, b%n_col )
      a%m(:,:) = dble(b%m(:,:))
    type is ( DMat )
      call a%ini( b%n_row, b%n_col )
      a%m(:,:) = b%m(:,:)
    type is ( CMat )
      call a%ini( b%n_row, b%n_col )
      a%m(:,:) = dble(b%m(:,:))
    end select
  end subroutine copy_dmat

  function product_dmat( a, b ) result(c)
    type(DMat), intent(in) :: a, b
    type(DMat) :: c
    integer(int32) :: m, k, n
    m = a%n_row
    k = a%n_col
    if(a%n_col /= b%n_row) then
      write(*, '(a)') 'Error in product_dmat'
      stop
    end if
    n = b%n_col
    call c%zeros(m,n)
    call dgemm('n','n',m,n,k,1.d0,a%m,m,b%m,k,0.d0,c%m,m)
  end function product_dmat

  function sum_dmat( a, b ) result(c)
    type(DMat), intent(in) :: a, b
    type(DMat) :: c
    integer(int32) :: m, n, i
    if(a%n_row /= b%n_row .or. a%n_col /= b%n_col) then
      write(*, '(a)') 'Error in sum_dmat'
      stop
    end if
    m = a%n_row
    n = a%n_col
    if(m < 1 .or. n < 1) return
    call copy_dmat(c, a)
    do i = 1, n
      call daxpy(m, 1.d0, b%m(:,i), 1, c%m(:,i), 1)
    end do
  end function sum_dmat

  function subtract_dmat( a, b ) result(c)
    type(DMat), intent(in) :: a, b
    type(DMat) :: c
    integer(int32) :: m, n, i
    if(a%n_row /= b%n_row .or. a%n_col /= b%n_col) then
      write(*, '(a)') 'Error in subtract_dmat'
      stop
    end if
    m = a%n_row
    n = a%n_col
    if(m < 1 .or. n < 1) return
    call copy_dmat(c,a)
    do i = 1, n
      call daxpy(m, -1.d0, b%m(:,i), 1, c%m(:,i), 1)
    end do
  end function subtract_dmat

  function scale_l_dmat( b, a ) result(c)
    type(DMat), intent(in) :: a
    real(real64), intent(in) :: b
    type(DMat) :: c
    integer(int32) :: m, n, i
    m = a%n_row
    n = a%n_col
    if(m < 1 .or. n < 1) return
    call copy_dmat(c, a)
    do i = 1, n
      call dscal(m, b, c%m(:,i), 1)
    end do
  end function scale_l_dmat

  function scale_r_dmat( a, b ) result(c)
    type(DMat), intent(in) :: a
    real(real64), intent(in) :: b
    type(DMat) :: c
    integer(int32) :: m, n, i
    m = a%n_row
    n = a%n_col
    if(m < 1 .or. n < 1) return
    call copy_dmat(c, a)
    do i = 1, n
      call dscal(m, b, c%m(:,i), 1)
    end do
  end function scale_r_dmat

  function scale_inv_dmat( a, b ) result(c)
    type(DMat), intent(in) :: a
    real(real64), intent(in) :: b
    type(DMat) :: c
    integer(int32) :: m, n, i
    m = a%n_row
    n = a%n_col
    if(m < 1 .or. n < 1) return
    call copy_dmat(c, a)
    do i = 1, n
      call dscal(m, 1.d0/b, c%m(:,i), 1)
    end do
  end function scale_inv_dmat

  function trans_dmat( a ) result(b)
    class(DMat), intent(in) :: a
    type(DMat) :: b
    integer(int32) :: n, m
    m = a%n_row
    n = a%n_col
    if(m < 1 .or. n < 1) return
    call b%ini(n,m)
    b%m = transpose(a%m)
  end function trans_dmat

  function inverse_dmat(r) result(s)
    class(DMat), intent(in) :: r
    type(DMat) :: s
    real(real64), allocatable :: a(:,:)
    real(real64), allocatable :: work(:)
    integer(int32), allocatable :: ipvt(:)
    integer(int32) :: info, n
    if( r%n_row /= r%n_col ) then
      write(*,*) "Error in inverse_dmat"
      return
    end if
    n = r%n_row
    if(n < 1) return
    call s%ini(n,n)
    allocate(work(n*n),ipvt(n))
    allocate(a(n,n))
    a = r%m
    call dgetrf(n,n,a,n,ipvt,info)
    call dgetri(n,a,n,ipvt,work,n**2,info)
    s%m = a
    deallocate(a,work,ipvt)
  end function inverse_dmat

  function determinant_dmat(r) result(d)
    class(DMat), intent(in) :: r
    real(real64) :: d
    integer(int32) :: n, i, info
    real(real64), allocatable :: a(:,:)
    integer(int32), allocatable :: ipiv(:)
    d = 0.0
    if( r%n_row /= r%n_col ) then
      write(*,*) "Error in determinant_dmat"
      return
    end if
    n = r%n_row
    if(n < 1) return
    allocate(ipiv(n), a(n,n))
    a = r%m
    call dgetrf(n, n, a, n, ipiv, info)
    if(info /= 0) then
      write(*,'(a, i3)') "error in det: info = ", info
      stop
    end if
    d = 1.0
    do i = 1, n
      if(ipiv(i) .ne. i) then
        d = -d * a(i, i)
      else
        d = d * a(i, i)
      end if
    end do
    deallocate(ipiv, a)
  end function determinant_dmat

  subroutine set_random_dmat(this, n_row, n_col)
    class(DMat), intent(inout) :: this
    integer(int32), intent(in) :: n_row, n_col
    integer(int32) :: i
    type(DVec) :: v
    call this%ini(n_row,n_col)
    if(n_row < 1 .or. n_col < 1) return
    do i = 1, n_col
      call v%rndm(n_row)
      this%m(:,i) = v%v(:)
      call v%fin()
    end do
  end subroutine set_random_dmat

  subroutine set_dmat_from_array( this, m )
    type(DMat), intent(inout) :: this
    real(real64), intent(in) :: m(:,:)
    call this%ini( size(m,1), size(m,2) )
    this%m(:,:) = m(:,:)
  end subroutine set_dmat_from_array

  subroutine set_diag_dmat(b, a)
    class(DMat), intent(inout) :: b
    type(DVec), intent(in) :: a
    integer(int32) :: n, i
    n = a%n_size
    if(n < 1) return
    call b%zeros(n,n)
    do i = 1, n
      b%M(i,i) = a%V(i)
    end do
  end subroutine set_diag_dmat

  function outer_product_dmat(a, b) result(c)
    type(DVec), intent(in) :: a, b
    type(DMat) :: c
    integer(int32) :: n, m
    n = a%n_size
    m = b%n_size
    call c%zeros(n,m)
    call dger(n, m, 1.d0, a%v, 1, b%v, 1, c%m, n)
  end function outer_product_dmat

  function dmat_dvec_product(a, b) result(c)
    type(DMat), intent(in) :: a
    type(DVec), intent(in) :: b
    type(DVec) :: c
    integer(int32) :: m, k
    m = a%n_row
    k = a%n_col
    if( k /= b%n_size ) then
      write(*, '(a)') 'Error in dmat_dvec_product'
      stop
    end if
    call c%ini(m)
    call dgemv('n',m,k,1.d0,a%m,m,b%v,1,0.d0,c%v,1)
  end function dmat_dvec_product

  function dvec_dmat_product(a, b) result(c)
    type(DMat), intent(in) :: b
    type(DVec), intent(in) :: a
    type(DVec) :: c
    integer(int32) :: m, k, n
    m = b%n_row
    k = b%n_col
    n = a%n_size
    if(m /= n) then
      write(*, '(a)') 'Error in dvec_dmat_product'
      stop
    end if
    call c%ini(k)
    call dgemv('t',m,k,1.d0,b%m,m,a%v,1,0.d0,c%v,1)
  end function dvec_dmat_product

  subroutine print_dmat(this, msg, iunit, binary)
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
  end subroutine print_dmat

  subroutine initialize_cmat( this, n_row, n_col )
    class(CMat), intent(inout) :: this
    integer(int32), intent(in) :: n_row, n_col
    if( allocated( this%m )) call this%fin()
    if( n_row < 1 .or. n_col < 1 ) return
    this%n_row = n_row
    this%n_col = n_col
    allocate( this%m( this%n_row, this%n_col ))
  end subroutine initialize_cmat

  subroutine finalize_cmat( this )
    class(CMat), intent(inout) :: this
    deallocate( this%m )
    this%n_row = 0
    this%n_col = 0
  end subroutine finalize_cmat

  subroutine zeros_cmat( this, n_row, n_col )
    class(CMat), intent(inout) :: this
    integer(int32), intent(in) :: n_row, n_col
    if(allocated( this%m )) call this%fin()
    if(n_row < 1 .or. n_col < 1) return
    call this%ini( n_row, n_col )
    this%m(:,:) = (0.d0,0.d0)
  end subroutine zeros_cmat

  subroutine eye_cmat( this, n )
    class(CMat), intent(inout) :: this
    integer(int32), intent(in) :: n
    integer(int32) :: i
    if(allocated( this%m )) call this%fin()
    if(n < 1) return
    call this%zeros( n, n )
    do i = 1, n
      this%m(i,i) = (1.d0,0.d0)
    end do
  end subroutine eye_cmat

  function get_size_cmat( this ) result( ns )
    class(CMat), intent(in) :: this
    integer(int32) :: ns(2)
    ns(1) = this%n_row
    ns(2) = this%n_col
  end function get_size_cmat

  subroutine copy_cmat( a, b )
    type(CMat), intent(inout) :: a
    class(*), intent(in) :: b
    select type( b )
    type is ( SMat )
      call a%ini( b%n_row, b%n_col )
      a%m(:,:) = dble(b%m(:,:))
    type is ( DMat )
      call a%ini( b%n_row, b%n_col )
      a%m(:,:) = b%m(:,:)
    type is ( CMat )
      call a%ini( b%n_row, b%n_col )
      a%m(:,:) = dble(b%m(:,:))
    end select
  end subroutine copy_cmat

  function product_cmat( a, b ) result(c)
    type(CMat), intent(in) :: a, b
    type(CMat) :: c
    integer(int32) :: m, k, n
    m = a%n_row
    k = a%n_col
    if(a%n_col /= b%n_row) then
      write(*, '(a)') 'Error in product_cmat'
      stop
    end if
    n = b%n_col
    call c%ini(m,n)
    call zgemm('n','n',m,n,k,(1.d0,0.d0),a%m,m,b%m,k,(0.d0,0.d0),c%m,m)
  end function product_cmat

  function sum_cmat( a, b ) result(c)
    type(CMat), intent(in) :: a, b
    type(CMat) :: c
    integer(int32) :: m, n, i
    if(a%n_row /= b%n_row .or. a%n_col /= b%n_col) then
      write(*, '(a)') 'Error in sum_cmat'
      stop
    end if
    m = a%n_row
    n = a%n_col
    if(m < 1 .or. n < 1) return
    call copy_cmat(c, a)
    do i = 1, n
      call zaxpy(m, (1.d0,0.d0), b%m(:,i), 1, c%m(:,i), 1)
    end do
  end function sum_cmat

  function subtract_cmat( a, b ) result(c)
    type(CMat), intent(in) :: a, b
    type(CMat) :: c
    integer(int32) :: m, n, i
    if(a%n_row /= b%n_row .or. a%n_col /= b%n_col) then
      write(*, '(a)') 'Error in subtract_cmat'
      stop
    end if
    m = a%n_row
    n = a%n_col
    if(m < 1 .or. n < 1) return
    call copy_cmat(c, a)
    do i = 1, n
      call zaxpy(m, (-1.d0,0.d0), b%m(:,i), 1, c%m(:,i), 1)
    end do
  end function subtract_cmat

  function scale_l_cmat( b, a ) result(c)
    type(CMat), intent(in) :: a
    complex(real64), intent(in) :: b
    type(CMat) :: c
    integer(int32) :: m, n, i
    m = a%n_row
    n = a%n_col
    if(m < 1 .or. n < 1) return
    call copy_cmat(c, a)
    do i = 1, n
      call zscal(m, b, c%m(:,i), 1)
    end do
  end function scale_l_cmat

  function scale_r_cmat( a, b ) result(c)
    type(CMat), intent(in) :: a
    complex(real64), intent(in) :: b
    type(CMat) :: c
    integer(int32) :: m, n, i
    m = a%n_row
    n = a%n_col
    if(m < 1 .or. n < 1) return
    call copy_cmat(c, a)
    do i = 1, n
      call zscal(m, b, c%m(:,i), 1)
    end do
  end function scale_r_cmat

  function scale_inv_cmat( a, b ) result(c)
    type(CMat), intent(in) :: a
    complex(real64), intent(in) :: b
    type(CMat) :: c
    integer(int32) :: m, n, i
    m = a%n_row
    n = a%n_col
    if(m < 1 .or. n < 1) return
    call copy_cmat(c, a)
    do i = 1, n
      call zscal(m, (1.d0,0.d0)/b, c%m(:,i), 1)
    end do
  end function scale_inv_cmat

  function trans_cmat( a ) result(b)
    class(CMat), intent(in) :: a
    type(CMat) :: b
    integer(int32) :: n, m
    m = a%n_row
    n = a%n_col
    if(m < 1 .or. n < 1) return
    call b%ini(n,m)
    b%m = transpose(a%m)
  end function trans_cmat

  function c_conjugate_cmat( a ) result(b)
    class(CMat), intent(in) :: a
    type(CMat) :: b
    integer(int32) :: n, m
    m = a%n_row
    n = a%n_col
    if(m < 1 .or. n < 1) return
    call b%ini(n,m)
    b%m = conjg(a%m)
  end function c_conjugate_cmat

  function h_conjugate_cmat( a ) result(b)
    class(CMat), intent(in) :: a
    type(CMat) :: b
    if( a%n_row<1 .or.  a%n_col<1) return
    b = a%star()
    b = b%t()
  end function h_conjugate_cmat

  function inverse_cmat(r) result(s)
    class(CMat), intent(in) :: r
    type(CMat) :: s
    complex(real64), allocatable :: a(:,:)
    complex(real64), allocatable :: work(:)
    integer(int32), allocatable :: ipvt(:)
    integer(int32) :: info, n
    if( r%n_row /= r%n_col ) then
      write(*,*) "Error in inverse_cmat"
      return
    end if
    n = r%n_row
    if(n < 1) return
    call s%ini(n,n)
    allocate(work(n*n),ipvt(n))
    allocate(a(n,n))
    a = r%m
    call zgetrf(n,n,a,n,ipvt,info)
    call zgetri(n,a,n,ipvt,work,n**2,info)
    s%m = a
    deallocate(a,work,ipvt)
  end function inverse_cmat

  function determinant_cmat(r) result(d)
    class(CMat), intent(in) :: r
    complex(real64) :: d
    integer(int32) :: n, i, info
    complex(real64), allocatable :: a(:,:)
    integer(int32), allocatable :: ipiv(:)
    d = 0.0
    if( r%n_row /= r%n_col ) then
      write(*,*) "Error in determinant_cmat"
      return
    end if
    n = r%n_row
    if(n < 1) return
    allocate(ipiv(n), a(n,n))
    a = r%m
    call zgetrf(n, n, a, n, ipiv, info)
    if(info /= 0) then
      write(*,'(a, i3)') "error in det: info = ", info
      stop
    end if
    d = 1.0
    do i = 1, n
      if(ipiv(i) .ne. i) then
        d = -d * a(i, i)
      else
        d = d * a(i, i)
      end if
    end do
    deallocate(ipiv, a)
  end function determinant_cmat

  subroutine set_random_cmat(this, n_row, n_col)
    class(CMat), intent(inout) :: this
    integer(int32), intent(in) :: n_row, n_col
    integer(int32) :: i
    type(CVec) :: v
    call this%ini(n_row,n_col)
    if(n_row < 1 .or. n_col < 1) return
    do i = 1, n_col
      call v%rndm(n_row)
      this%m(:,i) = v%v(:)
      call v%fin()
    end do
  end subroutine set_random_cmat

  subroutine set_cmat_from_array( this, m )
    type(CMat), intent(inout) :: this
    complex(real64), intent(in) :: m(:,:)
    call this%ini( size(m,1), size(m,2) )
    this%m(:,:) = m(:,:)
  end subroutine set_cmat_from_array

  subroutine set_diag_cmat(b, a)
    class(CMat), intent(inout) :: b
    type(CVec), intent(in) :: a
    integer(int32) :: n, i
    n = a%n_size
    if(n < 1) return
    call b%zeros(n,n)
    do i = 1, n
      b%M(i,i) = a%V(i)
    end do
  end subroutine set_diag_cmat

  function outer_product_cmat(a, b) result(c)
    type(CVec), intent(in) :: a, b
    type(CMat) :: c
    integer(int32) :: n, m
    n = a%n_size
    m = b%n_size
    call c%zeros(n,m)
    call zgeru(n, m, (1.d0,0.d0), a%v, 1, b%v, 1, c%m, n)
  end function outer_product_cmat

  function cmat_cvec_product(a, b) result(c)
    type(CMat), intent(in) :: a
    type(CVec), intent(in) :: b
    type(CVec) :: c
    integer(int32) :: m, k
    m = a%n_row
    k = a%n_col
    if( k /= b%n_size ) then
      write(*, '(a)') 'Error in cmat_cvec_product'
      stop
    end if
    call c%ini(m)
    call zgemv('n',m,k,(1.d0,0.d0),a%m,m,b%v,1,(0.d0,0.d0),c%v,1)
  end function cmat_cvec_product

  function cvec_cmat_product(a, b) result(c)
    type(CMat), intent(in) :: b
    type(CVec), intent(in) :: a
    type(CVec) :: c
    integer(int32) :: m, k, n
    m = b%n_row
    k = b%n_col
    n = a%n_size
    if(m /= n) then
      write(*, '(a)') 'Error in cvec_cmat_product'
      stop
    end if
    call c%ini(k)
    call zgemv('t',m,k,(1.d0,0.d0),b%m,m,a%v,1,(0.d0,0.d0),c%v,1)
  end function cvec_cmat_product

  subroutine print_cmat(this, msg, iunit, binary)
    class(cmat), intent(in) :: this
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
        write(unt,'(a)') 'Real:'
        do i=1,n
          write(unt,cfmt) dble(this%m(i,:))
        end do
        write(unt,'(a)') 'Imag:'
        do i=1,n
          write(unt,cfmt) aimag(this%m(i,:))
        end do
      else
        do i = 1, n
          do j = 1, m
            write(unt,'(2i8,2f16.6)') i,j,this%m(i,j)
          end do
        end do
      end if
    end if
  end subroutine print_cmat
end module matrix_definitions
