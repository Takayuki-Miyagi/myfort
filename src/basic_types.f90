module basic_types
  use, intrinsic :: iso_fortran_env
  implicit none
  private :: copy_int32
  private :: copy_int64
  private :: copy_real32
  private :: copy_real64
  private :: copy_complex32
  private :: copy_complex64
  private :: copy_char
  private :: copy_logical
  integer(int32) :: iseed(4) = (/3239, 4241, 1903, 1093/) ! seed of random numbers
  interface assignment(=)
    module procedure :: copy_int32
    module procedure :: copy_int64
    module procedure :: copy_real32
    module procedure :: copy_real64
    module procedure :: copy_complex32
    module procedure :: copy_complex64
    module procedure :: copy_char
    module procedure :: copy_logical
  end interface assignment(=)
contains
  subroutine copy_int32( r, this )
    class(*), intent(in) :: this
    integer(int32) :: r
    select type( this )
    type is ( integer(int32) )
      r = this
      class default
      write(*,*) "disagreement of type"
    end select
  end subroutine copy_int32

  subroutine copy_int64( r, this )
    class(*), intent(in) :: this
    integer(int64) :: r
    select type( this )
    type is ( integer(int64) )
      r = this
      class default
      write(*,*) "disagreement of type"
    end select
  end subroutine copy_int64

  subroutine copy_real32( r, this )
    class(*), intent(in) :: this
    real(real32) :: r
    select type( this )
    type is ( real(real32) )
      r = this
      class default
      write(*,*) "disagreement of type"
    end select
  end subroutine copy_real32

  subroutine copy_real64( r, this )
    class(*), intent(in) :: this
    real(real64) :: r
    select type( this )
    type is ( real(real64) )
      r = this
      class default
      write(*,*) "disagreement of type"
    end select
  end subroutine copy_real64

  subroutine copy_complex32( r, this )
    class(*), intent(in) :: this
    complex(real32) :: r
    select type( this )
    type is ( complex(real32) )
      r = this
      class default
      write(*,*) "disagreement of type"
    end select
  end subroutine copy_complex32

  subroutine copy_complex64( r, this )
    class(*), intent(in) :: this
    complex(real64) :: r
    select type( this )
    type is ( complex(real64) )
      r = this
      class default
      write(*,*) "disagreement of type"
    end select
  end subroutine copy_complex64

  subroutine copy_char( r, this )
    class(*), intent(in) :: this
    character(:), allocatable :: r
    select type( this )
    type is ( character(LEN=*) )
      allocate(r, source=this)
      class default
      write(*,*) "disagreement of type"
    end select
  end subroutine copy_char

  subroutine copy_logical( r, this )
    class(*), intent(in) :: this
    logical, intent(out) :: r
    select type( this )
    type is ( logical )
      r = this
      class default
      write(*,*) "disagreement of type"
    end select
  end subroutine copy_logical
end module basic_types
