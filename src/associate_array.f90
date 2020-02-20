module associate_array
  use basic_types
  implicit none

  type :: list
    integer(int32) :: n_size
    class(*), allocatable :: elements(:)
  end type list

  type :: tree
    type(tree), pointer :: left => null()
    type(tree), pointer :: right => null()
    class(*), allocatable :: key
    class(*), allocatable :: val
  end  type tree
contains
  subroutine add_element( this, key, val)
    class(tree), intent(inout) :: this
    class(*), intent(in) :: key
    class(*), intent(in) :: val
    allocate(this%key, source=key)
    allocate(this%val, source=val)
  end subroutine add_element

  subroutine print_key( this, iunit )
    class(tree), intent(in), target :: this
    integer(int32), intent(in), optional :: iunit
    call print_basic_class( this%key, iunit )
  end subroutine print_key

  subroutine print_val( this, iunit )
    class(tree), intent(in), target :: this
    integer(int32), intent(in), optional :: iunit
    call print_basic_class( this%val, iunit )
  end subroutine print_val

  subroutine print_basic_class( p, iunit )
    class(*), intent(in) :: p
    integer(int32), intent(in), optional :: iunit
    integer(int32) :: write_unit = 6
    if(present(iunit)) write_unit = iunit
    select type( p )
    type is ( integer(int32) )
      write(write_unit,*) p
    type is ( integer(int64) )
      write(write_unit,*) p
    type is ( real(real32) )
      write(write_unit,*) p
    type is ( real(real64) )
      write(write_unit,*) p
    type is ( complex(real32) )
      write(write_unit,*) p
    type is ( complex(real64) )
      write(write_unit,*) p
    type is ( character(LEN=*) )
      write(write_unit,*) p
    type is ( logical )
      write(write_unit,*) p
    end select
  end subroutine print_basic_class
end module associate_array
