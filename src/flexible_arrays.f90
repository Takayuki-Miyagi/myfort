module flexible_arrays
  use basic_types
  use linear_algebra

  type :: list
    type(list), pointer :: next => null()
    class(*), allocatable :: element
  end type list

contains
  subroutine init_list( this )
    class(list), intent(inout), pointer :: this
    allocate( this )
    nullify( this%next )
  end subroutine init_list

  subroutine free_list( this )
    class(list), intent(inout) :: this
  end subroutine free_list

  subroutine push_back( this )
    class(list), intent(inout) :: this
  end subroutine push_back

end module flexible_arrays
