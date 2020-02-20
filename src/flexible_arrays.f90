module flexible_arrays
  use basic_types
  use vector_single
  use vector_double
  use vector_complex
  use matrix_single
  use matrix_double
  use matrix_complex

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
