module mat_vec_double
  use General
  use vector_double, only: DVec
  use matrix_double, only: DMat
  implicit none
  public :: OuterProductD, MVProductD, VMProductD
contains
  type(DMat) function OuterProductD(a, b) result(c)
    type(DVec), intent(in) :: a, b
    integer(int32) :: n, m
    n = size(a%v)
    m = size(b%v)
    call c%ini(n,m)
    call dger(n, m, 1.d0, a%v, 1, b%v, 1, c%m, n)
  end function OuterProductD

  type(DVec) function MVProductD(a, b) result(c)
    type(DMat), intent(in) :: a
    type(DVec), intent(in) :: b
    integer(int32) :: m, k
    m = size(a%m, 1)
    k = size(a%m, 2)
    if(size(a%m, 2) /= size(b%v)) then
      write(*, '(a)') 'Error in MVProduct'
      stop
    end if
    call c%Ini(m)
    call dgemv('n',m,k,1.d0,a%m,m,b%v,1,0.d0,c%v,1)
  end function MVProductD

  type(DVec) function VMProductD(a, b) result(c)
    type(DMat), intent(in) :: b
    type(DVec), intent(in) :: a
    integer(int32) :: m, k, n
    m = size(b%m, 1)
    k = size(b%m, 2)
    n = size(a%v, 1)
    if(m /= n) then
      write(*, '(a)') 'Error in VMProduct'
      stop
    end if
    call c%Ini(k)
    call dgemv('t',m,k,1.d0,b%m,m,a%v,1,0.d0,c%v,1)
  end function VMProductD
end module mat_vec_double
