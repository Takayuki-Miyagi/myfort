module mat_vec_single
  use basic_types
  use vector_single, only: SVec
  use matrix_single, only: SMat
  implicit none
  public :: OuterProductS, MVProductS, VMProductS
contains
  type(SMat) function OuterProductS(a, b) result(c)
    type(SVec), intent(in) :: a, b
    integer(int32) :: n, m
    n = size(a%v)
    m = size(b%v)
    call c%ini(n,m)
    call sger(n, m, 1.0, a%v, 1, b%v, 1, c%m, n)
  end function OuterProductS

  type(SVec) function MVProductS(a, b) result(c)
    type(SMat), intent(in) :: a
    type(SVec), intent(in) :: b
    integer(int32) :: m, k
    m = size(a%m, 1)
    k = size(a%m, 2)
    if(size(a%m, 2) /= size(b%v)) then
      write(*, '(a)') 'Error in MVProduct'
      stop
    end if
    call c%Ini(m)
    call sgemv('n',m,k,1.0,a%m,m,b%v,1,0.0,c%v,1)
  end function MVProductS

  type(SVec) function VMProductS(a, b) result(c)
    type(SMat), intent(in) :: b
    type(SVec), intent(in) :: a
    integer(int32) :: m, k, n
    m = size(b%m, 1)
    k = size(b%m, 2)
    n = size(a%v, 1)
    if(m /= n) then
      write(*, '(a)') 'Error in VMProduct'
      stop
    end if
    call c%Ini(k)
    call sgemv('t',m,k,1.0,b%m,m,a%v,1,0.0,c%v,1)
  end function VMProductS
end module mat_vec_single
