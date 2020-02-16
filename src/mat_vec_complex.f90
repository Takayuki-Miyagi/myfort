module mat_vec_complex
  use General
  use vector_complex, only: CVec
  use matrix_complex, only: CMat
  implicit none
  public :: OuterProductC, MVProductC, VMProductC
contains
  type(CMat) function OuterProductC(a, b) result(c)
    type(CVec), intent(in) :: a, b
    integer(int32) :: n, m
    n = size(a%v)
    m = size(b%v)
    call c%ini(n,m)
    call zgeru(n, m, 1.d0, a%v, 1, b%v, 1, c%m, n)
  end function OuterProductC

  type(CVec) function MVProductC(a, b) result(c)
    type(CMat), intent(in) :: a
    type(CVec), intent(in) :: b
    integer(int32) :: m, k
    m = size(a%m, 1)
    k = size(a%m, 2)
    if(size(a%m, 2) /= size(b%v)) then
      write(*, '(a)') 'Error in MVProduct'
      stop
    end if
    call c%Ini(m)
    call zgemv('n',m,k,1.d0,a%m,m,b%v,1,0.d0,c%v,1)
  end function MVProductC

  type(CVec) function VMProductC(a, b) result(c)
    type(CMat), intent(in) :: b
    type(CVec), intent(in) :: a
    integer(int32) :: m, k, n
    m = size(b%m, 1)
    k = size(b%m, 2)
    n = size(a%v, 1)
    if(m /= n) then
      write(*, '(a)') 'Error in VMProduct'
      stop
    end if
    call c%Ini(k)
    call zgemv('t',m,k,1.d0,b%m,m,a%v,1,0.d0,c%v,1)
  end function VMProductC
end module mat_vec_complex
