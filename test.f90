program test
  use myfort
  implicit none
  type(sys) :: s
  call timer%init()

  !call test_SVec()
  !call test_SMat()

  !call test_DVec()
  !call test_DMat()

  call test_wfs()
  call test_iteration_methods()

  call timer%prt()
contains
  subroutine test_SVec()
    type(SVec) :: a, b, c
    real(4) :: v
    integer :: n = 2
    write(*,*) "# test single-vector"
    call a%Random(n)
    call b%Random(n)

    call a%prt('a')
    call b%prt('b')
    write(*,'(a, f12.6)') 'a * b = ', a * b
  end subroutine test_SVec

  subroutine test_SMat()
    type(SMat) :: a, b, c
    integer :: n = 2

    write(*,*) "# test single-matrix"
    call a%Random(n,n)
    call b%Random(n,n)
    c = a + b
    call c%prt(iunit=6,msg='a + b')
  end subroutine test_SMat

  subroutine test_DVec()
    type(DVec) :: a, b, c
    integer :: n = 2
    write(*,*) "# test double-vector"
    call a%Random(n)
    call b%Random(n)
    write(*,'(a, f12.6)') 'a * b = ', a * b
  end subroutine test_DVec

  subroutine test_DMat()
    type(DMat) :: a, b, c
    integer :: n = 2

    write(*,*) "# test double-matrix"
    call a%Random(n,n)
    call b%Random(n,n)
    c = a + b
    call c%prt(iunit=6,msg='a + b')
  end subroutine test_DMat

  subroutine test_wfs()
    integer :: i, nmesh = 1000
    real(real64) :: rmax=1.d-3
    integer :: N=20
    integer :: L=40
    write(*,*) "# test wave function "
    open(15,file="wf.txt")
    do i = 1, nmesh
      !write(15,*) i, rmax/dble(nmesh) * dble(i), laguerre_radial_wf_norm(n,l,1.d0, rmax/dble(nmesh) * dble(i))**2
      !write(15,*) i, rmax/dble(nmesh) * dble(i), mom_laguerre_radial_wf_norm(n,l,1.d0, rmax/dble(nmesh) * dble(i))**2
      write(15,*) i, rmax/dble(nmesh) * dble(i), spherical_bessel( L, rmax/dble(nmesh)*dble(i) )
    end do
  end subroutine test_wfs

  subroutine test_iteration_methods()
    integer :: i, n = 1000, ndim
    real(8), allocatable :: fk(:), xk(:)
    type(Optimizer) :: opt
    write(*,*) "# test iteration method"
    ndim = 2
    allocate(fk(ndim), xk(ndim))
    !call opt%init(n=ndim, method="broyden", a = 0.9d0, m = 10)
    !call opt%init(n=ndim, method="bfgs", a = 0.1d0, m = 10)
    call opt%init(n=ndim, method="mbroyden", a = 0.7d0, m = 10)
    !call opt%init(n=ndim, method="l-bfgs", a = 0.7d0, m = 10)
    call opt%set_init([1.d0,1.d0])
    do i = 1, n
      call func(xk, fk)
      call opt%GetNextInput(fk, i)
      xk = opt%xj%v
      call func(xk, fk)
      write(*,'(i4, 4f12.5)') i, xk(:), fk(:) - xk(:)
      if(opt%r < 1.d-8) exit
    end do
    call func(xk, fk)
    write(*,'(2f12.5)') fk(:)
    call opt%Fin()
  end subroutine test_iteration_methods

  subroutine func(x, f)
    real(8), intent(in) :: x(:)
    real(8), intent(out) :: f(:)
    f(:) = cos(x(:))
    if(size(f) < 2) return
    f(2) = sin(x(2))
  end subroutine func
end program test
