program test
  use myfort
  implicit none
  type(sys) :: s
  type(tree) :: map
  class(*), pointer :: a
  character(:), allocatable :: b
  type(DVec) :: v,vv
  call timer%init()

  call test_LA()

  !call test_wfs()
  !call test_iteration_methods()

  !call add_element(map, "aaaaa", "bbbbb")
  !call print_key(map)
  !call print_val(map)
  !v = [1.d0,2.d0,3.d0]
  !allocate(a, source=v)
  !vv = a
  !call vv%prnt()

  call timer%prt()
contains

  subroutine test_LA()
    type(SVec) :: sv1, sv2
    type(DVec) :: dv1, dv2
    type(CVec) :: cv1, cv2
    type(SMat) :: sm1, sm2
    type(DMat) :: dm1, dm2
    type(CMat) :: cm1, cm2


    write(*,*) " test linear-algebra library "
    write(*,*) " single precision"
    write(*,*)
    call sv1%rndm(3)
    call sv2%rndm(3)
    call sv1%prnt("1st svec")
    call sv2%prnt("2nd svec")
    sm1 = sv1 .x. sv2
    sm2 = sv2 .x. sv1
    call sm1%prnt("v1 x v2")
    call sm2%prnt("v2 x v1")
    sm1 = sm1 * sm2
    call sm1%prnt("m1 * m2")
    sv1 = sm2 * sv2
    call sv1%prnt("(v2 x v1) * v2")
    sv1 = sv2 * sm2
    call sv1%prnt("v2 * (v2 x v1)")

    write(*,*)
    write(*,*) " double precision"
    write(*,*)
    call dv1%rndm(3)
    call dv2%rndm(3)
    call dv1%prnt("1st dvec")
    call dv2%prnt("2nd dvec")
    dm1 = dv1 .x. dv2
    dm2 = dv2 .x. dv1
    call dm1%prnt("v1 x v2")
    call dm2%prnt("v2 x v1")
    dm1 = dm1 * dm2
    call dm1%prnt("m1 * m2")
    dv1 = dm2 * dv2
    call dv1%prnt("(v2 x v1) * v2")
    dv1 = dv2 * dm2
    call dv1%prnt("v2 * (v2 x v1)")
  end subroutine test_LA

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
