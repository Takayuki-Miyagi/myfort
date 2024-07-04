program test
  use myfort
  implicit none
  call timer%init()

  !call test_LA()
  call test_renormalization()
  !call test_myfort_types()

  !call test_wfs()
  !call test_iteration_methods()

  call timer%prt()
contains

!  subroutine test_myfort_types()
!    type(sys) :: s
!    type(tree) :: map
!    class(*), pointer :: a
!    character(:), allocatable :: b
!    type(DVec) :: v,vv
!    v = [1.d0,2.d0,3.d0]
!    allocate(a, source=v)
!    vv = a
!    call vv%prnt()
!  end subroutine test_myfort_types

  subroutine test_renormalization()
    type(DMat) :: m, gen
    type(SRGSolver) :: srg
    type(LeeSuzukiSolver) :: ls
    integer :: n, i
    real(8) :: alpha=1.d0
    n = 10
    call m%rndm(n,n)
    call gen%zeros(n,n)
    m = m + m%t()
    do i = 1, n
      gen%m(i,i) = m%m(i,i)
    end do
    call srg%init( m, "Hflow" )
    call srg%h%prnt(" before flow ")
    call gen%prnt(" generator ")
    call srg%SRGFlow( m, gen, alpha )
    call srg%h%prnt(" after flow ")
    call srg%fin()

    call m%rndm(n,n)
    m = m + m%t()
    call ls%init( m )
    call m%prnt( "before LS decoupling" )
    call ls%LeeSuzuki( 2, 4, n-2-4, m )
    call ls%H%prnt( "after LS decoupling" )
    call ls%fin()

  end subroutine test_renormalization

  subroutine test_LA()
    type(SVec) :: sv1, sv2
    type(DVec) :: dv1, dv2
    type(CVec) :: cv1, cv2
    type(SMat) :: sm1, sm2
    type(DMat) :: dm1, dm2
    type(CMat) :: cm1, cm2
    type(EigenSolSymD) :: dsol


    write(*,*) " test linear-algebra library "
    write(*,*) " single precision"
    write(*,*)
    call sv1%rndm(3)
    call sv2%rndm(3)
    call sv1%prnt("1st svec")
    call sv2%prnt("2nd svec")
    sm1 = sv1 .x. sv2
    call sm1%prnt("v1 x v2")
    sm2 = sv2 .x. sv1
    call sm2%prnt("v2 x v1")
    sm1 = sm1 * sm2
    call sm1%prnt("m1 * m2")
    sv1 = sm2 * sv2
    call sv1%prnt("(v2 x v1) * v2")
    sv1 = sv2 * sm2
    call sv1%prnt("v2 * (v2 x v1)")

    call sm1%rndm(3,3)
    sm2 = sm1%inv()
    sm1 = sm1 * sm2
    call sm1%prnt(" Mat * Mat-1")

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

    call dm1%rndm(3,3)
    dm2 = dm1 + dm1%t()
    call dm2%prnt(" Mat ")
    call dsol%init( dm2 )
    call dsol%DiagSym( dm2 )
    call dsol%eig%prnt("EigenValues")
    call dsol%vec%prnt("EigenVectors")
    call dm1%diag( dsol%eig )
    dm2 = dsol%vec * dm1 * dsol%vec%T()
    call dm2%prnt(" U D U-1 ")
    call dm1%rndm(3,3)
    dm2 = dm1%inv()
    dm1 = dm1 * dm2
    call dm1%prnt(" Mat * Mat-1")
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
