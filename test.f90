program test
  use myfort
  implicit none
  type(sys) :: s
  call timer%init()

  write(*,*) "# test single-vector"
  call test_SVec()
  write(*,*) "# test single-matrix"
  call test_SMat()

  write(*,*) "# test double-vector"
  call test_DVec()
  write(*,*) "# test double-matrix"
  call test_DMat()

  !call test_wfs()

  call timer%prt()
contains
  subroutine test_SVec()
    type(SVec) :: a, b, c
    real(4) :: v
    integer :: n = 2
    call a%Random(n)
    call b%Random(n)

    call a%prt('a')
    call b%prt('b')
    c = a + b
    call c%prt('a + b')
    c = a - b
    call c%prt('a - b')

    write(*,'(a, f12.6)') 'a * b = ', a * b
    write(*,'(a, f12.6)') 'Norm a = ', a%Nrm()
    write(*,'(a, f12.6)') 'Norm b = ', b%Nrm()
    write(*,'(a, f12.6)') 'Norm**2 a = ', a%Nrm2()
    write(*,'(a, f12.6)') 'Norm**2 b = ', b%Nrm2()

    c = a%blk(1,1)
    call c%prt('a(1,1)')
    c = a%blk(2,2)
    call c%prt('a(2,2)')

  end subroutine test_SVec

  subroutine test_SMat()
    type(SMat) :: a, b, c
    integer :: n = 2

    call a%Random(n,n)
    call b%Random(n,n)
    c = a + a%T()
    call c%prt(iunit=6,msg='a + a^t')

    c = a + b
    call c%prt(iunit=6,msg='a + b')

    c = a - b
    call c%prt(iunit=6,msg='a - b')

    c = a * a%inv()
    call c%prt(iunit=6,msg='a * a^-1')

    c = b * a%inv()
    call c%prt(iunit=6,msg='b * a^-1')

    c = a%blk(1,1,1,n)
    call c%prt(iunit=6,msg='1st column vector')

    c = a%blk(1,n,1,1)
  end subroutine test_SMat

  subroutine test_DVec()
    type(DVec) :: a, b, c
    integer :: n = 2
    call a%Random(n)
    call b%Random(n)

    call a%prt('a')
    call b%prt('b')
    c = a + b
    call c%prt('a + b')
    c = a - b
    call c%prt('a - b')

    write(*,'(a, f12.6)') 'a * b = ', a * b
    write(*,'(a, f12.6)') 'Norm a = ', a%Nrm()
    write(*,'(a, f12.6)') 'Norm b = ', b%Nrm()
    write(*,'(a, f12.6)') 'Norm**2 a = ', a%Nrm2()
    write(*,'(a, f12.6)') 'Norm**2 b = ', b%Nrm2()

    c = a%blk(1,1)
    call c%prt('a(1,1)')
    c = a%blk(2,2)
    call c%prt('a(2,2)')
  end subroutine test_DVec

  subroutine test_DMat()
    type(DMat) :: a, b, c
    integer :: n = 2

    call a%Random(n,n)
    call b%Random(n,n)

    call a%prt(iunit=6,msg='a')
    call b%prt(iunit=6,msg='b')

    c = a + a%T()
    call c%prt(iunit=6,msg='a + a^t')

    c = a + b
    call c%prt(iunit=6,msg='a + b')

    c = a - b
    call c%prt(iunit=6,msg='a - b')

    c = a * a%inv()
    call c%prt(iunit=6,msg='a * a^-1')

    c = b * a%inv()
    call c%prt(iunit=6,msg='b * a^-1')

    c = a%blk(1,1,1,n)
    call c%prt(iunit=6,msg='1st row vector')
  end subroutine test_DMat

  subroutine test_wfs()
    integer :: i, nmesh = 1000
    real(real64) :: rmax=1.d3
    integer :: N=20
    integer :: L=20
    open(15,file="wf.txt")
    do i = 1, nmesh
      !write(15,*) i, rmax/dble(nmesh) * dble(i), laguerre_radial_wf_norm(n,l,1.d0, rmax/dble(nmesh) * dble(i))**2
      write(15,*) i, rmax/dble(nmesh) * dble(i), mom_laguerre_radial_wf_norm(n,l,1.d0, rmax/dble(nmesh) * dble(i))**2
    end do
  end subroutine test_wfs
end program test
