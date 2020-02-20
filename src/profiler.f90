module profiler
  use omp_lib
  use basic_types
  use general
  implicit none

  public :: timer
  private :: InitProf
  private :: StartProf
  private :: FinProf
  private :: AddProf
  private :: PrintSummary
  private :: GetCurrentMemory
  private :: GetTargetMemory
  private :: CountUpMemory
  private :: CountDownMemory
  private :: PrintMemory

  type, private :: prof
    type(imap) :: counter
    type(dmap) :: timer
    type(imap) :: nmem
    type(dmap) :: memories
    real(real64) :: start_time
    real(real64) :: total_memory = 0.d0
    real(real64) :: max_memory = 0.d0
  contains
    procedure :: InitProf
    procedure :: StartProf
    procedure :: FinProf
    procedure :: AddProf
    procedure :: PrintSummary
    procedure :: GetCurrentMemory
    procedure :: GetTargetMemory
    procedure :: CountUpMemory
    procedure :: CountDownMemory
    procedure :: PrintMemory
    generic :: init => InitProf
    generic :: start => StartProf
    generic :: fin => FinProf
    generic :: add => AddProf
    generic :: prt => PrintSummary
    generic :: get_memory => GetCurrentMemory, GetTargetMemory
    generic :: countup_memory => CountUpMemory
    generic :: countdw_memory => CountDownMemory
    generic :: prt_memory => PrintMemory
  end type prof

  type(prof) :: timer

contains
  subroutine InitProf(this)
    class(prof), intent(inout) :: this

    this%start_time = omp_get_wtime()
    call this%counter%append('Total',1)
    call this%timer%append('Total',0.d0)
  end subroutine InitProf

  subroutine StartProf(this,key)
    class(prof), intent(inout) :: this
    character(*), intent(in) :: key
    call this%counter%append(key,0)
    call this%timer%append(key,0.d0)
  end subroutine StartProf

  subroutine AddProf(this,key,time_in)
    class(prof), intent(inout) :: this
    character(*), intent(in) :: key
    real(real64), intent(in) :: time_in
    integer(int32) :: n
    real(real64) :: time
    if(.not. this%counter%Search(key)) call this%start(key)
    n = this%counter%Get(key)
    time = this%timer%Get(key)
    call this%counter%append(key,n+1)
    call this%timer%append(key,time+time_in)
  end subroutine AddProf

  subroutine PrintSummary(this)
    class(prof), intent(inout) :: this
    integer(int32) :: n, i
    real(real64) :: r, ttotal
    call this%timer%append('Total',omp_get_wtime() - this%start_time)
    write(*,*)
    ttotal = this%timer%get('Total')
    n = int(ttotal)
    write(*,'(a,i5,a,i2.2,a,i2.2)') '    summary of time, total = ', n/3600, ':', mod(n, 3600)/60, ':', mod(n, 60)
    write(*,*)
    write(*,'(37x,a)') "time,    ncall, time/ncall,   ratio "
    r = ttotal
    do i = 1, this%counter%GetLen()
      if(this%counter%key(i) == 'Total') cycle
      call PrintEach(this%counter%key(i), this%counter%val(i), this%timer%val(i))
      r = r - this%timer%val(i)
    end do
    write(*,'(1a30, 1f12.3, 22x, 1f9.4)') "misc", r, r/ttotal
    write(*,*)
  contains
    subroutine PrintEach(title, ncall, time)
      character(*), intent(in) :: title
      integer(int32), intent(in) :: ncall
      real(real64), intent(in) :: time
      write(*,'(1a30, 1f12.3, 1i10, 1f12.5, 1f9.4)') title,  &
          time, ncall, time/ncall, &
          time/ttotal
    end subroutine PrintEach

  end subroutine PrintSummary

  function GetCurrentMemory(this) result(r)
    class(prof), intent(inout) :: this
    real(real64) :: r
    r = this%total_memory
  end function GetCurrentMemory

  function GetTargetMemory(this,key) result(r)
    class(prof), intent(inout) :: this
    character(*), intent(in) :: key
    real(real64) :: r
    r = this%memories%get(key)
  end function GetTargetMemory

  subroutine CountUpMemory(this, key, memory_in)
    class(prof), intent(inout) :: this
    character(*), intent(in) :: key
    real(real64), intent(in) :: memory_in
    integer(int32) :: n
    real(real64) :: mem

    if(.not. this%memories%Search(key)) then
      call this%nmem%append(key,0)
      call this%memories%append(key,0.d0)
    end if
    n = this%nmem%Get(key)
    mem = this%memories%Get(key)
    call this%nmem%append(key,n+1)
    call this%memories%append(key,mem+memory_in)
  end subroutine CountUpMemory

  ! call timer%cmemory before using
  ! msg is object name
  subroutine CountDownMemory(this, key, memory_in)
    class(prof), intent(inout) :: this
    character(*), intent(in) :: key
    real(real64), intent(in) :: memory_in
    integer(int32) :: n
    real(real64) :: mem

    n = this%nmem%Get(key)
    mem = this%memories%Get(key)
    call this%nmem%append(key,n-1)
    call this%memories%append(key,mem-memory_in)
  end subroutine CountDownMemory

  subroutine PrintMemory(this)
    class(prof), intent(inout) :: this
    integer(int32) :: i
    real(real64) :: mtotal, a

    mtotal = 0.d0
    mtotal = this%max_memory / dble(1024) ** 2

    write(*,'(a,f12.6,a)') '      Max used memory = ', mtotal, ' GB'
    write(*,*)
    write(*,'(30x,a)') "memory (GB),    ncall,   memory/ncall (GB)"
    do i = 1, this%nmem%GetLen()
      if(this%nmem%key(i) == 'Total') cycle
      a = this%memories%val(i) / dble(1024) ** 2
      call PrintEach(this%nmem%key(i), this%nmem%val(i), a)
    end do
    write(*,*)
  contains
    subroutine PrintEach(title, ncall, amnt)
      character(*), intent(in) :: title
      integer(int32), intent(in) :: ncall
      real(real64), intent(in) :: amnt

      write(*,'(1a30, 1f12.3, 1i10, 1f20.5)') title,  amnt, ncall, amnt/ncall
    end subroutine PrintEach

  end subroutine PrintMemory

  subroutine FinProf(this)
    class(prof), intent(inout) :: this
    call this%prt_memory()
    call this%prt()
  end subroutine FinProf

end module profiler
