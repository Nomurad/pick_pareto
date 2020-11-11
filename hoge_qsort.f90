program hoge
  implicit none
  interface
    function comp(x, y) result(r)
      integer, intent(in) :: x(:), y
      logical,allocatable :: r(:)
    end function comp
  end interface
  integer :: x(35)
 
  x = ([3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3, 2, 3, 8, 4, 6, 2, 6, 4, 3, 3, 8, 3, 2, 7, 9, 5, 0, 2, 8])
  write(6,'(a, 35i2)') "qsortf",qsortf(x, ltcomp)
 
contains
  recursive function qsortf(a, f) result(r)
    integer, intent(in) :: a(:)
    procedure(comp) :: f
    integer, allocatable :: r(:)
    logical, allocatable :: mask(:)
 
    if(size(a) < 2)then
      r = a
    else
      mask = f(a(2:), a(1))
      r = [qsortf(pack(a(2:), mask), f), a(1), qsortf(pack(a(2:), .not. mask), f)]
    endif
  end function qsortf
 
  function ltcomp(x, y) result(r)
      integer, intent(in) :: x(:), y
      logical,allocatable :: r(:)
      r = x < y
  end function ltcomp
end program hoge
