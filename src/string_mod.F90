module string_mod

  implicit none

  private

  public to_string
  public string_split

  interface to_string
    module procedure integer_to_string
    module procedure real4_to_string
    module procedure real8_to_string
    module procedure logical_to_string
  end interface to_string

contains

  function integer_to_string(x) result(res)

    integer, intent(in) :: x
    character(:), allocatable :: res

    character(range(x)+2) tmp

    write(tmp, '(i0)') x
    res = tmp

  end function integer_to_string

  function real4_to_string(x, decimal_width) result(res)

    real(4), intent(in) :: x
    integer, intent(in), optional :: decimal_width
    character(range(x)+2) res

    integer w, y
    character(10) fmt
    character(range(x)+2) tmp

    w = merge(decimal_width, 4, present(decimal_width))
    write(fmt, "('(G', I0, '.', I0, ')')") w, w-6
    write(res, fmt) x

  end function real4_to_string

  function real8_to_string(x, decimal_width) result(res)

    real(8), intent(in) :: x
    integer, intent(in), optional :: decimal_width
    character(range(x)+2) res

    integer w, y
    character(10) fmt

    w = merge(decimal_width, 4, present(decimal_width))
    write(fmt, "('(E', I0, '.', I0, ')')") w, w-6
    write(res, fmt) x

  end function real8_to_string

  function logical_to_string(x) result(res)

    logical, intent(in) :: x
    character(5), allocatable :: res

    res = merge('true ', 'false', x)

  end function logical_to_string

  function string_split(x, at, delim) result(res)

    character(*), intent(in) :: x
    integer, intent(in) :: at 
    character(*), intent(in), optional :: delim
    character(:), allocatable :: res

    character(:), allocatable :: delim_
    integer start_pos, end_pos, count

    if (present(delim)) then
      delim_ = delim
    else
      delim_ = ' '
    end if
    start_pos = 1
    end_pos = 1
    count = 0
    do
      end_pos = index(trim(x(start_pos:)), delim_)
      if (end_pos == 0) then
        end_pos = len(x) + 1
        exit
      end if
      count = count + 1
      end_pos = end_pos + start_pos - 1
      if (count == at) exit
      start_pos = end_pos + 1
    end do
    if (count == 0) then
      res = ''
    else
      res = x(start_pos:end_pos-1)
    end if

  end function string_split

end module string_mod
