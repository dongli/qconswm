module log_mod

  use map_mod
  use time_mod
  use string_mod

  implicit none

  type(map_type) diags

contains

  subroutine log_add_diag(name, value)

    character(*), intent(in) :: name
    real, intent(in) :: value

    call diags%insert(name, value)

  end subroutine log_add_diag

  subroutine log_notice(message, file, line)

    character(*), intent(in) :: message
    character(*), intent(in), optional :: file
    character(*), intent(in), optional :: line

    if (present(file) .and. present(line)) then
      write(6, *) '[Notice]: ' // trim(file) // ': ' // trim(line) // ': ' // trim(message)
    else
      write(6, *) '[Notice]: ' // trim(message)
    end if

  end subroutine log_notice

  subroutine log_warning(message, file, line)

    character(*), intent(in) :: message
    character(*), intent(in), optional :: file
    character(*), intent(in), optional :: line

    if (present(file) .and. present(line)) then
      write(6, *) '[Warning]: ' // trim(file) // ': ' // trim(line) // ': ' // trim(message)
    else
      write(6, *) '[Warning]: ' // trim(message)
    end if

  end subroutine log_warning

  subroutine log_error(message, file, line)

    character(*), intent(in) :: message
    character(*), intent(in), optional :: file
    character(*), intent(in), optional :: line

    if (present(file) .and. present(line)) then
      write(6, *) '[Error]: ' // trim(file) // ': ' // trim(line) // ': ' // trim(message)
    else
      write(6, *) '[Error]: ' // trim(message)
    end if
    stop 1

  end subroutine log_error

  subroutine log_step()

    type(map_iterator_type) iter
    class(*), pointer :: value

    write(6, '(" [Step]: ", A)', advance='no') trim(curr_time_format)

    iter = map_iterator_type(diags)
    do while (.not. iter%at_end())
      value => iter%value()
      select type(value)
      type is (real)
        write(6, '(X, A)', advance='no') trim(to_string(value, 20))
      end select
      call iter%next()
    end do
    write(6, *)

  end subroutine log_step

end module log_mod