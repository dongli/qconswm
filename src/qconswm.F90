program qconswm

  use params_mod
  use io_mod
  use time_mod
  use dycore_mod
  use rossby_haurwitz_test_mod

  character(256) namelist_file_path

  if (command_argument_count() /= 1) then
    write(6, *) 'Usage: ./qconswm <namelist_file_path>'
    stop 1
  end if

  call get_command_argument(1, namelist_file_path)

  call params_read(namelist_file_path)

  call dycore_init()

  if (is_restart_run) then
    call dycore_restart()
  else
    call rossby_haurwitz_test_set_initial_condition()
  end if

  call dycore_run()

  call dycore_final()

end program qconswm
