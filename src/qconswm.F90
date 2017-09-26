program qconswm

  use params_mod
  use dycore_mod

  character(256) namelist_file_path

  if (command_argument_count() /= 1) then
    write(6, *) 'Usage: ./qconswm <namelist_file_path>'
    stop 1
  end if

  call get_command_argument(1, namelist_file_path)

  call params_read(namelist_file_path)

  call dycore_init()

  call dycore_final()

end program qconswm
