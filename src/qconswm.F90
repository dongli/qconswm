program qconswm

  use params_mod

  character(256) namelist_file_path

  call get_command_argument(1, namelist_file_path)

  call params_read(namelist_file_path)

end program qconswm
