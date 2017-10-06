program qconswm

  use params_mod
  use io_mod
  use dycore_mod
  use rossby_haurwitz_test_mod

  character(256) namelist_file_path

  if (command_argument_count() /= 1) then
    write(6, *) 'Usage: ./qconswm <namelist_file_path>'
    stop 1
  end if

  call get_command_argument(1, namelist_file_path)

  call params_read(namelist_file_path)

  call io_create_dataset(desc='Rossby-Haurwitz test', author='Li Dong <dongli@lasg.iap.ac.cn>', file_prefix='rh_test')

  call dycore_init()

  call rossby_haurwitz_test_set_initial_condition()

  call dycore_run()

  call dycore_final()

end program qconswm
