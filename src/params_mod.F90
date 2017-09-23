module params_mod

  implicit none

  integer num_lon
  integer num_lat

  namelist /qconswm_params/ &
    num_lon, &
    num_lat

contains

  subroutine params_read(file_path)

    character(*), intent(in) :: file_path

    open(10, file=file_path)
    read(10, nml=qconswm_params)
    close(10)

  end subroutine params_read

end module params_mod
