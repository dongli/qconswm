module params_mod

  implicit none

  ! Contant parameters
  real, parameter :: pi = atan(1.0) * 4.0
  real, parameter :: rad_to_deg = 180.0 / pi
  real, parameter :: deg_to_rad = pi / 180.0
  real, parameter :: omega = 2.0 * pi / 86400.0
  real, parameter :: radius = 6.371e6

  integer num_lon
  integer num_lat

  ! Options:
  ! - predict-correct
  ! - runge-kutta
  character(30) time_scheme ! Time integration scheme
  integer time_order ! Time integration order (different schemes will have different meanings)
  logical qcon_modified ! Switch whether quadratic conservation modification is added

  namelist /qconswm_params/ &
    num_lon, &
    num_lat, &
    time_scheme, &
    time_order, &
    qcon_modified

contains

  subroutine params_read(file_path)

    character(*), intent(in) :: file_path

    open(10, file=file_path)
    read(10, nml=qconswm_params)
    close(10)

  end subroutine params_read

end module params_mod
