module params_mod

  implicit none

  ! Contant parameters
  real, parameter :: pi = atan(1.0) * 4.0
  real, parameter :: rad_to_deg = 180.0 / pi
  real, parameter :: deg_to_rad = pi / 180.0
  real, parameter :: omega = 2.0 * pi / 86400.0
  real, parameter :: radius = 6.371e6
  real, parameter :: g = 9.8

  integer num_lon
  integer num_lat
  real(8) time_step_size
 
  integer :: year_range(2) = [0, 0]
  integer :: month_range(2) = [1, 1]
  integer :: day_range(2) = [1, 1]
  integer :: hour_range(2) = [0, 0]
  integer :: minute_range(2) = [0, 0]
  integer :: second_range(2) = [0, 0]
  character(30) time_units

  ! Options:
  ! - predict-correct
  ! - runge-kutta
  character(30) time_scheme ! Time integration scheme
  integer time_order ! Time integration order (different schemes will have different meanings)
  logical qcon_modified ! Switch whether quadratic conservation modification is added

  namelist /qconswm_params/ &
    num_lon, &
    num_lat, &
    year_range, &
    month_range, &
    day_range, &
    hour_range, &
    month_range, &
    second_range, &
    time_units, &
    time_step_size, &
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
