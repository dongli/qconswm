module time_mod

  use datetime_mod
  use timedelta_mod
  use map_mod
  use params_mod, time_step_size_in => time_step_size, time_units_in => time_units

  implicit none

  private

  public time_init
  public time_advance
  public time_ended
  public time_add_alert
  public time_is_alerted

  public curr_time_format
  public time_step
  public old_time_idx
  public new_time_idx
  public time_units

  type alert_type
    type(timedelta_type) period
    type(datetime_type) last_time
  end type alert_type

  type(datetime_type) start_time
  type(datetime_type) end_time
  type(datetime_type) curr_time
  type(timedelta_type) time_step_size
  type(map_type) alerts
  integer time_step
  integer old_time_idx
  integer new_time_idx
  character(30) time_units
  character(30) curr_time_format

contains

  subroutine time_init()

    start_time = datetime(year=year_range(1), month=month_range(1), day=day_range(1), &
      hour=hour_range(1), minute=minute_range(1), second=second_range(1))
    end_time = datetime(year=year_range(2), month=month_range(2), day=day_range(2), &
      hour=hour_range(2), minute=minute_range(2), second=second_range(2))

    time_step = 0
    old_time_idx = 1
    new_time_idx = 2
    time_step_size = timedelta(seconds=time_step_size_in)
    time_units = time_units_in

    curr_time = start_time
    curr_time_format = curr_time%isoformat()

  end subroutine time_init

  subroutine time_advance()

    integer tmp

    tmp = old_time_idx
    old_time_idx = new_time_idx
    new_time_idx = tmp
    time_step = time_step + 1

    curr_time = curr_time + time_step_size
    curr_time_format = curr_time%isoformat()

  end subroutine time_advance

  logical function time_ended() result(res)

    res = curr_time > end_time

  end function time_ended

  subroutine time_add_alert(name, days, hours, minutes)

    character(*), intent(in) :: name
    class(*), intent(in), optional :: days
    class(*), intent(in), optional :: hours
    class(*), intent(in), optional :: minutes

    type(alert_type) alert

    alert%period = timedelta(days, hours, minutes)
    alert%last_time = start_time
    call alerts%insert(name, alert)

  end subroutine time_add_alert

  function time_is_alerted(name) result(res)

    character(*), intent(in) :: name
    logical res

    type(alert_type), pointer :: alert => null()
    type(datetime_type) time

    alert => get_alert(name)
    if (associated(alert)) then
      time = alert%last_time + alert%period
      if (time <= curr_time) then
        alert%last_time = curr_time
        res = .true.
      else
        res = .false.
      end if
    else
      res = .false.
    end if

  end function time_is_alerted

  function get_alert(name) result(res)

    character(*), intent(in) :: name
    type(alert_type), pointer :: res

    class(*), pointer :: value

    value => alerts%value(name)
    select type (value)
    type is (alert_type)
      res => value
    end select

  end function get_alert

end module time_mod