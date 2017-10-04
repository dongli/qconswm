module dycore_mod

  use params_mod
  use mesh_mod
  use parallel_mod

  implicit none

  private

  public dycore_init
  public dycore_run
  public dycore_final
  public state

  type coef_type
    ! Coriolis coefficient at full/half meridional grids
    real, allocatable :: full_cori(:)
    real, allocatable :: half_cori(:)
    ! Curvature coefficient at full/half meridional grids 
    real, allocatable :: full_curv(:)
    real, allocatable :: half_curv(:)
    ! Zonal difference coefficient at full/half meridional grids
    real, allocatable :: full_dlon(:)
    real, allocatable :: half_dlon(:)
    ! Meridional difference coefficient at full/half meridional grids
    real, allocatable :: full_dlat(:)
    real, allocatable :: half_dlat(:)
  end type coef_type

  type state_type
    real, allocatable :: u(:,:)
    real, allocatable :: v(:,:)
    real, allocatable :: gd(:,:) ! Geopotential depth
  end type state_type

  type tend_type
    real, allocatable :: u_adv_lon(:,:)
    real, allocatable :: u_adv_lat(:,:)
    real, allocatable :: v_adv_lon(:,:)
    real, allocatable :: v_adv_lat(:,:)
    real, allocatable :: fu(:,:)
    real, allocatable :: fv(:,:)
    real, allocatable :: u_pgf(:,:)
    real, allocatable :: v_pgf(:,:)
    real, allocatable :: mass_div_lon(:,:)
    real, allocatable :: mass_div_lat(:,:)
    real, allocatable :: du(:,:)
    real, allocatable :: dv(:,:)
    real, allocatable :: dgd(:,:)
  end type tend_type

  type control_type
    ! 1: predict-correct
    ! 2: runge-kutta
    ! 3: leap-frog
    integer time_scheme
    integer time_step
    integer last_time_step
    integer old_time_idx
    integer new_time_idx
    real time_step_size
  end type control_type

  ! IAP transformed variables
  type iap_type
    real, allocatable :: u(:,:)
    real, allocatable :: v(:,:)
    real, allocatable :: gd(:,:)
  end type iap_type

  real, allocatable :: ghs(:,:) ! Surface geopotential

  type(coef_type) coef
  type(state_type) state(2)
  type(tend_type) tend(2)
  type(control_type) control
  type(iap_type) iap(2)

contains

  subroutine dycore_init()

    integer i, j, time_idx

    call mesh_init()
    call parallel_init()

    allocate(coef%full_cori(mesh%num_full_lat))
    allocate(coef%half_cori(mesh%num_half_lat))
    allocate(coef%full_curv(mesh%num_full_lat))
    allocate(coef%half_curv(mesh%num_half_lat))
    allocate(coef%full_dlon(mesh%num_full_lat))
    allocate(coef%half_dlon(mesh%num_half_lat))
    allocate(coef%full_dlat(mesh%num_full_lat))
    allocate(coef%half_dlat(mesh%num_half_lat))

    do j = 1, mesh%num_full_lat
      coef%full_cori(j) = 2.0 * omega * mesh%full_sin_lat(j)
      coef%full_curv(j) = mesh%full_sin_lat(j) / mesh%full_cos_lat(j) / radius
      coef%full_dlon(j) = 2.0 * radius * mesh%dlon * mesh%full_cos_lat(j)
      coef%full_dlat(j) = 2.0 * radius * mesh%dlat * mesh%full_cos_lat(j)
    end do

    do j = 1, mesh%num_half_lat
      coef%half_cori(j) = 2.0 * omega * mesh%half_sin_lat(j)
      coef%half_curv(j) = mesh%half_sin_lat(j) / mesh%half_cos_lat(j) / radius
      coef%half_dlon(j) = 2.0 * radius * mesh%dlon * mesh%half_cos_lat(j)
      coef%half_dlat(j) = 2.0 * radius * mesh%dlat * mesh%half_cos_lat(j)
    end do

    do time_idx = 1, 2
      call parallel_allocate(state(time_idx)%u, half_lon=.true.)
      call parallel_allocate(state(time_idx)%v, half_lat=.true.)
      call parallel_allocate(state(time_idx)%gd)
      call parallel_allocate(tend(time_idx)%u_adv_lon, half_lon=.true.)
      call parallel_allocate(tend(time_idx)%u_adv_lat, half_lon=.true.)
      call parallel_allocate(tend(time_idx)%fv, half_lon=.true.)
      call parallel_allocate(tend(time_idx)%u_pgf, half_lon=.true.)
      call parallel_allocate(tend(time_idx)%v_adv_lon, half_lat=.true.)
      call parallel_allocate(tend(time_idx)%v_adv_lat, half_lat=.true.)
      call parallel_allocate(tend(time_idx)%fu, half_lat=.true.)
      call parallel_allocate(tend(time_idx)%v_pgf, half_lat=.true.)
      call parallel_allocate(tend(time_idx)%mass_div_lon)
      call parallel_allocate(tend(time_idx)%mass_div_lat)
      call parallel_allocate(tend(time_idx)%du, half_lon=.true.)
      call parallel_allocate(tend(time_idx)%dv, half_lat=.true.)
      call parallel_allocate(tend(time_idx)%dgd)
      call parallel_allocate(iap(time_idx)%u, half_lon=.true.)
      call parallel_allocate(iap(time_idx)%v, half_lat=.true.)
      call parallel_allocate(iap(time_idx)%gd)
    end do

    call parallel_allocate(ghs)

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        ghs(i,j) = i + j
      end do
    end do
    call parallel_fill_halo(ghs, left_halo=.true., right_halo=.true., top_halo=.true., bottom_halo=.true.)
    do j = lbound(ghs, 2), ubound(ghs, 2)
      do i = lbound(ghs, 1), ubound(ghs, 1)
        if (i == ubound(ghs, 1)) then
          write(6, '(F5.1)') ghs(i,j)
        else
          write(6, '(F5.1)', advance='no') ghs(i,j)
        end if
      end do
    end do

    select case (time_scheme)
    case ('predict-correct')
      control%time_scheme = 1
    case ('runge-kutta')
      control%time_scheme = 2
    case ('leap-frog')
      control%time_scheme = 3
    case default
      write(6, *) '[Error]: Unknown time_scheme ' // trim(time_scheme) // '!'
      stop 1
    end select

    control%time_step = 0
    control%time_step_size = time_step_size
    control%old_time_idx = 1
    control%new_time_idx = 2

    write(6, *) '[Notice]: Dycore module is initialized.'

  end subroutine dycore_init

  subroutine dycore_run()

    do while (control%time_step /= control%last_time_step)
      call time_integrate()
      call advance_time()
    end do

  end subroutine dycore_run

  subroutine dycore_final()

    integer time_idx

    call mesh_final()
    call parallel_final()

    if (allocated(coef%full_cori)) deallocate(coef%full_cori)
    if (allocated(coef%half_cori)) deallocate(coef%half_cori)
    if (allocated(coef%full_curv)) deallocate(coef%full_curv)
    if (allocated(coef%half_curv)) deallocate(coef%half_curv)
    if (allocated(coef%full_dlon)) deallocate(coef%full_dlon)
    if (allocated(coef%half_dlon)) deallocate(coef%half_dlon)
    if (allocated(coef%full_dlat)) deallocate(coef%full_dlat)
    if (allocated(coef%half_dlat)) deallocate(coef%half_dlat)
    do time_idx = 1, 2
      if (allocated(state(time_idx)%u)) deallocate(state(time_idx)%u)
      if (allocated(state(time_idx)%v)) deallocate(state(time_idx)%v)
      if (allocated(state(time_idx)%gd)) deallocate(state(time_idx)%gd)
      if (allocated(tend(time_idx)%u_adv_lon)) deallocate(tend(time_idx)%u_adv_lon)
      if (allocated(tend(time_idx)%u_adv_lat)) deallocate(tend(time_idx)%u_adv_lat)
      if (allocated(tend(time_idx)%v_adv_lon)) deallocate(tend(time_idx)%v_adv_lon)
      if (allocated(tend(time_idx)%v_adv_lat)) deallocate(tend(time_idx)%v_adv_lat)
      if (allocated(tend(time_idx)%fu)) deallocate(tend(time_idx)%fu)
      if (allocated(tend(time_idx)%fv)) deallocate(tend(time_idx)%fv)
      if (allocated(tend(time_idx)%u_pgf)) deallocate(tend(time_idx)%u_pgf)
      if (allocated(tend(time_idx)%v_pgf)) deallocate(tend(time_idx)%v_pgf)
      if (allocated(tend(time_idx)%mass_div_lon)) deallocate(tend(time_idx)%mass_div_lon)
      if (allocated(tend(time_idx)%mass_div_lat)) deallocate(tend(time_idx)%mass_div_lat)
      if (allocated(tend(time_idx)%du)) deallocate(tend(time_idx)%du)
      if (allocated(tend(time_idx)%dv)) deallocate(tend(time_idx)%dv)
      if (allocated(tend(time_idx)%dgd)) deallocate(tend(time_idx)%dgd)
      if (allocated(iap(time_idx)%u)) deallocate(iap(time_idx)%u)
      if (allocated(iap(time_idx)%v)) deallocate(iap(time_idx)%v)
      if (allocated(iap(time_idx)%gd)) deallocate(iap(time_idx)%gd)
    end do
    if (allocated(ghs)) deallocate(ghs)

    write(6, *) '[Notice]: Dycore module is finalized.'

  end subroutine dycore_final

  subroutine advance_time()

    integer tmp

    tmp = control%old_time_idx
    control%old_time_idx = control%new_time_idx
    control%new_time_idx = tmp
    control%time_step = control%time_step + 1

  end subroutine advance_time

  subroutine iap_transform(state, iap)

    type(state_type), intent(in) :: state
    type(iap_type), intent(out) :: iap

    integer i, j

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        iap%gd(i,j) = sqrt(state%gd(i,j))
      end do
    end do

    ! TODO: Should we only fill right and top halo?
    call parallel_fill_halo(iap%gd(:,:))

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        iap%u(i,j) = 0.5 * (iap%gd(i,j) + iap%gd(i+1,j)) * state%u(i,j)
      end do
    end do

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        iap%v(i,j) = 0.5 * (iap%gd(i,j) + iap%gd(i,j+1)) * state%v(i,j)
      end do
    end do

  end subroutine iap_transform

  subroutine space_operators(state, iap, tend)

    integer i, j

    call calc_momentum_advection_operator(state, iap, tend)
    call calc_coriolis_and_curvature_operator(iap, tend)
    call calc_pressure_gradient_force_operator(state, iap, tend)
    call calc_mass_divergence_operator(iap, tend)

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%du(i,j) = - tend%u_adv_lon(i,j) - tend%u_adv_lat(i,j) + tend%fv(i,j) - tend%u_pgf(i,j)
      end do
    end do

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%dv(i,j) = - tend%v_adv_lon(i,j) - tend%v_adv_lat(i,j) - tend%fu(i,j) - tend%v_pgf(i,j)
      end do
    end do

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%dgd(i,j) = - tend%mass_div_lon(i,j) - tend%mass_div_lat(i,j)
      end do
    end do

  end subroutine space_operators

  subroutine calc_momentum_advection_operator(state, iap, tend)

    type(state_type), intent(in) :: state
    type(iap_type), intent(in) :: iap
    type(tend_type), intent(out) :: tend

    real up1, um1, vp1, vm1
    integer i, j

    ! U
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        up1 = state%u(i,j) + state%u(i+1,j)
        um1 = state%u(i,j) + state%u(i-1,j)
        tend%u_adv_lon(i,j) = 0.5 * coef%full_dlon(j) * (up1 * iap%u(i+1,j) - um1 * iap%u(i-1,j))

        vp1 = (state%v(i,j) + state%v(i+1,j)) * mesh%half_cos_lat(j)
        vm1 = (state%v(i,j-1) + state%v(i+1,j-1)) * mesh%half_cos_lat(j-1)
        tend%u_adv_lat(i,j) = 0.5 * coef%full_dlat(j) * (vp1 * iap%u(i,j+1) - vm1 * iap%ut(i,j-1))
      end do
    end do

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        up1 = state%u(i,j) + state%u(i,j+1)
        um1 = state%u(i-1,j) + state%u(i-1,j+1)
        tend%v_adv_lon(i,j) = 0.5 * coef%half_dlon(j) * (up1 * iap%v(i+1,j) - um1 * iap%vt(i-1,j))
      end do
    end do

    ! V
    do j = parallel%half_lat_start_idx_no_pole, parallel%half_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        vp1 = state%v(i,j) * mesh%half_cos_lat(j) + state%v(i,j+1) * mesh%half_cos_lat(j+1)
        vm1 = state%v(i,j) * mesh%half_cos_lat(j) + state%v(i,j-1) * mesh%half_cos_lat(j-1)
        tend%v_adv_lat(i,j) = 0.5 * coef%half_dlat(j) * (vp1 * iap%v(i,j+1) - vm1 * iap%v(i,j-1))
      end do
    end do

    ! Handle meridional advection at North Pole.
    if (parallel%half_lat_start_idx == parallel%half_lat_south_pole_idx) then
      j = parallel%half_lat_south_pole_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        vp1 = state%v(i,j) * mesh%half_cos_lat(j) + state%v(i,j+1) * mesh%half_cos_lat(j+1)
        tend%v_adv_lat(i,j) = 0.5 * coef%half_dlat(j) * vp1 * iap%v(i,j+1)
      end do
    end if

    ! Handle meridional advection at South Pole.
    if (parallel%half_lat_end_idx == parallel%half_lat_north_pole_idx) then
      j = parallel%half_lat_north_pole_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        vm1 = state%v(i,j) * mesh%half_cos_lat(j) + state%v(i,j-1) * mesh%half_cos_lat(j-1)
        tend%v_adv_lat(i,j) = - 0.5 * coef%half_dlat(j) * vm1 * iap%v(i,j-1)
      end do
    end if

  end subroutine calc_momentum_advection_operator

  subroutine calc_coriolis_and_curvature_operator(iap, tend)

    type(iap_type), intent(in) :: iap
    type(tend_type), intent(out) :: tend

    integer i, j

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%fv(i,j) = 0.25 * (coef%half_cori(j  ) * (iap%v(i,j  ) + iap%v(i+1,j  )) + &
                                coef%half_cori(j-1) * (iap%v(i,j-1) + iap%v(i+1,j-1)))
      end do
    end do

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%fu(i,j) = 0.25 * (coef%full_cori(j  ) * (iap%u(i,j  ) + iap%u(i-1,j  )) + &
                                coef%full_cori(j+1) * (iap%u(i,j+1) + iap%u(i-1,j+1)))
      end do
    end do

  end subroutine calc_coriolis_and_curvature_operator

  subroutine calc_pressure_gradient_force_operator(state, iap, tend)

    type(state_type), intent(in) :: state
    type(iap_type), intent(in) :: iap
    type(tend_type), intent(out) :: tend

    integer i, j

    ! U
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%u_pgf(i,j) = coef%full_dlon(j) * (iap%gd(i,j) + iap%gd(i+1,j)) * &
          (state%gd(i+1,j) - state%gd(i,j) + ghs(i+1,j) - ghs(i,j))
      end do
    end do

    ! V
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%half_lon_end_idx
        tend%v_pgf(i,j) = coef%half_dlat(j) * (iap%gd(i,j) + iap%gd(i,j+1)) * &
          (state%gd(i,j+1) - state%gd(i,j) + ghs(i,j+1) - ghs(i,j)) * mesh%half_cos_lat(j)
      end do
    end do

  end subroutine calc_pressure_gradient_force_operator

  subroutine calc_mass_divergence_operator(iap, tend)

    type(iap_type), intent(in) :: iap
    type(tend_type), intent(out) :: tend

    integer i, j
    real sp, np, sum

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%mass_div_lon(i,j) = coef%full_dlon(j) * ( &
          (iap%gd(i,j) + iap%gd(i+1,j)) * iap%u(i,j) - &
          (iap%gd(i,j) + iap%gd(i-1,j)) * iap%u(i-1,j))
        tend%mass_div_lat(i,j) = coef%full_dlat(j) * ( &
          (iap%gd(i,j) + iap%gd(i,j+1)) * iap%v(i,j) * mesh%half_cos_lat(j) - &
          (iap%gd(i,j) + iap%gd(i,j-1)) * iap%v(i,j-1) * mesh%half_cos_lat(j-1))
      end do
    end do

    if (parallel%full_lat_start_idx == parallel%full_lat_north_pole_idx) then
      j = parallel%full_lat_north_pole_idx
      sp = 0.0
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        sp = sp + (iap%gd(i,j) + iap%gd(i,j+1)) * iap%v(i,j) * mesh%half_cos_lat(j)
      end do
      call parallel_zonal_sum(sp, sum)
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%mass_div_lat(i,j) = sum
      end do
    end if

    if (parallel%full_lat_end_idx == parallel%full_lat_south_pole_idx) then
      j = parallel%full_lat_south_pole_idx
      np = 0.0
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        np = np - (iap%gd(i,j) + iap%gd(i,j-1)) * iap%v(i,j-1) * mesh%half_cos_lat(j-1)
      end do
      call parallel_zonal_sum(np, sum)
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%mass_div_lat(i,j) = sum
      end do
    end if

  end subroutine calc_mass_divergence_operator

  subroutine update_state(dt, tend, old_state, new_state)

    real, intent(in) :: dt
    type(tend_type), intent(in) :: tend
    type(state_type), intent(out) :: state

    integer i, j

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        new_state%u(i,j) = old_state%u(i,j) + dt * tend%du(i,j)
      end do
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        new_state%v(i,j) = old_state%v(i,j) + dt * tend%dv(i,j)
      end do
    end do
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        new_state%gd(i,j) = old_state%gd(i,j) + dt * tend%dgd(i,j)
      end do
    end do

  end subroutine update_state

  real function inner_product(tend1, tend2) result(res)

    type(tend_type), intent(in) :: tend1
    type(tend_type), intent(in) :: tend2

    integer i, j

    res = 0.0

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        res = res + dtend1%du(i,j) * dtend2%du(i,j)
      end do
    end do

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        res = res + dtend1%dv(i,j) * dtend2%dv(i,j)
      end do
    end do

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        res = res + dtend1%dgd(i,j) * dtend2%dgd(i,j)
      end do
    end do

  end function inner_product

  subroutine time_integrate()

    select case (control%time_scheme)
    case (1)
      call predict_correct()
    case (2)
      call runge_kutta()
    case (3)
      call leap_frog()
    end select

  end subroutine time_integrate

  subroutine predict_correct()

    integer old, new
    real dt, dt05

    old = control%old_time_idx
    new = control%new_time_idx
    dt05 = control%time_step_size * 0.5

    ! Do first predict step.
    call iap_transform(state(old), iap(old))
    call space_operators(state(old), iap(old), tend(old))
    call update_state(dt05, tend(old), state(old), state(new))

    ! Do second predict step.
    call iap_transform(state(new), iap(new))
    call space_operators(state(new), iap(new), tend(old))
    call update_state(dt05, tend(old), state(old), state(new))

    ! Do correct step.
    call iap_transform(state(new), iap(new))
    call space_operators(state(new), iap(new), tend(new))
    ! Calculate modified time step size.
    dt = control%time_step_size * inner_product(tend(old), tend(new)) / inner_product(tend(new), tend(new))
    call update_state(dt, tend(new), state(old), state(new))

  end subroutine predict_correct

  subroutine runge_kutta()

  end subroutine runge_kutta

  subroutine leap_frog()

  end subroutine leap_frog

end module dycore_mod
