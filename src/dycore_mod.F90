module dycore_mod

  use log_mod
  use params_mod, time_scheme_in => time_scheme, split_scheme_in => split_scheme
  use mesh_mod
  use time_mod
  use parallel_mod
  use io_mod
  use types_mod

  implicit none

  private

  public dycore_init
  public dycore_run
  public dycore_final
  public state
  public static

  ! 1: predict-correct
  ! 2: runge-kutta
  ! 3: leap-frog
  integer time_scheme
  ! 1: csp1: first order conservative split
  ! 2: csp2: second order conservative split
  ! 3: isp: inproved second order conservative split
  integer split_scheme
  integer, parameter :: all_pass = 0
  integer, parameter :: fast_pass = 1
  integer, parameter :: slow_pass = 2

  type(coef_type) coef
  type(state_type) state(2)
  type(static_type) static
  type(tend_type) tend(2)
  type(iap_type) iap(2)

contains

  subroutine dycore_init()

    integer i, j, time_idx

    call mesh_init()
    call time_init()
    call parallel_init()
    call io_init()

    allocate(coef%cori(mesh%num_full_lat))
    allocate(coef%curv(mesh%num_full_lat))
    allocate(coef%full_dlon(mesh%num_full_lat))
    allocate(coef%half_dlon(mesh%num_half_lat))
    allocate(coef%full_dlat(mesh%num_full_lat))
    allocate(coef%half_dlat(mesh%num_half_lat))

    do j = 1, mesh%num_full_lat
      coef%cori(j) = 2.0 * omega * mesh%full_sin_lat(j)
      if (j == 1 .or. j == mesh%num_full_lat) then
        coef%curv(j) = 0.0
      else
        coef%curv(j) = mesh%full_sin_lat(j) / mesh%full_cos_lat(j) / radius
      end if
      coef%full_dlon(j) = 2.0 * radius * mesh%dlon * mesh%full_cos_lat(j)
      coef%full_dlat(j) = 2.0 * radius * mesh%dlat * mesh%full_cos_lat(j)
    end do

    do j = 1, mesh%num_half_lat
      coef%half_dlon(j) = 2.0 * radius * mesh%dlon * mesh%half_cos_lat(j)
      coef%half_dlat(j) = 2.0 * radius * mesh%dlat * mesh%half_cos_lat(j)
    end do

    do time_idx = 1, 2
      call parallel_allocate(state(time_idx)%u, half_lon=.true.)
      call parallel_allocate(state(time_idx)%v, half_lat=.true.)
      call parallel_allocate(state(time_idx)%gd)
      call parallel_allocate(state(time_idx)%ua)
      call parallel_allocate(state(time_idx)%va)
      call parallel_allocate(tend(time_idx)%u_adv_lon, half_lon=.true.)
      call parallel_allocate(tend(time_idx)%u_adv_lat, half_lon=.true.)
      call parallel_allocate(tend(time_idx)%fv, half_lon=.true.)
      call parallel_allocate(tend(time_idx)%cv, half_lon=.true.)
      call parallel_allocate(tend(time_idx)%u_pgf, half_lon=.true.)
      call parallel_allocate(tend(time_idx)%v_adv_lon, half_lat=.true.)
      call parallel_allocate(tend(time_idx)%v_adv_lat, half_lat=.true.)
      call parallel_allocate(tend(time_idx)%fu, half_lat=.true.)
      call parallel_allocate(tend(time_idx)%cu, half_lat=.true.)
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
    call parallel_allocate(static%ghs)

    select case (time_scheme_in)
    case ('predict-correct')
      time_scheme = 1
    case ('runge-kutta')
      time_scheme = 2
    case ('leap-frog')
      time_scheme = 3
    case default
      call log_error('Unknown time_scheme ' // trim(time_scheme_in) // '!')
    end select

    select case (split_scheme_in)
    case ('csp1')
      split_scheme = 1
    case ('csp2')
      split_scheme = 2
    case ('isp')
      split_scheme = 3
    case default
      split_scheme = 0
      call log_notice('No fast-slow split.')
    end select

    call io_add_dim('lon', size=mesh%num_full_lon)
    call io_add_dim('lat', size=mesh%num_full_lat)
    call io_add_dim('time')

    call io_add_var('u', long_name='u wind component', units='m s-1', dim_names=['lon ', 'lat ', 'time'])
    call io_add_var('v', long_name='v wind component', units='m s-1', dim_names=['lon ', 'lat ', 'time'])
    call io_add_var('gd', long_name='geopotential depth', units='m2 s-2', dim_names=['lon ', 'lat ', 'time'])

    call log_notice('Dycore module is initialized.')

  end subroutine dycore_init

  subroutine dycore_run()

    call reset_cos_lat_at_poles()

    call iap_transform(state(old_time_idx), iap(old_time_idx))

    call output(state(old_time_idx))
    call log_add_diag('total_mass', total_mass(state(old_time_idx)))
    call log_add_diag('total_energy', total_energy(iap(old_time_idx)))
    call log_step()

    do while (.not. time_ended())
      call time_integrate()
      call time_advance()
      call output(state(old_time_idx))
      call log_add_diag('total_mass', total_mass(state(old_time_idx)))
      call log_add_diag('total_energy', total_energy(iap(old_time_idx)))
      call log_step()
    end do

  end subroutine dycore_run

  subroutine dycore_final()

    integer time_idx

    call mesh_final()
    call parallel_final()

    if (allocated(coef%cori)) deallocate(coef%cori)
    if (allocated(coef%curv)) deallocate(coef%curv)
    if (allocated(coef%full_dlon)) deallocate(coef%full_dlon)
    if (allocated(coef%half_dlon)) deallocate(coef%half_dlon)
    if (allocated(coef%full_dlat)) deallocate(coef%full_dlat)
    if (allocated(coef%half_dlat)) deallocate(coef%half_dlat)
    do time_idx = 1, 2
      if (allocated(state(time_idx)%u)) deallocate(state(time_idx)%u)
      if (allocated(state(time_idx)%v)) deallocate(state(time_idx)%v)
      if (allocated(state(time_idx)%gd)) deallocate(state(time_idx)%gd)
      if (allocated(state(time_idx)%ua)) deallocate(state(time_idx)%ua)
      if (allocated(state(time_idx)%va)) deallocate(state(time_idx)%va)
      if (allocated(tend(time_idx)%u_adv_lon)) deallocate(tend(time_idx)%u_adv_lon)
      if (allocated(tend(time_idx)%u_adv_lat)) deallocate(tend(time_idx)%u_adv_lat)
      if (allocated(tend(time_idx)%v_adv_lon)) deallocate(tend(time_idx)%v_adv_lon)
      if (allocated(tend(time_idx)%v_adv_lat)) deallocate(tend(time_idx)%v_adv_lat)
      if (allocated(tend(time_idx)%fu)) deallocate(tend(time_idx)%fu)
      if (allocated(tend(time_idx)%fv)) deallocate(tend(time_idx)%fv)
      if (allocated(tend(time_idx)%cu)) deallocate(tend(time_idx)%cu)
      if (allocated(tend(time_idx)%cv)) deallocate(tend(time_idx)%cv)
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
    if (allocated(static%ghs)) deallocate(static%ghs)

    call log_notice('Dycore module is finalized.')

  end subroutine dycore_final

  subroutine reset_cos_lat_at_poles()

    integer j

    j = parallel%full_lat_south_pole_idx
    mesh%full_cos_lat(j) = mesh%half_cos_lat(parallel%half_lat_south_pole_idx) * 0.25
    coef%full_dlon(j) = 2.0 * radius * mesh%dlon * mesh%full_cos_lat(j)
    coef%full_dlat(j) = 2.0 * radius * mesh%dlat * mesh%full_cos_lat(j)

    j = parallel%full_lat_north_pole_idx
    mesh%full_cos_lat(j) = mesh%half_cos_lat(parallel%half_lat_north_pole_idx) * 0.25
    coef%full_dlon(j) = 2.0 * radius * mesh%dlon * mesh%full_cos_lat(j)
    coef%full_dlat(j) = 2.0 * radius * mesh%dlat * mesh%full_cos_lat(j)

  end subroutine reset_cos_lat_at_poles

  subroutine output(state)

    type(state_type), intent(inout) :: state

    integer i, j

    if (.not. time_is_alerted('output #1')) return

    ! Convert wind from C grid to A grid.
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        state%ua(i,j) = 0.5 * (state%u(i,j) + state%u(i+1,j))
      end do
    end do

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%va(i,j) = 0.5 * (state%v(i,j) + state%v(i,j-1))
      end do
    end do

    call io_start_output()
    call io_output('lon', mesh%lon_deg(:))
    call io_output('lat', mesh%lat_deg(:))
    call io_output('u', state%ua(parallel%full_lon_start_idx:parallel%full_lon_end_idx, &
                                 parallel%full_lat_start_idx:parallel%full_lat_end_idx))
    call io_output('v', state%va(parallel%full_lon_start_idx:parallel%full_lon_end_idx, &
                                 parallel%full_lat_start_idx:parallel%full_lat_end_idx))
    call io_output('gd', state%gd(parallel%full_lon_start_idx:parallel%full_lon_end_idx, &
                                  parallel%full_lat_start_idx:parallel%full_lat_end_idx))
    call io_end_output()

  end subroutine output

  subroutine iap_transform(state, iap)

    type(state_type), intent(in) :: state
    type(iap_type), intent(inout) :: iap

    integer i, j

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        iap%gd(i,j) = sqrt(state%gd(i,j))
      end do
    end do

    call parallel_fill_halo(iap%gd(:,:), all_halo=.true.)

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

    call parallel_fill_halo(iap%u(:,:), all_halo=.true.)
    call parallel_fill_halo(iap%v(:,:), all_halo=.true.)

  end subroutine iap_transform

  subroutine space_operators(state, iap, tend)

    type(state_type), intent(in) :: state
    type(iap_type), intent(in) :: iap
    type(tend_type), intent(inout) :: tend

    integer i, j

    call momentum_advection_operator(state, iap, tend)
    call coriolis_operator(iap, tend)
    call curvature_operator(state, iap, tend)
    call pressure_gradient_force_operator(state, iap, tend)
    call mass_divergence_operator(iap, tend)

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%du(i,j) = - tend%u_adv_lon(i,j) - tend%u_adv_lat(i,j) + tend%fv(i,j) + tend%cv(i,j) - tend%u_pgf(i,j)
      end do
    end do

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%dv(i,j) = - tend%v_adv_lon(i,j) - tend%v_adv_lat(i,j) - tend%fu(i,j) - tend%cu(i,j) - tend%v_pgf(i,j)
      end do
    end do

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%dgd(i,j) = - tend%mass_div_lon(i,j) - tend%mass_div_lat(i,j)
      end do
    end do

  end subroutine space_operators

  subroutine momentum_advection_operator(state, iap, tend)

    type(state_type), intent(in) :: state
    type(iap_type), intent(in) :: iap
    type(tend_type), intent(inout) :: tend

    real up1, um1, vp1, vm1
    integer i, j

    ! U
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        up1 = state%u(i,j) + state%u(i+1,j)
        um1 = state%u(i,j) + state%u(i-1,j)
        tend%u_adv_lon(i,j) = 0.5 / coef%full_dlon(j) * (up1 * iap%u(i+1,j) - um1 * iap%u(i-1,j))

        vp1 = (state%v(i,j  ) + state%v(i+1,j  )) * mesh%half_cos_lat(j  )
        vm1 = (state%v(i,j-1) + state%v(i+1,j-1)) * mesh%half_cos_lat(j-1)
        tend%u_adv_lat(i,j) = 0.5 / coef%full_dlat(j) * (vp1 * iap%u(i,j+1) - vm1 * iap%u(i,j-1))
      end do
    end do

    ! V
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        up1 = state%u(i,  j) + state%u(i,  j+1)
        um1 = state%u(i-1,j) + state%u(i-1,j+1)
        tend%v_adv_lon(i,j) = 0.5 / coef%half_dlon(j) * (up1 * iap%v(i+1,j) - um1 * iap%v(i-1,j))
      end do
    end do

    do j = parallel%half_lat_start_idx_no_pole, parallel%half_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        vp1 = state%v(i,j) * mesh%half_cos_lat(j) + state%v(i,j+1) * mesh%half_cos_lat(j+1)
        vm1 = state%v(i,j) * mesh%half_cos_lat(j) + state%v(i,j-1) * mesh%half_cos_lat(j-1)
        tend%v_adv_lat(i,j) = 0.5 / coef%half_dlat(j) * (vp1 * iap%v(i,j+1) - vm1 * iap%v(i,j-1))
      end do
    end do

    ! Handle meridional advection at South Pole.
    if (parallel%has_south_pole) then
      j = parallel%half_lat_south_pole_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        vp1 = state%v(i,j) * mesh%half_cos_lat(j) + state%v(i,j+1) * mesh%half_cos_lat(j+1)
        tend%v_adv_lat(i,j) = 0.5 / coef%half_dlat(j) * vp1 * iap%v(i,j+1)
      end do
    end if

    ! Handle meridional advection at North Pole.
    if (parallel%has_north_pole) then
      j = parallel%half_lat_north_pole_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        vm1 = state%v(i,j) * mesh%half_cos_lat(j) + state%v(i,j-1) * mesh%half_cos_lat(j-1)
        tend%v_adv_lat(i,j) = - 0.5 / coef%half_dlat(j) * vm1 * iap%v(i,j-1)
      end do
    end if

  end subroutine momentum_advection_operator

  subroutine coriolis_operator(iap, tend)

    type(iap_type), intent(in) :: iap
    type(tend_type), intent(inout) :: tend

    real c1, c2
    integer i, j

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      c1 = mesh%half_cos_lat(j  ) / mesh%full_cos_lat(j)
      c2 = mesh%half_cos_lat(j-1) / mesh%full_cos_lat(j)
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%fv(i,j) = 0.25 * coef%cori(j) * (c1 * (iap%v(i,j  ) + iap%v(i+1,j  )) + &
                                              c2 * (iap%v(i,j-1) + iap%v(i+1,j-1)))
      end do
    end do

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%fu(i,j) = 0.25 * (coef%cori(j  ) * (iap%u(i,j  ) + iap%u(i-1,j  )) + &
                               coef%cori(j+1) * (iap%u(i,j+1) + iap%u(i-1,j+1)))
      end do
    end do

#ifndef NDEBUG
    if (parallel%has_south_pole) then
      j = parallel%full_lat_south_pole_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        if (iap%u(i,j) /= 0.0) then
          call log_error('U at South Pole should be zero!')
        end if
      end do
    end if
    if (parallel%has_north_pole) then
      j = parallel%full_lat_north_pole_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        if (iap%u(i,j) /= 0.0) then
          call log_error('U at North Pole should be zero!')
        end if
      end do
    end if
#endif

  end subroutine coriolis_operator

  subroutine curvature_operator(state, iap, tend)

    type(state_type), intent(in) :: state
    type(iap_type), intent(in) :: iap
    type(tend_type), intent(inout) :: tend

    real c1, c2
    integer i, j

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      c1 = mesh%half_cos_lat(j  ) / mesh%full_cos_lat(j)
      c2 = mesh%half_cos_lat(j-1) / mesh%full_cos_lat(j)
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%cv(i,j) = 0.25 * coef%curv(j) * state%u(i,j) * &
          (c1 * (iap%v(i,j  ) + iap%v(i+1,j  )) + &
           c2 * (iap%v(i,j-1) + iap%v(i+1,j-1)))
      end do
    end do

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%cu(i,j) = 0.25 * &
          (coef%curv(j  ) * state%u(i,  j  ) * iap%u(i,  j  ) + &
           coef%curv(j+1) * state%u(i,  j+1) * iap%u(i,  j+1) + &
           coef%curv(j  ) * state%u(i-1,j  ) * iap%u(i-1,j  ) + &
           coef%curv(j+1) * state%u(i-1,j+1) * iap%u(i-1,j+1))
      end do
    end do

  end subroutine curvature_operator

  subroutine pressure_gradient_force_operator(state, iap, tend)

    type(state_type), intent(in) :: state
    type(iap_type), intent(in) :: iap
    type(tend_type), intent(inout) :: tend

    integer i, j

    ! U
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%u_pgf(i,j) = (iap%gd(i,j) + iap%gd(i+1,j)) / coef%full_dlon(j) * &
          (state%gd(i+1,j) - state%gd(i,j) + static%ghs(i+1,j) - static%ghs(i,j))
      end do
    end do

    ! V
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%v_pgf(i,j) = (iap%gd(i,j) + iap%gd(i,j+1)) / coef%half_dlat(j) * &
          (state%gd(i,j+1) - state%gd(i,j) + static%ghs(i,j+1) - static%ghs(i,j)) * mesh%half_cos_lat(j)
      end do
    end do

  end subroutine pressure_gradient_force_operator

  subroutine mass_divergence_operator(iap, tend)

    type(iap_type), intent(in) :: iap
    type(tend_type), intent(inout) :: tend

    integer i, j
    real sp, np, sum

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%mass_div_lon(i,j) = 1.0 / coef%full_dlon(j) * ( &
          (iap%gd(i,j) + iap%gd(i+1,j)) * iap%u(i,  j) - &
          (iap%gd(i,j) + iap%gd(i-1,j)) * iap%u(i-1,j))
        tend%mass_div_lat(i,j) = 1.0 / coef%full_dlat(j) * ( &
          (iap%gd(i,j) + iap%gd(i,j+1)) * iap%v(i,j  ) * mesh%half_cos_lat(j  ) - &
          (iap%gd(i,j) + iap%gd(i,j-1)) * iap%v(i,j-1) * mesh%half_cos_lat(j-1))
      end do
    end do

    if (parallel%full_lat_start_idx == parallel%full_lat_south_pole_idx) then
      j = parallel%full_lat_south_pole_idx
      sp = 0.0
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        sp = sp + (iap%gd(i,j) + iap%gd(i,j+1)) * iap%v(i,j) * mesh%half_cos_lat(j)
      end do
      call parallel_zonal_sum(sp, sum)
      sum = sum / mesh%num_full_lon / coef%full_dlat(j)
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%mass_div_lat(i,j) = sum
      end do
    end if

    if (parallel%full_lat_end_idx == parallel%full_lat_north_pole_idx) then
      j = parallel%full_lat_north_pole_idx
      np = 0.0
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        np = np - (iap%gd(i,j) + iap%gd(i,j-1)) * iap%v(i,j-1) * mesh%half_cos_lat(j-1)
      end do
      call parallel_zonal_sum(np, sum)
      sum = sum / mesh%num_full_lon / coef%full_dlat(j)
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%mass_div_lat(i,j) = sum
      end do
    end if

  end subroutine mass_divergence_operator

  subroutine update_state(dt, tend, old_state, old_iap, new_state, new_iap)

    real, intent(in) :: dt
    type(tend_type), intent(in) :: tend
    type(state_type), intent(in) :: old_state
    type(iap_type), intent(in) :: old_iap
    type(state_type), intent(inout) :: new_state
    type(iap_type), intent(inout) :: new_iap

    integer i, j

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        new_state%gd(i,j) = old_state%gd(i,j) + dt * tend%dgd(i,j)
        new_iap%gd(i,j) = sqrt(new_state%gd(i,j))
      end do
    end do

    call parallel_fill_halo(new_iap%gd(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%gd(:,:), all_halo=.true.)

    ! Update IAP wind state.
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        new_iap%u(i,j) = old_iap%u(i,j) + dt * tend%du(i,j)
      end do
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        new_iap%v(i,j) = old_iap%v(i,j) + dt * tend%dv(i,j)
      end do
    end do

    ! Transform from IAP to normal state.
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        new_state%u(i,j) = new_iap%u(i,j) * 2.0 / (new_iap%gd(i,j) + new_iap%gd(i+1,j))
      end do
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        new_state%v(i,j) = new_iap%v(i,j) * 2.0 / (new_iap%gd(i,j) + new_iap%gd(i,j+1))
      end do
    end do

    call parallel_fill_halo(new_iap%u(:,:), all_halo=.true.)
    call parallel_fill_halo(new_iap%v(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%u(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%v(:,:), all_halo=.true.)

  end subroutine update_state

  real function inner_product(tend1, tend2) result(res)

    type(tend_type), intent(in) :: tend1
    type(tend_type), intent(in) :: tend2

    integer i, j

    res = 0.0

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        res = res + tend1%du(i,j) * tend2%du(i,j) * mesh%full_cos_lat(j)
      end do
    end do

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        res = res + tend1%dv(i,j) * tend2%dv(i,j) * mesh%half_cos_lat(j)
      end do
    end do

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        res = res + tend1%dgd(i,j) * tend2%dgd(i,j) * mesh%full_cos_lat(j)
      end do
    end do

  end function inner_product

  real function total_mass(state)

    type(state_type), intent(in) :: state

    integer i, j

    total_mass = 0.0
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        total_mass = total_mass + mesh%full_cos_lat(j) * mesh%dlon * mesh%dlat * state%gd(i,j)
      end do
    end do
    total_mass = total_mass * radius**2

  end function total_mass

  real function total_energy(iap)

    type(iap_type), intent(in) :: iap

    integer i, j

    total_energy = 0.0
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        total_energy = total_energy + iap%u(i,j)**2 * mesh%full_cos_lat(j)
      end do
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        total_energy = total_energy + iap%v(i,j)**2 * mesh%half_cos_lat(j)
      end do
    end do
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        total_energy = total_energy + (iap%gd(i,j)**2 + static%ghs(i,j))**2 * mesh%full_cos_lat(j)
      end do
    end do

  end function total_energy

  subroutine time_integrate()

    select case (time_scheme)
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
    real dt, ip1, ip2

    old = old_time_idx
    new = new_time_idx
    dt = time_step_size * 0.5

    ! Do first predict step.
    call space_operators(state(old), iap(old), tend(old))
    call update_state(dt, tend(old), state(old), iap(old), state(new), iap(new))

    ! Do second predict step.
    call space_operators(state(new), iap(new), tend(old))
    call update_state(dt, tend(old), state(old), iap(old), state(new), iap(new))

    ! Do correct step.
    call space_operators(state(new), iap(new), tend(new))
    ip1 = inner_product(tend(old), tend(new))
    ip2 = inner_product(tend(new), tend(new))
    call log_add_diag('beta', ip1 / ip2)
    if (qcon_modified) then
      dt = time_step_size * ip1 / ip2
    else
      dt = time_step_size
    end if
    call update_state(dt, tend(new), state(old), iap(old), state(new), iap(new))

  end subroutine predict_correct

  subroutine runge_kutta()

  end subroutine runge_kutta

  subroutine leap_frog()

  end subroutine leap_frog

end module dycore_mod
