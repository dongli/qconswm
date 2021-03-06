module dycore_mod

  use ieee_arithmetic
  use params_mod, time_scheme_in => time_scheme, split_scheme_in => split_scheme
  use log_mod
  use types_mod
  use mesh_mod
  use time_mod
  use parallel_mod
  use io_mod
  use history_mod
  use restart_mod

  implicit none

  private

  public dycore_init
  public dycore_restart
  public dycore_run
  public dycore_final
  public state
  public static

  ! 1: predict_correct
  ! 2: runge_kutta
  ! 3: leap_frog
  ! 4: middle_point
  integer time_scheme
  ! 1: csp1: first order conservative split
  ! 2: csp2: second order conservative split
  ! 3: isp: inproved second order conservative split
  integer split_scheme
  integer, parameter :: all_pass = 0
  integer, parameter :: fast_pass = 1
  integer, parameter :: slow_pass = 2

  type(coef_type) coef
  type(state_type) state(-1:2)
  type(static_type) static
  type(tend_type) tend(0:2)

  integer, parameter :: half_time_idx = 0

contains

  subroutine dycore_init()

    integer i, j, time_idx

    if (case_name == '') then
      call log_error('case_name is not set!')
    end if

    call mesh_init()
    call time_init()
    call parallel_init()
    call io_init()
    call history_init()
    call restart_init()

    call allocate_data(coef)
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

    do time_idx = 0, 2
      call allocate_data(state(time_idx))
    end do
    do time_idx = 0, 2
      call allocate_data(tend(time_idx))
    end do
    call allocate_data(static)

    select case (time_scheme_in)
    case ('predict_correct')
      time_scheme = 1
    case ('runge_kutta')
      time_scheme = 2
    case ('leap_frog')
      time_scheme = 3
    case ('middle_point')
      time_scheme = 4
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

    call log_notice('Dycore module is initialized.')

  end subroutine dycore_init

  subroutine dycore_restart()

    call restart_read(state(old_time_idx))

  end subroutine dycore_restart

  subroutine dycore_run()

    call reset_cos_lat_at_poles()

    call iap_transform(state(old_time_idx))

    call output(state(old_time_idx))
    call log_add_diag('total_mass', total_mass(state(old_time_idx)))
    call log_add_diag('total_energy', total_energy(state(old_time_idx)))
    call log_step()

    do while (.not. time_is_finished())
      call time_integrate()
      call time_advance()
      call output(state(old_time_idx))
      call log_add_diag('total_mass', total_mass(state(old_time_idx)))
      call log_add_diag('total_energy', total_energy(state(old_time_idx)))
      call log_step()
    end do

  end subroutine dycore_run

  subroutine dycore_final()

    integer time_idx

    call mesh_final()
    call parallel_final()
    call history_final()

    call deallocate_data(coef)
    do time_idx = lbound(state, 1), ubound(state, 1)
      call deallocate_data(state(time_idx))
    end do
    do time_idx = lbound(tend, 1), ubound(tend, 1)
      call deallocate_data(tend(time_idx))
    end do
    call deallocate_data(static)

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

    type(state_type), intent(in) :: state

    if (time_is_alerted('hist0.output')) call history_write(state)
    if (time_is_alerted('restart.output')) call restart_write(state)

  end subroutine output

  subroutine iap_transform(state)

    type(state_type), intent(inout) :: state

    integer i, j

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%iap%gd(i,j) = sqrt(state%gd(i,j))
      end do
    end do

    call parallel_fill_halo(state%iap%gd(:,:), all_halo=.true.)

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        state%iap%u(i,j) = 0.5 * (state%iap%gd(i,j) + state%iap%gd(i+1,j)) * state%u(i,j)
      end do
    end do

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%iap%v(i,j) = 0.5 * (state%iap%gd(i,j) + state%iap%gd(i,j+1)) * state%v(i,j)
      end do
    end do

    call parallel_fill_halo(state%iap%u(:,:), all_halo=.true.)
    call parallel_fill_halo(state%iap%v(:,:), all_halo=.true.)

  end subroutine iap_transform

  subroutine space_operators(state, tend, dt, pass)

    use omp_lib

    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend
    real, intent(in) :: dt
    integer, intent(in) :: pass

    real cfl, mean_dlon
    integer i, j, k

    ! Calculate maximum CFL along each zonal circle.
    state%coarse_factor(:) = 1
    mean_dlon = sum(coef%full_dlon) / mesh%num_full_lat
!$omp parallel do private(cfl) schedule(static)
    if (use_zonal_coarse) then
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        state%max_cfl(j) = 0.0
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          cfl = dt * state%iap%gd(i,j) / coef%full_dlon(j)
          if (state%max_cfl(j) < cfl) state%max_cfl(j) = cfl
        end do
        ! Calculate coarsed state.
        if (state%max_cfl(j) > 0.2) then
          ! Find coarse_factor based on excess of CFL.
          do k = 1, size(zonal_coarse_factors)
            !if (state%max_cfl(j) < zonal_coarse_factors(k)) then
            if (mean_dlon / coef%full_dlon(j) < zonal_coarse_factors(k)) then
              state%coarse_factor(j) = zonal_coarse_factors(k)
              exit
            end if
            if (zonal_coarse_factors(k) == 0) then
              ! call log_warning('Insufficient zonal_coarse_factors or time step size is too large!')
              state%coarse_factor(j) = zonal_coarse_factors(k - 1)
              exit
            end if
          end do
          call fine_array_to_coarse_array(state%gd(:,j), state%coarse_gd(:,j), state%coarse_factor(j))
          call fine_array_to_coarse_array(state%iap%u(:,j), state%iap%coarse_u(:,j), state%coarse_factor(j))
          call fine_array_to_coarse_array(static%ghs(:,j), static%coarse_ghs(:,j), state%coarse_factor(j))
          state%iap%coarse_gd(:,j) = sqrt(state%coarse_gd(:,j))
          ! print *, j, state%max_cfl(j)
        end if
      end do
    end if
!$omp end parallel do
    call log_add_diag('num_coarse_zonal', count(state%coarse_factor /= 1))

    select case (pass)
    case (all_pass)
      call zonal_momentum_advection_operator(state, tend)
      call meridional_momentum_advection_operator(state, tend)
      call coriolis_operator(state, tend)
      call curvature_operator(state, tend)
      call zonal_pressure_gradient_force_operator(state, tend)
      call meridional_pressure_gradient_force_operator(state, tend)
      call zonal_mass_divergence_operator(state, tend)
      call meridional_mass_divergence_operator(state, tend)

!$omp parallel do collapse(2) schedule(static)
      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%du(i,j) = - tend%u_adv_lon(i,j) - tend%u_adv_lat(i,j) + tend%fv(i,j) + tend%cv(i,j) - tend%u_pgf(i,j)
        end do
      end do
!$omp end parallel do

!$omp parallel do collapse(2) schedule(static)
      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dv(i,j) = - tend%v_adv_lon(i,j) - tend%v_adv_lat(i,j) - tend%fu(i,j) - tend%cu(i,j) - tend%v_pgf(i,j)
        end do
      end do
!$omp end parallel do

!$omp parallel do collapse(2) schedule(static)
      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dgd(i,j) = - tend%mass_div_lon(i,j) - tend%mass_div_lat(i,j)
        end do
      end do
!$omp end parallel do
    case (slow_pass)
#ifndef NDEBUG
      tend%fv = 0.0
      tend%cv = 0.0
      tend%u_pgf = 0.0
      tend%fu = 0.0
      tend%cu = 0.0
      tend%v_pgf = 0.0
      tend%mass_div_lon = 0.0
      tend%mass_div_lat = 0.0
#endif
      call zonal_momentum_advection_operator(state, tend)
      call meridional_momentum_advection_operator(state, tend)

!$omp parallel do collapse(2) schedule(static)
      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%du(i,j) = - tend%u_adv_lon(i,j) - tend%u_adv_lat(i,j)
        end do
      end do
!$omp end parallel do

!$omp parallel do collapse(2) schedule(static)
      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dv(i,j) = - tend%v_adv_lon(i,j) - tend%v_adv_lat(i,j)
        end do
      end do
!$omp end parallel do

      tend%dgd = 0.0
    case (fast_pass)
#ifndef NDEBUG
      tend%u_adv_lon = 0.0
      tend%u_adv_lat = 0.0
      tend%v_adv_lon = 0.0
      tend%v_adv_lat = 0.0
#endif
      call coriolis_operator(state, tend)
      call curvature_operator(state, tend)
      call zonal_pressure_gradient_force_operator(state, tend)
      call meridional_pressure_gradient_force_operator(state, tend)
      call zonal_mass_divergence_operator(state, tend)
      call meridional_mass_divergence_operator(state, tend)

!$omp parallel do collapse(2) schedule(static)
      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%du(i,j) =    tend%fv(i,j) + tend%cv(i,j) - tend%u_pgf(i,j)
        end do
      end do
!$omp end parallel do

!$omp parallel do collapse(2) schedule(static)
      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dv(i,j) =  - tend%fu(i,j) - tend%cu(i,j) - tend%v_pgf(i,j)
        end do
      end do
!$omp end parallel do

!$omp parallel do collapse(2) schedule(static)
      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dgd(i,j) = - tend%mass_div_lon(i,j) - tend%mass_div_lat(i,j)
        end do
      end do
!$omp end parallel do
    end select

    ! call check_antisymmetry(tend, state)

  end subroutine space_operators

  subroutine zonal_u_momentum_advection(j, coarse_factor, lb, ub, u, iap_u, du)

    integer, intent(in) :: j
    integer, intent(in) :: coarse_factor
    integer, intent(in) :: lb
    integer, intent(in) :: ub
    real, intent(in) :: u(lb:ub)
    real, intent(in) :: iap_u(lb:ub)
    real, intent(out) :: du(lb:ub)

    integer i

!$omp parallel do schedule(static)
    do i = lb + parallel%lon_halo_width, ub - parallel%lon_halo_width
      du(i) = 0.5 / coef%full_dlon(j) / coarse_factor * ((u(i) + u(i+1)) * iap_u(i+1) - (u(i) + u(i-1)) * iap_u(i-1))
    end do
!$omp end parallel do

  end subroutine zonal_u_momentum_advection

  subroutine zonal_v_momentum_advection(j, coarse_factor, lb, ub, u1, u2, iap_v, dv)

    integer, intent(in) :: j
    integer, intent(in) :: coarse_factor
    integer, intent(in) :: lb
    integer, intent(in) :: ub
    real, intent(in) :: u1(lb:ub)
    real, intent(in) :: u2(lb:ub)
    real, intent(in) :: iap_v(lb:ub)
    real, intent(out) :: dv(lb:ub)

    integer i

!$omp parallel do schedule(static)
    do i = lb + parallel%lon_halo_width, ub - parallel%lon_halo_width
      dv(i) = 0.5 / coef%half_dlon(j) / coarse_factor * ((u1(i) + u2(i)) * iap_v(i+1) - (u1(i-1) + u2(i-1)) * iap_v(i-1))
    end do
!$omp end parallel do

  end subroutine zonal_v_momentum_advection

  subroutine zonal_momentum_advection_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    real coarse_tend(parallel%half_lon_lb:parallel%half_lon_ub)
    integer j, factor

    ! U
!$omp parallel do schedule(static)
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      call zonal_u_momentum_advection( &
        j, 1, parallel%half_lon_lb, parallel%half_lon_ub, &
        state%u(:,j), state%iap%u(:,j), tend%u_adv_lon(:,j))
    end do
!$omp end parallel do

    ! V
!$omp parallel do schedule(static)
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      call zonal_v_momentum_advection( &
        j, 1, parallel%full_lon_lb, parallel%full_lon_ub, &
        state%u(:,j), state%u(:,j+1), state%iap%v(:,j), tend%v_adv_lon(:,j))
    end do
!$omp end parallel do

  end subroutine zonal_momentum_advection_operator

  subroutine meridional_momentum_advection_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    real vp1, vm1
    integer i, j

    ! U
!$omp parallel do collapse(2) private(vp1, vm1) schedule(static)
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        vp1 = (state%v(i,j  ) + state%v(i+1,j  )) * mesh%half_cos_lat(j  )
        vm1 = (state%v(i,j-1) + state%v(i+1,j-1)) * mesh%half_cos_lat(j-1)
        tend%u_adv_lat(i,j) = 0.5 / coef%full_dlat(j) * (vp1 * state%iap%u(i,j+1) - vm1 * state%iap%u(i,j-1))
      end do
    end do
!$omp end parallel do

    ! V
!$omp parallel do collapse(2) private(vp1, vm1) schedule(static)
    do j = parallel%half_lat_start_idx_no_pole, parallel%half_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        vp1 = state%v(i,j) * mesh%half_cos_lat(j) + state%v(i,j+1) * mesh%half_cos_lat(j+1)
        vm1 = state%v(i,j) * mesh%half_cos_lat(j) + state%v(i,j-1) * mesh%half_cos_lat(j-1)
        tend%v_adv_lat(i,j) = 0.5 / coef%half_dlat(j) * (vp1 * state%iap%v(i,j+1) - vm1 * state%iap%v(i,j-1))
      end do
    end do
!$omp end parallel do

    ! Handle meridional advection at South Pole.
    if (parallel%has_south_pole) then
      j = parallel%half_lat_south_pole_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        vp1 = state%v(i,j) * mesh%half_cos_lat(j) + state%v(i,j+1) * mesh%half_cos_lat(j+1)
        tend%v_adv_lat(i,j) = 0.5 / coef%half_dlat(j) * vp1 * state%iap%v(i,j+1)
      end do
    end if

    ! Handle meridional advection at North Pole.
    if (parallel%has_north_pole) then
      j = parallel%half_lat_north_pole_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        vm1 = state%v(i,j) * mesh%half_cos_lat(j) + state%v(i,j-1) * mesh%half_cos_lat(j-1)
        tend%v_adv_lat(i,j) = - 0.5 / coef%half_dlat(j) * vm1 * state%iap%v(i,j-1)
      end do
    end if

  end subroutine meridional_momentum_advection_operator

  subroutine coriolis_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    real c1, c2
    integer i, j

!$omp parallel do private(c1, c2) schedule(static)
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      c1 = mesh%half_cos_lat(j  ) / mesh%full_cos_lat(j)
      c2 = mesh%half_cos_lat(j-1) / mesh%full_cos_lat(j)
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%fv(i,j) = 0.25 * coef%cori(j) * (c1 * (state%iap%v(i,j  ) + state%iap%v(i+1,j  )) + &
                                              c2 * (state%iap%v(i,j-1) + state%iap%v(i+1,j-1)))
      end do
    end do
!$omp end parallel do
!$omp parallel do collapse(2) schedule(static)
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%fu(i,j) = 0.25 * (coef%cori(j  ) * (state%iap%u(i,j  ) + state%iap%u(i-1,j  )) + &
                               coef%cori(j+1) * (state%iap%u(i,j+1) + state%iap%u(i-1,j+1)))
      end do
    end do
!$omp end parallel do

#ifndef NDEBUG
    if (parallel%has_south_pole) then
      j = parallel%full_lat_south_pole_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        if (state%iap%u(i,j) /= 0.0) then
          call log_error('U at South Pole should be zero!')
        end if
      end do
    end if
    if (parallel%has_north_pole) then
      j = parallel%full_lat_north_pole_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        if (state%iap%u(i,j) /= 0.0) then
          call log_error('U at North Pole should be zero!')
        end if
      end do
    end if
#endif

  end subroutine coriolis_operator

  subroutine curvature_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    real c1, c2
    integer i, j

!$omp parallel do private(c1, c2) schedule(static)
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      c1 = mesh%half_cos_lat(j  ) / mesh%full_cos_lat(j)
      c2 = mesh%half_cos_lat(j-1) / mesh%full_cos_lat(j)
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%cv(i,j) = 0.25 * coef%curv(j) * state%u(i,j) * &
          (c1 * (state%iap%v(i,j  ) + state%iap%v(i+1,j  )) + &
           c2 * (state%iap%v(i,j-1) + state%iap%v(i+1,j-1)))
      end do
    end do
!$omp end parallel do

!$omp parallel do collapse(2) schedule(static)
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%cu(i,j) = 0.25 * &
          (coef%curv(j  ) * state%u(i,  j  ) * state%iap%u(i,  j  ) + &
           coef%curv(j+1) * state%u(i,  j+1) * state%iap%u(i,  j+1) + &
           coef%curv(j  ) * state%u(i-1,j  ) * state%iap%u(i-1,j  ) + &
           coef%curv(j+1) * state%u(i-1,j+1) * state%iap%u(i-1,j+1))
      end do
    end do
!$omp end parallel do

  end subroutine curvature_operator

  subroutine zonal_pressure_gradient_force(j, coarse_factor, lb, ub, gd, iap_gd, ghs, du)

    integer, intent(in) :: j
    integer, intent(in) :: coarse_factor
    integer, intent(in) :: lb
    integer, intent(in) :: ub
    real, intent(in) :: gd(lb:ub)
    real, intent(in) :: iap_gd(lb:ub)
    real, intent(in) :: ghs(lb:ub)
    real, intent(out) :: du(lb:ub)

    integer i

!$omp parallel do schedule(static)
    do i = lb + parallel%lon_halo_width, ub - parallel%lon_halo_width
      du(i) = (iap_gd(i) + iap_gd(i+1)) / coef%full_dlon(j) / coarse_factor * (gd(i+1) - gd(i) + ghs(i+1) - ghs(i))
    end do
!$omp end parallel do

  end subroutine zonal_pressure_gradient_force

  subroutine zonal_pressure_gradient_force_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    real coarse_tend(parallel%half_lon_lb:parallel%half_lon_ub)
    integer j, factor

!$omp parallel do private(factor, coarse_tend) schedule(static)
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      if (state%coarse_factor(j) == 1) then
        call zonal_pressure_gradient_force( &
          j, 1, parallel%half_lon_lb, parallel%half_lon_ub, &
          state%gd(:,j), state%iap%gd(:,j), static%ghs(:,j), tend%u_pgf(:,j))
      else
        factor = state%coarse_factor(j)
        call zonal_pressure_gradient_force( &
          j, factor, coarse_lb(factor), coarse_ub(factor), &
          state%coarse_gd(:,j), state%iap%coarse_gd(:,j), static%coarse_ghs(:,j), coarse_tend)
        call coarse_tend_to_fine_tend(coarse_tend, tend%u_pgf(:,j), factor)
      end if
    end do
!$omp end parallel do

  end subroutine zonal_pressure_gradient_force_operator

  subroutine meridional_pressure_gradient_force_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer i, j

!$omp parallel do collapse(2) schedule(static)
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%v_pgf(i,j) = (state%iap%gd(i,j) + state%iap%gd(i,j+1)) / coef%half_dlat(j) * &
          (state%gd(i,j+1) - state%gd(i,j) + static%ghs(i,j+1) - static%ghs(i,j)) * mesh%half_cos_lat(j)
      end do
    end do
!$omp end parallel do

  end subroutine meridional_pressure_gradient_force_operator

  subroutine zonal_mass_divergence(j, coarse_factor, lb, ub, iap_gd, iap_u, dgd)

    integer, intent(in) :: j
    integer, intent(in) :: coarse_factor
    integer, intent(in) :: lb
    integer, intent(in) :: ub
    real, intent(in) :: iap_gd(lb:ub)
    real, intent(in) :: iap_u(lb:ub)
    real, intent(out) :: dgd(lb:ub)

    integer i

!$omp parallel do schedule(static)
    do i = lb + parallel%lon_halo_width, ub - parallel%lon_halo_width
      dgd(i) = 1.0 / coef%full_dlon(j) / coarse_factor * ( &
        (iap_gd(i) + iap_gd(i+1)) * iap_u(i) - (iap_gd(i) + iap_gd(i-1)) * iap_u(i-1))
    end do
!$omp end parallel do

  end subroutine zonal_mass_divergence

  subroutine zonal_mass_divergence_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    real coarse_tend(parallel%half_lon_lb:parallel%half_lon_ub)
    integer j, factor

!$omp parallel do private(factor, coarse_tend) schedule(static)
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      if (state%coarse_factor(j) == 1) then
        call zonal_mass_divergence( &
          j, 1, parallel%full_lon_lb, parallel%full_lon_ub, &
          state%iap%gd(:,j), state%iap%u(:,j), tend%mass_div_lon(:,j))
      else
        factor = state%coarse_factor(j)
        call zonal_mass_divergence( &
          j, factor, coarse_lb(factor), coarse_ub(factor), &
          state%iap%coarse_gd(:,j), state%iap%coarse_u(:,j), coarse_tend)
        call coarse_tend_to_fine_tend(coarse_tend, tend%mass_div_lon(:,j), factor)
      end if
    end do
!$omp end parallel do

  end subroutine zonal_mass_divergence_operator

  subroutine meridional_mass_divergence_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer i, j
    real sp, np, sum

!$omp parallel do collapse(2) schedule(static)
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%mass_div_lat(i,j) = 1.0 / coef%full_dlat(j) * ( &
          (state%iap%gd(i,j) + state%iap%gd(i,j+1)) * state%iap%v(i,j  ) * mesh%half_cos_lat(j  ) - &
          (state%iap%gd(i,j) + state%iap%gd(i,j-1)) * state%iap%v(i,j-1) * mesh%half_cos_lat(j-1))
      end do
    end do
!$omp end parallel do

    if (parallel%full_lat_start_idx == parallel%full_lat_south_pole_idx) then
      j = parallel%full_lat_south_pole_idx
      sp = 0.0
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        sp = sp + (state%iap%gd(i,j) + state%iap%gd(i,j+1)) * state%iap%v(i,j) * mesh%half_cos_lat(j)
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
        np = np - (state%iap%gd(i,j) + state%iap%gd(i,j-1)) * state%iap%v(i,j-1) * mesh%half_cos_lat(j-1)
      end do
      call parallel_zonal_sum(np, sum)
      sum = sum / mesh%num_full_lon / coef%full_dlat(j)
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%mass_div_lat(i,j) = sum
      end do
    end if

  end subroutine meridional_mass_divergence_operator

  integer function coarse_lb(coarse_factor) result(lb)

    integer, intent(in) :: coarse_factor

    lb = 1 - parallel%lon_halo_width

  end function coarse_lb

  integer function coarse_ub(coarse_factor) result(ub)

    integer, intent(in) :: coarse_factor

    ub = mesh%num_full_lon / coarse_factor + parallel%lon_halo_width

  end function coarse_ub

  subroutine fine_array_to_coarse_array(fine_array, coarse_array, coarse_factor)

    ! NOTE: Here we allocate more-than-need space for coarse_array.
    real, intent(inout) :: fine_array(:)
    real, intent(out) :: coarse_array(:)
    integer, intent(in) :: coarse_factor

    integer i, j, count, m, n

    n = parallel%lon_halo_width
    coarse_array(:) = 0.0
    j = n + 1
    count = 0
    do i = 1 + n, size(fine_array) - n
      count = count + 1
      ! write(*, '(F20.5)', advance='no') fine_array(i)
      coarse_array(j) = coarse_array(j) + fine_array(i)
      if (count == coarse_factor) then
        coarse_array(j) = coarse_array(j) / coarse_factor
        ! write(*, '(" | ", F20.5)') coarse_array(j)
        j = j + 1
        count = 0
      end if
    end do

    ! Fill halo for coarse_array.
    m = (size(fine_array) - 2 * n) / coarse_factor + 2 * n
    coarse_array(1 : n) = coarse_array(m - 2 * n + 1 : m - n)
    coarse_array(m - n + 1 : m) = coarse_array(1 + n : 2 * n)

  end subroutine fine_array_to_coarse_array

  subroutine coarse_tend_to_fine_tend(coarse_tend, fine_tend, coarse_factor)

    real, intent(in) :: coarse_tend(:)
    real, intent(out) :: fine_tend(:)
    integer, intent(in) :: coarse_factor

    integer i, j, count, n

    n = parallel%lon_halo_width
    j = n + 1
    count = 0
    do i = 1 + n, size(fine_tend) - n
      count = count + 1
      fine_tend(i) = coarse_tend(j)
      if (count == coarse_factor) then
        j = j + 1
        count = 0
      end if
    end do

  end subroutine coarse_tend_to_fine_tend

  subroutine update_state(dt, tend, old_state, new_state)

    real, intent(in) :: dt
    type(tend_type), intent(in) :: tend
    type(state_type), intent(in) :: old_state
    type(state_type), intent(inout) :: new_state

    integer i, j

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        new_state%gd(i,j) = old_state%gd(i,j) + dt * tend%dgd(i,j)
        new_state%iap%gd(i,j) = sqrt(new_state%gd(i,j))
      end do
    end do

    call parallel_fill_halo(new_state%iap%gd(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%gd(:,:), all_halo=.true.)

    ! Update IAP wind state.
!$omp parallel do collapse(2) schedule(static)
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        new_state%iap%u(i,j) = old_state%iap%u(i,j) + dt * tend%du(i,j)
      end do
    end do
!$omp end parallel do
!$omp parallel do collapse(2) schedule(static)
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        new_state%iap%v(i,j) = old_state%iap%v(i,j) + dt * tend%dv(i,j)
      end do
    end do
!$omp end parallel do

    ! Transform from IAP to normal state.
!$omp parallel do collapse(2) schedule(static)
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        new_state%u(i,j) = new_state%iap%u(i,j) * 2.0 / (new_state%iap%gd(i,j) + new_state%iap%gd(i+1,j))
      end do
    end do
!$omp end parallel do
!$omp parallel do collapse(2) schedule(static)
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        new_state%v(i,j) = new_state%iap%v(i,j) * 2.0 / (new_state%iap%gd(i,j) + new_state%iap%gd(i,j+1))
      end do
    end do
!$omp end parallel do

    call parallel_fill_halo(new_state%iap%u(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%iap%v(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%u(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%v(:,:), all_halo=.true.)

  end subroutine update_state

  real function inner_product(tend1, tend2) result(res)

    type(tend_type), intent(in) :: tend1
    type(tend_type), intent(in) :: tend2

    integer i, j

    res = 0.0

!$omp parallel do collapse(2) reduction(+:res) schedule(static)
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        res = res + tend1%du(i,j) * tend2%du(i,j) * mesh%full_cos_lat(j)
      end do
    end do
!$omp end parallel do
!$omp parallel do collapse(2) reduction(+:res) schedule(static)
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        res = res + tend1%dv(i,j) * tend2%dv(i,j) * mesh%half_cos_lat(j)
      end do
    end do
!$omp end parallel do
!$omp parallel do collapse(2) reduction(+:res) schedule(static)
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        res = res + tend1%dgd(i,j) * tend2%dgd(i,j) * mesh%full_cos_lat(j)
      end do
    end do
!$omp end parallel do

  end function inner_product

  real function total_mass(state)

    type(state_type), intent(in) :: state

    integer i, j

    total_mass = 0.0
!$omp parallel do collapse(2) reduction(+:total_mass) schedule(static)
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        total_mass = total_mass + mesh%full_cos_lat(j) * mesh%dlon * mesh%dlat * state%gd(i,j)
      end do
    end do
!$omp end parallel do
    total_mass = total_mass * radius**2

    if (ieee_is_nan(total_mass)) then
      call log_error('Total mass is NaN!')
    end if

  end function total_mass

  real function total_energy(state)

    type(state_type), intent(in) :: state

    integer i, j

    total_energy = 0.0
!$omp parallel do collapse(2) reduction(+:total_energy) schedule(static)
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        total_energy = total_energy + state%iap%u(i,j)**2 * mesh%full_cos_lat(j)
      end do
    end do
!$omp end parallel do
!$omp parallel do collapse(2) reduction(+:total_energy) schedule(static)
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        total_energy = total_energy + state%iap%v(i,j)**2 * mesh%half_cos_lat(j)
      end do
    end do
!$omp end parallel do
!$omp parallel do collapse(2) reduction(+:total_energy) schedule(static)
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        total_energy = total_energy + (state%gd(i,j) + static%ghs(i,j))**2 * mesh%full_cos_lat(j)
      end do
    end do
!$omp end parallel do

    if (ieee_is_nan(total_energy)) then
      call log_error('Total energy is NaN!')
    end if

  end function total_energy

  subroutine time_integrate()

    real subcycle_time_step_size
    integer subcycle, time_idx1, time_idx2

    subcycle_time_step_size = time_step_size / subcycles
    time_idx1 = 0
    time_idx2 = new_time_idx

    select case (time_scheme)
    case (1) ! predict_correct
      select case (split_scheme)
      case (2) ! csp2
        call predict_correct(0.5 * time_step_size, old_time_idx, time_idx1, slow_pass)
        do subcycle = 1, subcycles
          call predict_correct(subcycle_time_step_size, time_idx1, time_idx2, fast_pass)
          call time_swap_indices(time_idx1, time_idx2)
        end do
        call predict_correct(0.5 * time_step_size, time_idx1, new_time_idx, slow_pass)
      case default
        call predict_correct(time_step_size)
      end select
    case (2) ! runge_kutta
      call runge_kutta()
    case (3) ! leap_frog
      call leap_frog()
    case (4) ! middle_point
      call middle_point(time_step_size)
    end select

  end subroutine time_integrate

  subroutine middle_point(time_step_size, old_time_idx_, new_time_idx_, pass_)

    real, intent(in) :: time_step_size
    integer, intent(in), optional :: old_time_idx_
    integer, intent(in), optional :: new_time_idx_
    integer, intent(in), optional :: pass_

    integer old, new, half, pass, iteration
    real dt, e1, e2

    old = merge(old_time_idx_, old_time_idx, present(old_time_idx_))
    new = merge(new_time_idx_, new_time_idx, present(new_time_idx_))
    half = half_time_idx
    pass = merge(pass_, all_pass, present(pass_))
    dt = time_step_size

    call copy_state(state(old), state(new))

    e1 = total_energy(state(old))
    do iteration = 1, 8
      call average_state(state(old), state(new), state(half))
      call space_operators(state(half), tend(old), dt, pass)
      call update_state(dt, tend(old), state(old), state(new))
      e2 = total_energy(state(new))
      if (abs(e2 - e1) * 2 / (e2 + e1) < 5.0e-15) then
        exit
      end if
    end do

  end subroutine middle_point

  subroutine predict_correct(time_step_size, old_time_idx_, new_time_idx_, pass_)

    real, intent(in) :: time_step_size
    integer, intent(in), optional :: old_time_idx_
    integer, intent(in), optional :: new_time_idx_
    integer, intent(in), optional :: pass_

    integer old, new, pass
    real dt, ip1, ip2

    old = merge(old_time_idx_, old_time_idx, present(old_time_idx_))
    new = merge(new_time_idx_, new_time_idx, present(new_time_idx_))
    pass = merge(pass_, all_pass, present(pass_))
    dt = time_step_size * 0.5

    ! Do first predict step.
    call space_operators(state(old), tend(old), dt, pass)
    call update_state(dt, tend(old), state(old), state(new))

    ! Do second predict step.
    call space_operators(state(new), tend(old), dt, pass)
    call update_state(dt, tend(old), state(old), state(new))

    ! Do correct step.
    call space_operators(state(new), tend(new), dt, pass)
    ip1 = inner_product(tend(old), tend(new))
    ip2 = inner_product(tend(new), tend(new))
    call log_add_diag('beta', ip1 / ip2)
    dt = time_step_size * merge(ip1 / ip2, 1.0, qcon_modified)
    call update_state(dt, tend(new), state(old), state(new))

  end subroutine predict_correct

  subroutine runge_kutta()

  end subroutine runge_kutta

  subroutine leap_frog()

  end subroutine leap_frog

  subroutine check_antisymmetry(tend, state)

    type(tend_type), intent(in) :: tend
    type(state_type), intent(in) :: state

    integer i, j
    real ip_u_adv_lon
    real ip_u_adv_lat
    real ip_fv
    real ip_cv
    real ip_u_pgf
    real ip_v_adv_lon
    real ip_v_adv_lat
    real ip_fu
    real ip_cu
    real ip_v_pgf
    real ip_mass_div_lon
    real ip_mass_div_lat

    ip_u_adv_lon = 0.0
    ip_u_adv_lat = 0.0
    ip_fv = 0.0
    ip_cv = 0.0
    ip_u_pgf = 0.0
    ip_v_adv_lon = 0.0
    ip_v_adv_lat = 0.0
    ip_fu = 0.0
    ip_cu = 0.0
    ip_v_pgf = 0.0
    ip_mass_div_lon = 0.0
    ip_mass_div_lat = 0.0

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        ip_u_adv_lon = ip_u_adv_lon + tend%u_adv_lon(i,j) * state%iap%u(i,j) * mesh%full_cos_lat(j)
        ip_u_adv_lat = ip_u_adv_lat + tend%u_adv_lat(i,j) * state%iap%u(i,j) * mesh%full_cos_lat(j)
        ip_fv = ip_fv + tend%fv(i,j) * state%iap%u(i,j) * mesh%full_cos_lat(j)
        ip_cv = ip_cv + tend%cv(i,j) * state%iap%u(i,j) * mesh%full_cos_lat(j)
        ip_u_pgf = ip_u_pgf + tend%u_pgf(i,j) * state%iap%u(i,j) * mesh%full_cos_lat(j)
      end do
    end do

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        ip_v_adv_lon = ip_v_adv_lon + tend%v_adv_lon(i,j) * state%iap%v(i,j) * mesh%half_cos_lat(j)
        ip_v_adv_lat = ip_v_adv_lat + tend%v_adv_lat(i,j) * state%iap%v(i,j) * mesh%half_cos_lat(j)
        ip_fu = ip_fu + tend%fu(i,j) * state%iap%v(i,j) * mesh%half_cos_lat(j)
        ip_cu = ip_cu + tend%cu(i,j) * state%iap%v(i,j) * mesh%half_cos_lat(j)
        ip_v_pgf = ip_v_pgf + tend%v_pgf(i,j) * state%iap%v(i,j) * mesh%half_cos_lat(j)
      end do
    end do

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        ip_mass_div_lon = ip_mass_div_lon + tend%mass_div_lon(i,j) * state%gd(i,j) * mesh%full_cos_lat(j)
        ip_mass_div_lat = ip_mass_div_lat + tend%mass_div_lat(i,j) * state%gd(i,j) * mesh%full_cos_lat(j)
      end do
    end do

    print *, &
      ip_u_adv_lon + &
      ip_v_adv_lon + &
      ip_u_pgf + ip_mass_div_lon + &
      ip_fv - ip_fu + &
      ip_cv - ip_cu + &
      ip_u_adv_lat + &
      ip_v_adv_lat + &
      ip_v_pgf + ip_mass_div_lat

  end subroutine check_antisymmetry

end module dycore_mod
