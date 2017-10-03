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
    real, allocatable :: u(:,:,:)
    real, allocatable :: v(:,:,:)
    real, allocatable :: gd(:,:,:) ! Geopotential depth
    real, allocatable :: ghs(:,:)  ! Surface geopotential
  end type state_type

  type dtend_type
    real, allocatable :: u_adv_lon(:,:,:)
    real, allocatable :: u_adv_lat(:,:,:)
    real, allocatable :: v_adv_lon(:,:,:)
    real, allocatable :: v_adv_lat(:,:,:)
    real, allocatable :: fu(:,:,:)
    real, allocatable :: fv(:,:,:)
    real, allocatable :: u_pgf(:,:,:)
    real, allocatable :: v_pgf(:,:,:)
    real, allocatable :: mass_div_lon(:,:,:)
    real, allocatable :: mass_div_lat(:,:,:)
  end type dtend_type

  ! IAP transformed variables
  real, allocatable :: ut(:,:,:)
  real, allocatable :: vt(:,:,:)
  real, allocatable :: gdt(:,:,:)

  type(coef_type) coef
  type(state_type) state
  type(dtend_type) dtend

  ! 1: predict-correct
  ! 2: runge-kutta
  ! 3: leap-frog
  integer time_scheme_option

contains

  subroutine dycore_init()

    integer i, j

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

    call parallel_allocate(state%u, half_lon=.true.)
    call parallel_allocate(state%v, half_lat=.true.)
    call parallel_allocate(state%gd)
    call parallel_allocate(state%ghs)
    call parallel_allocate(ut, half_lon=.true.)
    call parallel_allocate(vt, half_lat=.true.)
    call parallel_allocate(gdt)

    call parallel_allocate(dtend%u_adv_lon, half_lon=.true.)
    call parallel_allocate(dtend%u_adv_lat, half_lon=.true.)
    call parallel_allocate(dtend%fv, half_lon=.true.)
    call parallel_allocate(dtend%u_pgf, half_lon=.true.)
    call parallel_allocate(dtend%v_adv_lon, half_lat=.true.)
    call parallel_allocate(dtend%v_adv_lat, half_lat=.true.)
    call parallel_allocate(dtend%fu, half_lat=.true.)
    call parallel_allocate(dtend%v_pgf, half_lat=.true.)
    call parallel_allocate(dtend%mass_div_lon)
    call parallel_allocate(dtend%mass_div_lat)

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%ghs(i,j) = i + j
      end do
    end do
    call parallel_fill_halo(state%ghs, left_halo=.true., right_halo=.true., top_halo=.true., bottom_halo=.true.)
    do j = lbound(state%ghs, 2), ubound(state%ghs, 2)
      do i = lbound(state%ghs, 1), ubound(state%ghs, 1)
        if (i == ubound(state%ghs, 1)) then
          write(6, '(F5.1)') state%ghs(i,j)
        else
          write(6, '(F5.1)', advance='no') state%ghs(i,j)
        end if
      end do
    end do

    select case (time_scheme)
    case ('predict-correct')
      time_scheme_option = 1
    case ('runge-kutta')
      time_scheme_option = 2
    case ('leap-frog')
      time_scheme_option = 3
    case default
      write(6, *) '[Error]: Unknown time_scheme ' // trim(time_scheme) // '!'
      stop 1
    end select

    write(6, *) '[Notice]: Dycore module is initialized.'

  end subroutine dycore_init

  subroutine dycore_run()

  end subroutine dycore_run

  subroutine dycore_final()

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
    if (allocated(state%u)) deallocate(state%u)
    if (allocated(state%v)) deallocate(state%v)
    if (allocated(state%gd)) deallocate(state%gd)
    if (allocated(state%ghs)) deallocate(state%ghs)
    if (allocated(ut)) deallocate(ut)
    if (allocated(vt)) deallocate(vt)
    if (allocated(gdt)) deallocate(gdt)
    if (allocated(dtend%u_adv_lon)) deallocate(dtend%u_adv_lon)
    if (allocated(dtend%u_adv_lat)) deallocate(dtend%u_adv_lat)
    if (allocated(dtend%v_adv_lon)) deallocate(dtend%v_adv_lon)
    if (allocated(dtend%v_adv_lat)) deallocate(dtend%v_adv_lat)
    if (allocated(dtend%fu)) deallocate(dtend%fu)
    if (allocated(dtend%fv)) deallocate(dtend%fv)
    if (allocated(dtend%u_pgf)) deallocate(dtend%u_pgf)
    if (allocated(dtend%v_pgf)) deallocate(dtend%v_pgf)
    if (allocated(dtend%mass_div_lon)) deallocate(dtend%mass_div_lon)
    if (allocated(dtend%mass_div_lat)) deallocate(dtend%mass_div_lat)

    write(6, *) '[Notice]: Dycore module is finalized.'

  end subroutine dycore_final

  subroutine iap_transform(time_idx)

    integer, intent(in) :: time_idx

    integer i, j

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        gdt(i,j,time_idx) = sqrt(state%gd(i,j,time_idx))
      end do
    end do

    ! TODO: Should we only fill right and top halo?
    call parallel_fill_halo(gdt(:,:,time_idx))

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        ut(i,j,time_idx) = 0.5 * (gdt(i,j,time_idx) + gdt(i+1,j,time_idx)) * state%u(i,j,time_idx)
      end do
    end do

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        vt(i,j,time_idx) = 0.5 * (gdt(i,j,time_idx) + gdt(i,j+1,time_idx)) * state%v(i,j,time_idx)
      end do
    end do

  end subroutine iap_transform

  subroutine calc_momentum_advection_terms(time_idx)

    integer, intent(in) :: time_idx

    real up1, um1, vp1, vm1
    integer i, j

    ! U
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        up1 = state%u(i,j,time_idx) + state%u(i+1,j,time_idx)
        um1 = state%u(i,j,time_idx) + state%u(i-1,j,time_idx)
        dtend%u_adv_lon(i,j,time_idx) = 0.5 * coef%full_dlon(j) * (up1 * ut(i+1,j,time_idx) - um1 * ut(i-1,j,time_idx))

        vp1 = (state%v(i,j,time_idx) + state%v(i+1,j,time_idx)) * mesh%half_cos_lat(j)
        vm1 = (state%v(i,j-1,time_idx) + state%v(i+1,j-1,time_idx)) * mesh%half_cos_lat(j-1)
        dtend%u_adv_lat(i,j,time_idx) = 0.5 * coef%full_dlat(j) * (vp1 * ut(i,j+1,time_idx) - vm1 * ut(i,j-1,time_idx))
      end do
    end do

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        up1 = state%u(i,j,time_idx) + state%u(i,j+1,time_idx)
        um1 = state%u(i-1,j,time_idx) + state%u(i-1,j+1,time_idx)
        dtend%v_adv_lon(i,j,time_idx) = 0.5 * coef%half_dlon(j) * (up1 * vt(i+1,j,time_idx) - um1 * vt(i-1,j,time_idx))
      end do
    end do

    ! V
    do j = parallel%half_lat_start_idx_no_pole, parallel%half_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        vp1 = (state%v(i,j,time_idx) * mesh%half_cos_lat(j) + state%v(i,j+1,time_idx) * mesh%half_cos_lat(j+1))
        vm1 = (state%v(i,j,time_idx) * mesh%half_cos_lat(j) + state%v(i,j-1,time_idx) * mesh%half_cos_lat(j-1))
        dtend%v_adv_lat(i,j,time_idx) = 0.5 * coef%half_dlat(j) * (vp1 * vt(i,j+1,time_idx) - vm1 * vt(i,j-1,time_idx))
      end do
    end do

    ! Handle meridional advection at North Pole.
    if (parallel%half_lat_start_idx == parallel%half_lat_south_pole_idx) then
      j = parallel%half_lat_south_pole_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        vp1 = state%v(i,j,time_idx) * mesh%half_cos_lat(j) + state%v(i,j+1,time_idx) * mesh%half_cos_lat(j+1)
        dtend%v_adv_lat(i,j,time_idx) = 0.5 * coef%half_dlat(j) * vp1 * vt(i,j+1,time_idx)
      end do
    end if

    ! Handle meridional advection at South Pole.
    if (parallel%half_lat_end_idx == parallel%half_lat_north_pole_idx) then
      j = parallel%half_lat_north_pole_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        vm1 = state%v(i,j,time_idx) * mesh%half_cos_lat(j) + state%v(i,j-1,time_idx) * mesh%half_cos_lat(j-1)
        dtend%v_adv_lat(i,j,time_idx) = - 0.5 * coef%half_dlat(j) * vm1 * vt(i,j-1,time_idx)
      end do
    end if

  end subroutine calc_momentum_advection_terms

  subroutine calc_coriolis_and_curvature_terms(time_idx)

    integer, intent(in) :: time_idx

    integer i, j

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        dtend%fv(i,j,time_idx) = 0.25 * (coef%half_cori(j) * (vt(i,j,time_idx) + vt(i+1,j,time_idx)) + &
                                         coef%half_cori(j-1) * (vt(i,j-1,time_idx) + vt(i+1,j-1,time_idx)))
      end do
    end do

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        dtend%fu(i,j,time_idx) = 0.25 * (coef%full_cori(j) * (ut(i,j,time_idx) + ut(i-1,j,time_idx)) + &
                                         coef%full_cori(j+1) * (ut(i,j+1,time_idx) + ut(i-1,j+1,time_idx)))
      end do
    end do

  end subroutine calc_coriolis_and_curvature_terms

  subroutine calc_pressure_gradient_force_terms(time_idx)

    integer, intent(in) :: time_idx

    integer i, j

    ! U
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        dtend%u_pgf(i,j,time_idx) = coef%full_dlon(j) * (gdt(i,j,time_idx) + gdt(i+1,j,time_idx)) * &
          (state%gd(i+1,j,time_idx) - state%gd(i,j,time_idx) + state%ghs(i+1,j) - state%ghs(i,j))
      end do
    end do

    ! V
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%half_lon_end_idx
        dtend%v_pgf(i,j,time_idx) = coef%half_dlat(j) * (gdt(i,j,time_idx) + gdt(i,j+1,time_idx)) * &
          (state%gd(i,j+1,time_idx) - state%gd(i,j,time_idx) + state%ghs(i,j+1) - state%ghs(i,j)) * mesh%half_cos_lat(j)
      end do
    end do

  end subroutine calc_pressure_gradient_force_terms

  subroutine calc_mass_divergence(time_idx)

      integer, intent(in) :: time_idx

      integer i, j
      real sp, np

      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          dtend%mass_div_lon(i,j,time_idx) = coef%full_dlon(j) * ( &
            (gdt(i,j,time_idx) + gdt(i+1,j,time_idx)) * ut(i,j,time_idx) - &
            (gdt(i,j,time_idx) + gdt(i-1,j,time_idx)) * ut(i-1,j,time_idx))
          dtend%mass_div_lat(i,j,time_idx) = coef%full_dlat(j) * ( &
            (gdt(i,j,time_idx) + gdt(i,j+1,time_idx)) * vt(i,j,time_idx) * mesh%half_cos_lat(j) - &
            (gdt(i,j,time_idx) + gdt(i,j-1,time_idx)) * vt(i,j-1,time_idx) * mesh%half_cos_lat(j-1))
        end do
      end do

      if (parallel%full_lat_start_idx == parallel%full_lat_north_pole_idx) then
        j = parallel%full_lat_north_pole_idx
        sp = 0.0
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          sp = sp + (gdt(i,j,time_idx) + gdt(i,j+1,time_idx)) * vt(i,j,time_idx) * mesh%half_cos_lat(j)
        end do
        call parallel_zonal_sum(sp, dtend%mass_div_lat(1,j,time_idx))
      end if

      if (parallel%full_lat_end_idx == parallel%full_lat_south_pole_idx) then
        j = parallel%full_lat_south_pole_idx
        np = 0.0
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          np = np - (gdt(i,j,time_idx) + gdt(i,j-1,time_idx)) * vt(i,j-1,time_idx) * mesh%half_cos_lat(j-1)
        end do
        call parallel_zonal_sum(np, dtend%mass_div_lat(1,j,time_idx))
      end if

  end subroutine calc_mass_divergence

  subroutine time_integrate()

    select case (time_scheme_option)
    case (1)
      call predict_correct()
    case (2)
      call runge_kutta()
    case (3)
      call leap_frog()
    end select

  end subroutine time_integrate

  subroutine predict_correct()

  end subroutine predict_correct

  subroutine runge_kutta()

  end subroutine runge_kutta

  subroutine leap_frog()

  end subroutine leap_frog

end module dycore_mod
