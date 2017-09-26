module dycore_mod

  use params_mod
  use mesh_mod
  use parallel_mod

  implicit none

  private

  public dycore_init
  public dycore_final

  type coef_type
    ! Coriolis coefficient at full/half meridional grids
    real, allocatable :: full_cori(:)
    ! Curvature coefficient at full/half meridional grids 
    real, allocatable :: half_cori(:)
    real, allocatable :: full_curv(:)
    real, allocatable :: half_curv(:)
    real, allocatable :: full_dlon(:)
    real, allocatable :: half_dlon(:)
    real, allocatable :: full_dlat(:)
    real, allocatable :: half_dlat(:)
  end type coef_type

  type state_type
    real, allocatable :: u(:,:,:)
    real, allocatable :: v(:,:,:)
    real, allocatable :: gd(:,:,:) ! Geopotential depth
    real, allocatable :: ghs(:,:)  ! Surface geopotential
  end type state_type

  ! IAP transformed variables
  real, allocatable :: ut(:,:,:)
  real, allocatable :: vt(:,:,:)
  real, allocatable :: gdt(:,:,:)

  type(coef_type) coef
  type(state_type) state

contains

  subroutine dycore_init()

    integer j

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

    write(6, *) '[Notice]: Dycore module is initialized.'

  end subroutine dycore_init

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

end module dycore_mod
