module types_mod

  use mesh_mod
  use parallel_mod

  implicit none

  private

  public allocate_data
  public deallocate_data
  public coef_type
  public state_type
  public static_type
  public iap_type
  public tend_type

  type coef_type
    ! Coriolis coefficient at full meridional grids
    real, allocatable :: cori(:)
    ! Curvature coefficient at full meridional grids 
    real, allocatable :: curv(:)
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
    ! Wind on A grid
    real, allocatable :: ua(:,:)
    real, allocatable :: va(:,:)
  end type state_type

  type static_type
    real, allocatable :: ghs(:,:) ! Surface geopotential
  end type static_type

  type tend_type
    real, allocatable :: u_adv_lon(:,:)
    real, allocatable :: u_adv_lat(:,:)
    real, allocatable :: v_adv_lon(:,:)
    real, allocatable :: v_adv_lat(:,:)
    real, allocatable :: fu(:,:)
    real, allocatable :: fv(:,:)
    real, allocatable :: cu(:,:)
    real, allocatable :: cv(:,:)
    real, allocatable :: u_pgf(:,:)
    real, allocatable :: v_pgf(:,:)
    real, allocatable :: mass_div_lon(:,:)
    real, allocatable :: mass_div_lat(:,:)
    real, allocatable :: du(:,:)
    real, allocatable :: dv(:,:)
    real, allocatable :: dgd(:,:)

    real, allocatable :: a_u_adv_lon(:,:)
    real, allocatable :: a_u_adv_lat(:,:)
    real, allocatable :: a_v_adv_lon(:,:)
    real, allocatable :: a_v_adv_lat(:,:)
    real, allocatable :: a_fu(:,:)
    real, allocatable :: a_fv(:,:)
    real, allocatable :: a_cu(:,:)
    real, allocatable :: a_cv(:,:)
    real, allocatable :: a_u_pgf(:,:)
    real, allocatable :: a_v_pgf(:,:)
    real, allocatable :: a_du(:,:)
    real, allocatable :: a_dv(:,:)
  end type tend_type

  ! IAP transformed variables
  type iap_type
    real, allocatable :: u(:,:)
    real, allocatable :: v(:,:)
    real, allocatable :: gd(:,:)
  end type iap_type

  interface allocate_data
    module procedure allocate_coef_data
    module procedure allocate_static_data
    module procedure allocate_state_data
    module procedure allocate_iap_data
    module procedure allocate_tend_data
  end interface allocate_data

  interface deallocate_data
    module procedure deallocate_coef_data
    module procedure deallocate_static_data
    module procedure deallocate_state_data
    module procedure deallocate_iap_data
    module procedure deallocate_tend_data
  end interface deallocate_data

contains

  subroutine allocate_coef_data(coef)

    type(coef_type), intent(out) :: coef

    allocate(coef%cori(mesh%num_full_lat))
    allocate(coef%curv(mesh%num_full_lat))
    allocate(coef%full_dlon(mesh%num_full_lat))
    allocate(coef%half_dlon(mesh%num_half_lat))
    allocate(coef%full_dlat(mesh%num_full_lat))
    allocate(coef%half_dlat(mesh%num_half_lat))

  end subroutine allocate_coef_data

  subroutine allocate_static_data(static)

    type(static_type), intent(out) :: static

    if (.not. allocated(static%ghs)) call parallel_allocate(static%ghs)

  end subroutine allocate_static_data

  subroutine allocate_state_data(state)

    type(state_type), intent(out) :: state

    if (.not. allocated(state%u)) call parallel_allocate(state%u, half_lon=.true.)
    if (.not. allocated(state%v)) call parallel_allocate(state%v, half_lat=.true.)
    if (.not. allocated(state%gd)) call parallel_allocate(state%gd)
    if (.not. allocated(state%ua)) call parallel_allocate(state%ua)
    if (.not. allocated(state%va)) call parallel_allocate(state%va)

  end subroutine allocate_state_data

  subroutine allocate_iap_data(iap)

    type(iap_type), intent(out) :: iap

    if (.not. allocated(iap%u)) call parallel_allocate(iap%u, half_lon=.true.)
    if (.not. allocated(iap%v)) call parallel_allocate(iap%v, half_lat=.true.)
    if (.not. allocated(iap%gd)) call parallel_allocate(iap%gd)

  end subroutine allocate_iap_data

  subroutine allocate_tend_data(tend)

    type(tend_type), intent(out) :: tend

    if (.not. allocated(tend%u_adv_lon)) call parallel_allocate(tend%u_adv_lon, half_lon=.true.)
    if (.not. allocated(tend%u_adv_lat)) call parallel_allocate(tend%u_adv_lat, half_lon=.true.)
    if (.not. allocated(tend%fv)) call parallel_allocate(tend%fv, half_lon=.true.)
    if (.not. allocated(tend%cv)) call parallel_allocate(tend%cv, half_lon=.true.)
    if (.not. allocated(tend%u_pgf)) call parallel_allocate(tend%u_pgf, half_lon=.true.)
    if (.not. allocated(tend%v_adv_lon)) call parallel_allocate(tend%v_adv_lon, half_lat=.true.)
    if (.not. allocated(tend%v_adv_lat)) call parallel_allocate(tend%v_adv_lat, half_lat=.true.)
    if (.not. allocated(tend%fu)) call parallel_allocate(tend%fu, half_lat=.true.)
    if (.not. allocated(tend%cu)) call parallel_allocate(tend%cu, half_lat=.true.)
    if (.not. allocated(tend%v_pgf)) call parallel_allocate(tend%v_pgf, half_lat=.true.)
    if (.not. allocated(tend%mass_div_lon)) call parallel_allocate(tend%mass_div_lon)
    if (.not. allocated(tend%mass_div_lat)) call parallel_allocate(tend%mass_div_lat)
    if (.not. allocated(tend%du)) call parallel_allocate(tend%du, half_lon=.true.)
    if (.not. allocated(tend%dv)) call parallel_allocate(tend%dv, half_lat=.true.)
    if (.not. allocated(tend%dgd)) call parallel_allocate(tend%dgd)

  end subroutine allocate_tend_data

  subroutine deallocate_coef_data(coef)

    type(coef_type), intent(inout) :: coef

    if (allocated(coef%cori)) deallocate(coef%cori)
    if (allocated(coef%curv)) deallocate(coef%curv)
    if (allocated(coef%full_dlon)) deallocate(coef%full_dlon)
    if (allocated(coef%half_dlon)) deallocate(coef%half_dlon)
    if (allocated(coef%full_dlat)) deallocate(coef%full_dlat)
    if (allocated(coef%half_dlat)) deallocate(coef%half_dlat)

  end subroutine deallocate_coef_data

  subroutine deallocate_static_data(static)

    type(static_type), intent(inout) :: static

    if (allocated(static%ghs)) deallocate(static%ghs)

  end subroutine deallocate_static_data

  subroutine deallocate_state_data(state)

    type(state_type), intent(inout) :: state

    if (allocated(state%u)) deallocate(state%u)
    if (allocated(state%v)) deallocate(state%v)
    if (allocated(state%gd)) deallocate(state%gd)
    if (allocated(state%ua)) deallocate(state%ua)
    if (allocated(state%va)) deallocate(state%va)

  end subroutine deallocate_state_data

  subroutine deallocate_iap_data(iap)

    type(iap_type), intent(inout) :: iap

    if (allocated(iap%u)) deallocate(iap%u)
    if (allocated(iap%v)) deallocate(iap%v)
    if (allocated(iap%gd)) deallocate(iap%gd)

  end subroutine deallocate_iap_data

  subroutine deallocate_tend_data(tend)

    type(tend_type), intent(inout) :: tend

    if (allocated(tend%u_adv_lon)) deallocate(tend%u_adv_lon)
    if (allocated(tend%u_adv_lat)) deallocate(tend%u_adv_lat)
    if (allocated(tend%v_adv_lon)) deallocate(tend%v_adv_lon)
    if (allocated(tend%v_adv_lat)) deallocate(tend%v_adv_lat)
    if (allocated(tend%fu)) deallocate(tend%fu)
    if (allocated(tend%fv)) deallocate(tend%fv)
    if (allocated(tend%cu)) deallocate(tend%cu)
    if (allocated(tend%cv)) deallocate(tend%cv)
    if (allocated(tend%u_pgf)) deallocate(tend%u_pgf)
    if (allocated(tend%v_pgf)) deallocate(tend%v_pgf)
    if (allocated(tend%mass_div_lon)) deallocate(tend%mass_div_lon)
    if (allocated(tend%mass_div_lat)) deallocate(tend%mass_div_lat)
    if (allocated(tend%du)) deallocate(tend%du)
    if (allocated(tend%dv)) deallocate(tend%dv)
    if (allocated(tend%dgd)) deallocate(tend%dgd)

  end subroutine deallocate_tend_data

end module types_mod