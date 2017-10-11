module types_mod

  implicit none

  
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

end module types_mod