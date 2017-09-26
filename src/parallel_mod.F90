module parallel_mod

  ! use mpi
  use mesh_mod

  implicit none

  private

  public parallel
  public parallel_init
  public parallel_final

  type parallel_params_type
    integer full_lon_start_idx
    integer full_lon_end_idx
    integer half_lon_start_idx
    integer half_lon_end_idx
    integer full_lat_start_idx
    integer full_lat_end_idx
    integer half_lat_start_idx
    integer half_lat_end_idx
    integer lon_halo_width
    integer lat_halo_width
  end type parallel_params_type

  type(parallel_params_type) parallel

contains

  subroutine parallel_init()

    parallel%full_lon_start_idx = 1
    parallel%full_lat_end_idx = mesh%num_full_lon
    parallel%half_lon_start_idx = 1
    parallel%half_lat_end_idx = mesh%num_half_lon
    parallel%full_lat_start_idx = 1
    parallel%full_lat_end_idx = mesh%num_full_lat
    parallel%half_lat_start_idx = 1
    parallel%half_lat_end_idx = mesh%num_half_lat

    parallel%lon_halo_width = 1
    parallel%lat_halo_width = 1

  end subroutine parallel_init

  subroutine parallel_final()

  end subroutine parallel_final

end module parallel_mod