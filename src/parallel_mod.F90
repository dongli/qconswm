module parallel_mod

#ifndef NO_MPI
  use mpi
#endif
  use mesh_mod

  implicit none

  private

  public parallel
  public parallel_init
  public parallel_final
  public parallel_allocate
  public parallel_fill_halo

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
    integer full_lon_lb, half_lon_lb
    integer full_lon_ub, half_lon_ub
    integer full_lat_lb, half_lat_lb
    integer full_lat_ub, half_lat_ub
  end type parallel_params_type

  type(parallel_params_type) parallel

  interface parallel_allocate
    module procedure parallel_allocate_1
    module procedure parallel_allocate_2
  end interface parallel_allocate

  interface parallel_fill_halo
    module procedure parallel_fill_halo_1
    module procedure parallel_fill_halo_2
  end interface parallel_fill_halo

contains

  subroutine parallel_init()

    parallel%full_lon_start_idx = 1
    parallel%full_lon_end_idx = mesh%num_full_lon
    parallel%half_lon_start_idx = 1
    parallel%half_lon_end_idx = mesh%num_half_lon
    parallel%full_lat_start_idx = 1
    parallel%full_lat_end_idx = mesh%num_full_lat
    parallel%half_lat_start_idx = 1
    parallel%half_lat_end_idx = mesh%num_half_lat

    parallel%lon_halo_width = 1
    parallel%lat_halo_width = 1

    parallel%full_lon_lb = parallel%full_lon_start_idx - parallel%lon_halo_width
    parallel%full_lon_ub = parallel%full_lon_end_idx + parallel%lon_halo_width
    parallel%half_lon_lb = parallel%half_lon_start_idx - parallel%lon_halo_width
    parallel%half_lon_ub = parallel%half_lon_end_idx + parallel%lon_halo_width
    parallel%full_lat_lb = parallel%full_lat_start_idx - parallel%lat_halo_width
    parallel%full_lat_ub = parallel%full_lat_end_idx + parallel%lat_halo_width
    parallel%half_lat_lb = parallel%half_lat_start_idx - parallel%lat_halo_width
    parallel%half_lat_ub = parallel%half_lat_end_idx + parallel%lat_halo_width

    write(6, *) '[Notice]: Parallel module is initialized.'

  end subroutine parallel_init

  subroutine parallel_final()

    write(6, *) '[Notice]: Parallel module is finalized.'

  end subroutine parallel_final

  subroutine parallel_allocate_1(field, half_lon, half_lat)

    real, intent(out), allocatable ::  field(:,:)
    logical, intent(in), optional :: half_lon
    logical, intent(in), optional :: half_lat

    if (present(half_lon) .and. half_lon .and. present(half_lat) .and. half_lat) then
      allocate(field(parallel%half_lon_lb:parallel%half_lon_ub,parallel%half_lat_lb:parallel%half_lat_ub))
    else if (present(half_lon) .and. half_lon) then
      allocate(field(parallel%half_lon_lb:parallel%half_lon_ub,parallel%full_lat_lb:parallel%full_lat_ub))
    else if (present(half_lat) .and. half_lat) then
      allocate(field(parallel%full_lon_lb:parallel%full_lon_ub,parallel%half_lat_lb:parallel%half_lat_ub))
    else
      allocate(field(parallel%full_lon_lb:parallel%full_lon_ub,parallel%full_lat_lb:parallel%full_lat_ub))
    end if

  end subroutine parallel_allocate_1

  subroutine parallel_allocate_2(field, half_lon, half_lat)

    real, intent(out), allocatable ::  field(:,:,:)
    logical, intent(in), optional :: half_lon
    logical, intent(in), optional :: half_lat

    integer time_lb, time_ub

    ! Only support two time levels for now.
    time_lb = 1
    time_ub = 2

    if (present(half_lon) .and. half_lon .and. present(half_lat) .and. half_lat) then
      allocate(field(parallel%half_lon_lb:parallel%half_lon_ub,parallel%half_lat_lb:parallel%half_lat_ub,time_lb:time_ub))
    else if (present(half_lon) .and. half_lon) then
      allocate(field(parallel%half_lon_lb:parallel%half_lon_ub,parallel%full_lat_lb:parallel%full_lat_ub,time_lb:time_ub))
    else if (present(half_lat) .and. half_lat) then
      allocate(field(parallel%full_lon_lb:parallel%full_lon_ub,parallel%half_lat_lb:parallel%half_lat_ub,time_lb:time_ub))
    else
      allocate(field(parallel%full_lon_lb:parallel%full_lon_ub,parallel%full_lat_lb:parallel%full_lat_ub,time_lb:time_ub))
    end if

  end subroutine parallel_allocate_2

  subroutine parallel_fill_halo_1(field, left_halo, right_halo, top_halo, bottom_halo, all_halo)

    real, intent(inout) :: field(:,:)
    logical, intent(in), optional :: left_halo
    logical, intent(in), optional :: right_halo
    logical, intent(in), optional :: top_halo
    logical, intent(in), optional :: bottom_halo
    logical, intent(in), optional :: all_halo

    integer i, j, m, n

    if (present(left_halo) .and. left_halo) then
      m = lbound(field, 1) - 1
      n = ubound(field, 1) - 2 * parallel%lon_halo_width
      do j = lbound(field, 2), ubound(field, 2)
        do i = 1, parallel%lon_halo_width
          field(m+i,j) = field(n+i,j)
        end do
      end do
    end if

    ! |             |                             |              |              |
    ! lb            lb + w                        ub - 2w        ub - w         ub
    if (present(right_halo) .and. right_halo) then
      m = ubound(field, 1) - parallel%lon_halo_width
      n = lbound(field, 1) + parallel%lon_halo_width - 1
      do j = lbound(field, 2), ubound(field, 2)
        do i = 1, parallel%lon_halo_width
          field(m+i,j) = field(n+i,j)
        end do
      end do
    end if

    if (present(top_halo) .and. top_halo) then
      m = lbound(field, 2) - 1
      n = ubound(field, 2) - 2 * parallel%lat_halo_width
      do j = 1, parallel%lat_halo_width
        do i = lbound(field, 1), ubound(field, 1)
          field(i,m+j) = field(i,n+j)
        end do
      end do
    end if

    if (present(bottom_halo) .and. bottom_halo) then
      m = ubound(field, 2) - parallel%lat_halo_width
      n = lbound(field, 2) + parallel%lat_halo_width - 1
      do j = 1, parallel%lat_halo_width
        do i = lbound(field, 1), ubound(field, 1)
          field(i,m+j) = field(i,n+j)
        end do
      end do
    end if

  end subroutine parallel_fill_halo_1

  subroutine parallel_fill_halo_2(field, time_idx, left_halo, right_halo, top_halo, bottom_halo)

    real, intent(inout) :: field(:,:,:)
    integer, intent(in) :: time_idx
    logical, intent(in), optional :: left_halo
    logical, intent(in), optional :: right_halo
    logical, intent(in), optional :: top_halo
    logical, intent(in), optional :: bottom_halo

    call parallel_fill_halo_1(field(:,:,time_idx), left_halo, right_halo, top_halo, bottom_halo)

  end subroutine parallel_fill_halo_2

end module parallel_mod