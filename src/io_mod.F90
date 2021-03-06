module io_mod

  use netcdf
  use log_mod
  use map_mod
  use params_mod, start_time_in => start_time
  use time_mod
  use string_mod
  use parallel_mod

  implicit none

  private

  public io_init
  public io_create_dataset
  public io_add_meta
  public io_add_dim
  public io_add_var
  public io_start_output
  public io_output
  public io_end_output
  public io_start_input
  public io_get_meta
  public io_input
  public io_end_input

  type dataset_type
    integer :: id = -1
    character(30) name
    character(256) desc
    character(256) author
    character(256) file_prefix_or_path
    character(10) mode
    type(var_type), pointer :: time_var => null()
    type(map_type) metas
    type(map_type) dims
    type(map_type) vars
    real period
    integer :: time_step = 0
  contains
    procedure :: get_dim => get_dim_from_dataset
    procedure :: get_var => get_var_from_dataset
  end type dataset_type

  type dim_type
    integer id
    character(30) name
    character(256) long_name
    character(60) units
    integer size
  end type dim_type

  type var_dim_type
    type(dim_type), pointer :: ptr => null()
  end type var_dim_type

  type var_type
    integer id
    character(30) name
    character(256) long_name
    character(60) units
    integer data_type
    type(var_dim_type), allocatable :: dims(:)
  end type var_type

  type(map_type) datasets
  real time_units_in_seconds

  interface io_add_meta
    module procedure io_add_meta_integer
    module procedure io_add_meta_real
    module procedure io_add_meta_string
    module procedure io_add_meta_logical
    module procedure io_add_meta_integer_array
  end interface io_add_meta

  interface io_output
    module procedure io_output_real_1d
    module procedure io_output_real_2d
  end interface io_output

  interface io_get_meta
    module procedure io_get_meta_str
  end interface io_get_meta

  interface io_input
    module procedure io_input_real_2d
  end interface io_input

contains

  subroutine io_init()

    select case (time_units)
    case ('days')
      time_units_in_seconds = 86400.0
    case ('hours')
      time_units_in_seconds = 3600.0
    case ('seconds')
      time_units_in_seconds = 60.0
    case default
      call log_error('Invalid time_units ' // trim(time_units) // '!')
    end select

  end subroutine io_init

  subroutine io_create_dataset(name, desc, file_prefix, file_path, mode, period)

    character(*), intent(in), optional :: name
    character(*), intent(in), optional :: desc
    character(*), intent(in), optional :: file_prefix
    character(*), intent(in), optional :: file_path
    character(*), intent(in), optional :: mode
    character(*), intent(in), optional :: period

    character(30) name_, mode_, period_, period_value, period_unit
    character(256) desc_, file_prefix_, file_path_
    type(dataset_type) dataset
    logical is_exist
    integer i

    if (present(name)) then
      name_ = name
    else
      name_ = 'hist0'
    end if
    if (present(desc)) then
      desc_ = desc
    else
      desc_ = ''
    end if
    if (present(file_prefix)) then
      file_prefix_ = file_prefix
    else
      file_prefix_ = ''
    end if
    if (present(file_path)) then
      file_path_ = file_path
    else
      file_path_ = ''
    end if
    if (present(mode)) then
      mode_ = mode
    else
      mode_ = 'output'
    end if
    if (present(period)) then
      period_ = period
    else
      if (name_(1:4) == 'hist') then
        read(name_(5:len_trim(name_)), '(I' // to_string(len_trim(name_) - 4) // ')') i
        period_ = history_periods(i + 1)
      else
        period_ = 'once'
      end if
    end if
    if (present(file_path) .and. period_ /= 'once') then
      call log_warning('io_create_dataset: Set file_path to "' // trim(file_path_) // &
        '", but period is not once! Reset period.')
      period_ = 'once'
    end if

    if (mode_ == 'input') then
      inquire(file=file_path_, exist=is_exist)
      if (.not. is_exist) then
        call log_error('io_create_dataset: Input file "' // trim(file_path_) // '" does not exist!')
      end if
    end if

    if (datasets%mapped(name_)) then
      call log_error('Already created dataset ' // trim(name_) // '!')
    end if

    dataset%name = name_
    dataset%desc = desc_
    dataset%author = author
    if (file_prefix_ /= '' .and. file_path_ == '') then
      dataset%file_prefix_or_path = file_prefix_
    else if (file_prefix_ == '' .and. file_path_ /= '') then
      dataset%file_prefix_or_path = file_path_
    end if
    if (name_(1:4) == 'hist') then
      dataset%file_prefix_or_path = trim(dataset%file_prefix_or_path) // '.' // trim(string_delete(name_, 'ist'))
    end if
    dataset%mode = mode_
    if (period_ /= 'once') then
      period_value = string_split(period_, 1)
      period_unit = string_split(period_, 2)
      read(period_value, *) dataset%period 
    else
      period_unit = 'once'
    end if
    select case (period_unit)
    case ('days')
      dataset%period = dataset%period * 86400
    case ('hours')
      dataset%period = dataset%period * 3600
    case ('minutes')
      dataset%period = dataset%period * 60
    case ('seconds')
      dataset%period = dataset%period
    case ('steps')
      dataset%period = dataset%period * time_step_size
    case ('once')
      dataset%period = 0
    case default
      call log_error('Invalid IO period ' // trim(period_) // '!')
    end select

    call time_add_alert(trim(dataset%name) // '.' // trim(dataset%mode), seconds=dataset%period)

    call datasets%insert(trim(dataset%name) // '.' // trim(dataset%mode), dataset)

    call log_notice('Create ' // trim(dataset%mode) // ' dataset ' // trim(dataset%file_prefix_or_path) // '.')

  end subroutine io_create_dataset

  subroutine io_add_meta_integer(name, value, dataset_name)

    character(*), intent(in) :: name
    integer, intent(in) :: value
    character(*), intent(in), optional :: dataset_name

    type(dataset_type), pointer :: dataset

    dataset => get_dataset(dataset_name, 'output')

    call dataset%metas%insert(name, value)

  end subroutine io_add_meta_integer

  subroutine io_add_meta_real(name, value, dataset_name)

    character(*), intent(in) :: name
    real, intent(in) :: value
    character(*), intent(in), optional :: dataset_name

    type(dataset_type), pointer :: dataset

    dataset => get_dataset(dataset_name, 'output')

    call dataset%metas%insert(name, value)

  end subroutine io_add_meta_real

  subroutine io_add_meta_string(name, value, dataset_name)

    character(*), intent(in) :: name
    character(*), intent(in) :: value
    character(*), intent(in), optional :: dataset_name

    type(dataset_type), pointer :: dataset

    dataset => get_dataset(dataset_name, 'output')

    call dataset%metas%insert(name, value)

  end subroutine io_add_meta_string

  subroutine io_add_meta_logical(name, value, dataset_name)

    character(*), intent(in) :: name
    logical, intent(in) :: value
    character(*), intent(in), optional :: dataset_name

    type(dataset_type), pointer :: dataset

    dataset => get_dataset(dataset_name, 'output')

    call dataset%metas%insert(name, value)

  end subroutine io_add_meta_logical

  subroutine io_add_meta_integer_array(name, values, dataset_name)

    character(*), intent(in) :: name
    integer, intent(in) :: values(:)
    character(*), intent(in), optional :: dataset_name

    type(dataset_type), pointer :: dataset

    dataset => get_dataset(dataset_name, 'output')

    call dataset%metas%insert(name, to_string(values))

  end subroutine io_add_meta_integer_array

  subroutine io_add_dim(name, dataset_name, long_name, units, size)

    character(*), intent(in) :: name
    character(*), intent(in), optional :: dataset_name
    character(*), intent(in), optional :: long_name
    character(*), intent(in), optional :: units
    integer, intent(in), optional :: size

    type(dataset_type), pointer :: dataset
    type(dim_type) dim
    integer i

    dataset => get_dataset(dataset_name, 'output')

    if (dataset%dims%mapped(name)) then
      call log_error('Already added dimension ' // trim(name) // ' in dataset ' // trim(dataset%name) // '!')
    end if

    dim%name = name
    if (.not. present(long_name)) then
      select case (name)
      case ('lon')
        dim%long_name = 'longitude'
      case ('lat')
        dim%long_name = 'latitude'
      case ('time')
        dim%long_name = 'time'
      case default
        dim%long_name = long_name
      end select
    end if
    if (.not. present(units)) then
      select case (name)
      case ('lon')
        dim%units = 'degrees_east'
      case ('lat')
        dim%units = 'degrees_north'
      case ('time')
        write(dim%units, '(A, " since ", A)') trim(time_units), start_time_format
      case default
        dim%units = units
      end select
    end if
    if (.not. present(size)) then
      dim%size = NF90_UNLIMITED
    else
      dim%size = size
    end if

    call dataset%dims%insert(name, dim)

    ! Add corresponding dimension variable.
    call io_add_var(name, dataset_name, long_name=dim%long_name, units=dim%units, dim_names=[name])

  end subroutine io_add_dim

  subroutine io_add_var(name, dataset_name, long_name, units, data_type, dim_names)

    character(*), intent(in) :: name
    character(*), intent(in), optional :: dataset_name
    character(*), intent(in) :: long_name
    character(*), intent(in) :: units
    character(*), intent(in), optional :: data_type
    character(*), intent(in) :: dim_names(:)

    type(dataset_type), pointer :: dataset
    type(var_type) :: var
    type(map_iterator_type) iter
    type(dim_type), pointer :: dim
    integer i
    logical found
    real real

    dataset => get_dataset(dataset_name, 'output')

    if (dataset%vars%mapped(name)) then
      call log_error('Already added variable ' // trim(name) // ' in dataset ' // trim(dataset%name) // '!')
    end if

    var%name = name
    var%long_name = long_name
    var%units = units

    if (.not. present(data_type)) then
      select case (sizeof(real))
      case (4)
        var%data_type = NF90_FLOAT
      case (8)
        var%data_type = NF90_DOUBLE
      end select
    else
      select case (data_type)
      case ('real(4)')
        var%data_type = NF90_FLOAT
      case ('real(8)')
        var%data_type = NF90_DOUBLE
      case ('integer(4)')
        var%data_type = NF90_INT
      case ('integer(8)')
        var%data_type = NF90_INT64
      case default
        call log_error('Unknown data type ' // trim(data_type) // ' for variable ' // trim(name) // '!')
      end select
    end if

    allocate(var%dims(size(dim_names)))
    do i = 1, size(dim_names)
      found = .false.
      iter = map_iterator_type(dataset%dims)
      do while (.not. iter%at_end())
        dim => dataset%get_dim(iter%key())
        if (dim%name == trim(dim_names(i))) then
          var%dims(i)%ptr => dim
          found = .true.
        end if
        call iter%next()
      end do
      if (.not. found) then
        call log_error('Unknown dimension ' // trim(dim_names(i)) // ' for variable ' // trim(name) // '!')
      end if
    end do

    call dataset%vars%insert(name, var)

    if (name == 'time') dataset%time_var => dataset%get_var(name)

  end subroutine io_add_var

  subroutine io_start_output(dataset_name)

    character(*), intent(in), optional :: dataset_name

    character(256) file_path
    type(dataset_type), pointer :: dataset
    type(map_iterator_type) iter
    type(dim_type), pointer :: dim
    type(var_type), pointer :: var
    integer i, ierr, dimids(10)

    dataset => get_dataset(dataset_name, 'output')

    write(file_path, "(A, '.', A, '.nc')") trim(dataset%file_prefix_or_path), trim(curr_time_format)

    ierr = NF90_CREATE(file_path, NF90_CLOBBER, dataset%id)
    if (ierr /= NF90_NOERR) then
      call log_error('Failed to create NetCDF file to output!')
    end if
    ierr = NF90_PUT_ATT(dataset%id, NF90_GLOBAL, 'dataset', dataset%name)
    ierr = NF90_PUT_ATT(dataset%id, NF90_GLOBAL, 'desc', dataset%desc)
    ierr = NF90_PUT_ATT(dataset%id, NF90_GLOBAL, 'author', dataset%author)

    iter = map_iterator_type(dataset%metas)
    do while (.not. iter%at_end())
      select type (value => iter%value())
      type is (integer)
        ierr = NF90_PUT_ATT(dataset%id, NF90_GLOBAL, iter%key(), value)
      type is (real)
        ierr = NF90_PUT_ATT(dataset%id, NF90_GLOBAL, iter%key(), value)
      type is (character(*))
        ierr = NF90_PUT_ATT(dataset%id, NF90_GLOBAL, iter%key(), value)
      type is (logical)
        ierr = NF90_PUT_ATT(dataset%id, NF90_GLOBAL, iter%key(), to_string(value))
      end select
      call iter%next()
    end do

    iter = map_iterator_type(dataset%dims)
    do while (.not. iter%at_end())
      dim => dataset%get_dim(iter%key())
      ierr = NF90_DEF_DIM(dataset%id, dim%name, dim%size, dim%id)
      if (ierr /= NF90_NOERR) then
        call log_error('Failed to define dimension ' // trim(dim%name) // '!')
      end if
      call iter%next()
    end do

    iter = map_iterator_type(dataset%vars)
    do while (.not. iter%at_end())
      var => dataset%get_var(iter%key())
      do i = 1, size(var%dims)
        dimids(i) = var%dims(i)%ptr%id
      end do
      ierr = NF90_DEF_VAR(dataset%id, var%name, var%data_type, dimids(1:size(var%dims)), var%id)
      if (ierr /= NF90_NOERR) then
        call log_error('Failed to define variable ' // trim(var%name) // '!')
      end if
      ierr = NF90_PUT_ATT(dataset%id, var%id, 'long_name', trim(var%long_name))
      ierr = NF90_PUT_ATT(dataset%id, var%id, 'units', trim(var%units))
      call iter%next()
    end do

    ierr = NF90_ENDDEF(dataset%id)

    ! Write time dimension variable.
    if (associated(dataset%time_var)) then
      ierr = NF90_PUT_VAR(dataset%id, dataset%time_var%id, time_elapsed_seconds() / time_units_in_seconds)
      if (ierr /= NF90_NOERR) then
        call log_error('Failed to write variable time!')
      end if
      dataset%time_step = dataset%time_step + 1
    end if

  end subroutine io_start_output

  subroutine io_output_real_1d(name, array, dataset_name)

    character(*), intent(in) :: name
    real, intent(in) :: array(:)
    character(*), intent(in), optional :: dataset_name

    type(dataset_type), pointer :: dataset
    type(var_type), pointer :: var
    integer ierr    

    dataset => get_dataset(dataset_name, 'output')
    var => dataset%get_var(name)

    ierr = NF90_PUT_VAR(dataset%id, var%id, array)
    if (ierr /= NF90_NOERR) then
      call log_error('Failed to write variable ' // trim(name) // ' in dataset ' // trim(dataset%name) // '!')
    end if

  end subroutine io_output_real_1d

  subroutine io_output_real_2d(name, array, dataset_name)

    character(*), intent(in) :: name
    real, intent(in) :: array(:,:)
    character(*), intent(in), optional :: dataset_name

    type(dataset_type), pointer :: dataset
    type(var_type), pointer :: var
    integer lb1, ub1, lb2, ub2
    integer i, j, ierr, varid
    integer start(3), count(3)    
    real, allocatable :: buffer(:,:)

    dataset => get_dataset(dataset_name, 'output')
    var => dataset%get_var(name)

    lb1 = lbound(array, 1) + parallel%lon_halo_width
    ub1 = ubound(array, 1) - parallel%lon_halo_width
    lb2 = lbound(array, 2) + parallel%lat_halo_width
    ub2 = ubound(array, 2) - parallel%lat_halo_width
    allocate(buffer(lb1:ub1,lb2:ub2))

    do j = lb2, ub2
      do i = lb1, ub1
        buffer(i,j) = array(i,j)
      end do
    end do

    do i = 1, 3
      start(i) = 1
      if (var%dims(i)%ptr%size == NF90_UNLIMITED) then
        count(i) = 1
      else
        count(i) = var%dims(i)%ptr%size
      end if
    end do
    ierr = NF90_PUT_VAR(dataset%id, var%id, buffer, start, count)
    if (ierr /= NF90_NOERR) then
      call log_error('Failed to write variable ' // trim(name) // ' to ' // trim(dataset%name) // '!' // NF90_STRERROR(ierr))
    end if

    deallocate(buffer)

  end subroutine io_output_real_2d

  subroutine io_end_output(dataset_name)

    character(*), intent(in), optional :: dataset_name

    type(dataset_type), pointer :: dataset
    integer ierr

    dataset => get_dataset(dataset_name, 'output')

    ierr = NF90_CLOSE(dataset%id)
    if (ierr /= NF90_NOERR) then
      call log_error('Failed to close dataset ' // trim(dataset%name) // '!')
    end if

  end subroutine io_end_output

  subroutine io_start_input(dataset_name)

    character(*), intent(in), optional :: dataset_name

    type(dataset_type), pointer :: dataset
    integer ierr

    dataset => get_dataset(dataset_name, 'input')

    ierr = NF90_OPEN(dataset%file_prefix_or_path, NF90_NOWRITE, dataset%id)
    if (ierr /= NF90_NOERR) then
      call log_error('Failed to open NetCDF file to input! ' // trim(NF90_STRERROR(ierr)))
    end if

  end subroutine io_start_input

  function io_get_meta_str(name, dataset_name) result(res)

    character(*), intent(in) :: name
    character(*), intent(in), optional :: dataset_name
    character(:), allocatable :: res

    type(dataset_type), pointer :: dataset
    character(256) meta
    integer ierr

    dataset => get_dataset(dataset_name, 'input')

    ierr = NF90_GET_ATT(dataset%id, NF90_GLOBAL, name, meta)
    if (ierr /= NF90_NOERR) then
      call log_error('Failed to get meta "' // trim(name) // '" from file ' // trim(dataset%file_prefix_or_path) // '!')
    end if
    res = trim(meta)

  end function io_get_meta_str

  subroutine io_input_real_2d(name, array, dataset_name)

    character(*), intent(in) :: name
    real, intent(out) :: array(:,:)
    character(*), intent(in), optional :: dataset_name

    type(dataset_type), pointer :: dataset
    integer lb1, ub1, lb2, ub2
    integer i, j, ierr, varid
    real, allocatable :: buffer(:,:)

    dataset => get_dataset(dataset_name, 'input')

    lb1 = lbound(array, 1) + parallel%lon_halo_width
    ub1 = ubound(array, 1) - parallel%lon_halo_width
    lb2 = lbound(array, 2) + parallel%lat_halo_width
    ub2 = ubound(array, 2) - parallel%lat_halo_width
    allocate(buffer(lb1:ub1,lb2:ub2))

    ierr = NF90_INQ_VARID(dataset%id, name, varid)
    if (ierr /= NF90_NOERR) then
      call log_error('No variable "' // trim(name) // '" in dataset "' // trim(dataset%file_prefix_or_path) // '"!')
    end if
    ierr = NF90_GET_VAR(dataset%id, varid, buffer)

    do j = lb2, ub2
      do i = lb1, ub1
        array(i,j) = buffer(i,j)
      end do
    end do

    deallocate(buffer)

  end subroutine io_input_real_2d

  subroutine io_end_input(dataset_name)

    character(*), intent(in), optional :: dataset_name

    type(dataset_type), pointer :: dataset
    integer ierr

    dataset => get_dataset(dataset_name, 'input')

    ierr = NF90_CLOSE(dataset%id)

  end subroutine io_end_input

  function get_dataset(name, mode) result(res)

    character(*), intent(in), optional :: name
    character(*), intent(in), optional :: mode
    type(dataset_type), pointer :: res

    character(30) name_, mode_
    class(*), pointer :: value

    if (present(name)) then
      name_ = name
    else
      name_ = 'hist0'
    end if
    if (present(mode)) then
      mode_ = mode
    else
      mode_ = 'output'
    end if
    value => datasets%value(trim(name_) // '.' // trim(mode_))
    select type (value)
    type is (dataset_type)
      res => value
    end select

  end function get_dataset

  function get_dim_from_dataset(this, name) result(res)

    class(dataset_type), intent(in) :: this
    character(*), intent(in) :: name
    type(dim_type), pointer :: res

    class(*), pointer :: value

    value => this%dims%value(name)
    select type (value)
    type is (dim_type)
      res => value
    end select

  end function get_dim_from_dataset

  function get_var_from_dataset(this, name) result(res)

    class(dataset_type), intent(in) :: this
    character(*), intent(in) :: name
    type(var_type), pointer :: res

    class(*), pointer :: value

    value => this%vars%value(name)
    select type (value)
    type is (var_type)
      res => value
    end select

  end function get_var_from_dataset

end module io_mod
