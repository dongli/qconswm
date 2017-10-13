module io_mod

  use netcdf
  use log_mod
  use map_mod
  use params_mod, start_time_in => start_time
  use time_mod
  use string_mod

  implicit none

  private

  public io_init
  public io_create_dataset
  public io_add_dim
  public io_add_var
  public io_start_output
  public io_output
  public io_end_output

  type dataset_type
    integer id
    character(30) name
    character(256) desc
    character(256) author
    character(256) file_prefix
    type(var_type), pointer :: time_var => null()
    type(map_type) dims
    type(map_type) vars
    real output_period
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

  interface io_output
    module procedure io_output_real_1d
    module procedure io_output_real_2d
  end interface io_output

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

  subroutine io_create_dataset(name, desc, author, file_prefix)

    character(*), intent(in), optional :: name
    character(*), intent(in) :: desc
    character(*), intent(in) :: author
    character(*), intent(in) :: file_prefix

    character(30) name_
    type(dataset_type) dataset
    integer i

    if (present(name)) then
      name_ = name
    else
      name_ = 'master'
    end if

    if (datasets%mapped(name_)) then
      write(6, *) '[Error]: Already created dataset ' // trim(name_) // '!'
      stop 1
    end if

    dataset%name = name_
    dataset%desc = desc
    dataset%author = author
    dataset%file_prefix = file_prefix
    i = datasets%size() + 1
    name_ = string_split(output_periods(i), 1)
    read(name_, *) dataset%output_period 
    select case (string_split(output_periods(i), 2))
    case ('days')
      dataset%output_period = dataset%output_period * 86400
    case ('hours')
      dataset%output_period = dataset%output_period * 3600
    case ('minutes')
      dataset%output_period = dataset%output_period * 60
    case ('seconds')
      dataset%output_period = dataset%output_period
    case default
      call log_error('Invalid output period ' // trim(output_periods(i)))
    end select

    ! Add alert for output.
    call time_add_alert('output #' // to_string(i), seconds=dataset%output_period)

    call datasets%insert(dataset%name, dataset)

  end subroutine io_create_dataset

  subroutine io_add_dim(name, dataset_name, long_name, units, size)

    character(*), intent(in) :: name
    character(*), intent(in), optional :: dataset_name
    character(*), intent(in), optional :: long_name
    character(*), intent(in), optional :: units
    integer, intent(in), optional :: size

    type(dataset_type), pointer :: dataset
    type(dim_type) dim
    integer i

    dataset => get_dataset(dataset_name)

    if (dataset%dims%mapped(name)) then
      write(6, *) '[Error]: Already added dimension ' // trim(name) // ' in dataset ' // trim(dataset%name) // '!'
      stop 1
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

    dataset => get_dataset(dataset_name)

    if (dataset%vars%mapped(name)) then
      write(6, *) '[Error]: Already added variable ' // trim(name) // ' in dataset ' // trim(dataset%name) // '!'
      stop 1
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
        write(6, *) '[Error]: Unknown data type ' // trim(data_type) // ' for variable ' // trim(name) // '!'
        stop 1
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
        write(6, *) '[Error]: Unknown dimension ' // trim(dim_names(i)) // ' for variable ' // trim(name) // '!'
        stop 1
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

    dataset => get_dataset(dataset_name)

    write(file_path, "(A, '.', A, '.nc')") trim(dataset%file_prefix), trim(curr_time_format)

    ierr = NF90_CREATE(file_path, NF90_CLOBBER, dataset%id)
    if (ierr /= NF90_NOERR) then
      write(6, *) '[Error]: Failed to create NetCDF file to output!'
      stop 1
    end if
    ierr = NF90_PUT_ATT(dataset%id, NF90_GLOBAL, 'dataset', dataset%name)
    ierr = NF90_PUT_ATT(dataset%id, NF90_GLOBAL, 'desc', dataset%desc)
    ierr = NF90_PUT_ATT(dataset%id, NF90_GLOBAL, 'author', dataset%author)

    iter = map_iterator_type(dataset%dims)
    do while (.not. iter%at_end())
      dim => dataset%get_dim(iter%key())
      ierr = NF90_DEF_DIM(dataset%id, dim%name, dim%size, dim%id)
      if (ierr /= NF90_NOERR) then
        write(6, *) '[Error]: Failed to define dimension ' // trim(dim%name) // '!'
        stop 1
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
        write(6, *) '[Error]: Failed to define variable ' // trim(var%name) // '!'
        stop 1
      end if
      ierr = NF90_PUT_ATT(dataset%id, var%id, 'long_name', trim(var%long_name))
      ierr = NF90_PUT_ATT(dataset%id, var%id, 'units', trim(var%units))
      call iter%next()
    end do

    ierr = NF90_ENDDEF(dataset%id)

    ! Write time dimension variable.
    if (associated(dataset%time_var)) then
      ierr = NF90_PUT_VAR(dataset%id, dataset%time_var%id, &
        dataset%output_period * dataset%time_step / time_units_in_seconds)
      if (ierr /= NF90_NOERR) then
        write(6, *) '[Error]: Failed to write variable time!'
        stop 1
      end if
      dataset%time_step = dataset%time_step + 1
    end if

  end subroutine io_start_output

  subroutine io_output_real_1d(name, array, dataset_name)

    character(*), intent(in) :: name
    real, intent(in) :: array(:)
    character(30), intent(in), optional :: dataset_name

    type(dataset_type), pointer :: dataset
    type(var_type), pointer :: var
    integer ierr    

    dataset => get_dataset(dataset_name)
    var => dataset%get_var(name)

    ierr = NF90_PUT_VAR(dataset%id, var%id, array)
    if (ierr /= NF90_NOERR) then
      write(6, *) '[Error]: Failed to write variable ' // trim(name) // ' in dataset ' // trim(dataset%name) // '!'
      stop 1
    end if

  end subroutine io_output_real_1d

  subroutine io_output_real_2d(name, array, dataset_name)

    character(*), intent(in) :: name
    real, intent(in) :: array(:,:)
    character(30), intent(in), optional :: dataset_name

    type(dataset_type), pointer :: dataset
    type(var_type), pointer :: var
    integer start(3), count(3), i, ierr    

    dataset => get_dataset(dataset_name)
    var => dataset%get_var(name)

    do i = 1, 3
      start(i) = 1
      if (var%dims(i)%ptr%size == NF90_UNLIMITED) then
        count(i) = 1
      else
        count(i) = var%dims(i)%ptr%size
      end if
    end do
    ierr = NF90_PUT_VAR(dataset%id, var%id, array, start, count)
    if (ierr /= NF90_NOERR) then
      write(6, *) '[Error]: Failed to write variable ' // trim(name) // ' in dataset ' // trim(dataset%name) // '!'
      write(6, *) NF90_STRERROR(ierr)
      stop 1
    end if

  end subroutine io_output_real_2d

  subroutine io_end_output(dataset_name)

    character(*), intent(in), optional :: dataset_name

    type(dataset_type), pointer :: dataset
    integer ierr

    dataset => get_dataset(dataset_name)

    ierr = NF90_CLOSE(dataset%id)
    if (ierr /= NF90_NOERR) then
      write(6, *) '[Error]: Failed to close dataset ' // trim(dataset%name) // '!'
      stop 1
    end if

  end subroutine io_end_output

  function get_dataset(name) result(res)

    character(*), intent(in), optional :: name
    type(dataset_type), pointer :: res

    character(30) name_
    class(*), pointer :: value

    if (present(name)) then
      name_ = name
    else
      name_ = 'master'
    end if
    value => datasets%value(name_)
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