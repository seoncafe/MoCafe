module memory_mod
  use, intrinsic :: iso_c_binding,   only: C_PTR, C_F_POINTER
  use, intrinsic :: iso_fortran_env, only: character_storage_size, real32, real64, real128, int32, int64
  use mpi
  implicit none
  private
  integer, parameter :: wp = real64
  !-- Rewrote the code (2011-11-02). Now, "window" is internally stored.
  !
  !-- stroage_size returns the storage size in bits and is a standard in Fortran 2008.
  !-- character_storage_size is defined in module iso_fortran_env. (character_storage_size = 8 in most cases, but not always.)
  !-- sizeof returns the size in bytes, but is not a standard. (2017-06-10)
  !integer(kind=MPI_ADDRESS_KIND), parameter :: wp_size = int(sizeof(0.0_wp), MPI_ADDRESS_KIND)
  integer(kind=MPI_ADDRESS_KIND), parameter :: wp_size  = int(storage_size(0.0_wp)/character_storage_size, MPI_ADDRESS_KIND)
  integer(kind=MPI_ADDRESS_KIND), parameter :: int_size = int(storage_size(0)/character_storage_size, MPI_ADDRESS_KIND)
  integer, parameter :: MAX_NUM_WINDOWS = 1000
  integer :: windows(MAX_NUM_WINDOWS)
  integer :: num_windows = 0
  logical :: mpi_parameters_defined = .false.
  integer :: nprocs, p_rank, h_rank, hostcomm, same_hrank_comm, num_hosts

  public create_shared_mem, create_mem, reduce_mem, destroy_shared_mem, destroy_shared_mem_all
  interface create_shared_mem
     module procedure create_shared_mem_1D, create_shared_mem_2D, create_shared_mem_3D, &
                      create_shared_mem_4D, create_shared_mem_5D, &
                      create_shared_mem_1D_int
  end interface create_shared_mem
  interface create_mem
     module procedure create_mem_1D, create_mem_2D, create_mem_3D, create_mem_4D, create_mem_5D, &
                      create_mem_1D_int
  end interface create_mem
  interface reduce_mem
     module procedure reduce_mem_1D, reduce_mem_2D, reduce_mem_3D, reduce_mem_4D, reduce_mem_5D, &
                      reduce_mem_1D_int
  end interface reduce_mem

contains
  !---------------------------
  subroutine init_mpi_parameters()
  integer :: ierr
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, p_rank, ierr)
  call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, hostcomm, ierr)
  call MPI_COMM_RANK(hostcomm, h_rank, ierr)
  call MPI_COMM_SPLIT(MPI_COMM_WORLD, h_rank, p_rank, same_hrank_comm, ierr)
  call MPI_COMM_SIZE(same_hrank_comm, num_hosts, ierr)
  mpi_parameters_defined = .true.
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  end subroutine init_mpi_parameters
  !---------------------------
  subroutine create_shared_mem_1D_int(array,arrshape)
  implicit none
  integer, pointer, intent(inout) :: array(:)
  integer,          intent(in)    :: arrshape(1)
  integer :: ierr

  integer(kind=MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr
  integer     :: disp_unit

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (h_rank == 0) then
     wsize = int(arrshape(1), MPI_ADDRESS_KIND)*int_size
  else
     wsize = 0_MPI_ADDRESS_KIND
  endif

  disp_unit   = int_size
  num_windows = num_windows + 1
  if (num_windows > MAX_NUM_WINDOWS) then
     stop 'Please increase MAX_NUM_WINDOWS'
  endif
  call MPI_WIN_ALLOCATE_SHARED(wsize, disp_unit, MPI_INFO_NULL, hostcomm, ptr, windows(num_windows), ierr)
  if (hostcomm /= 0) then
     call MPI_WIN_SHARED_QUERY(windows(num_windows), 0, wsize, disp_unit, ptr, ierr)
  endif
  call C_F_POINTER(ptr, array, arrshape)
  call MPI_WIN_FENCE(0, windows(num_windows), ierr)

  if (h_rank == 0) then
     array = 0
  endif
  end subroutine create_shared_mem_1D_int
  !---------------------------
  subroutine create_shared_mem_1D(array,arrshape)
  implicit none
  real(kind=wp), pointer, intent(inout) :: array(:)
  integer,                intent(in)    :: arrshape(1)
  integer :: ierr
 
  integer(kind=MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr
  integer     :: disp_unit
  
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (h_rank == 0) then
     wsize = int(arrshape(1), MPI_ADDRESS_KIND)*wp_size
  else
     wsize = 0_MPI_ADDRESS_KIND
  endif
 
  disp_unit   = wp_size
  num_windows = num_windows + 1
  if (num_windows > MAX_NUM_WINDOWS) then
     stop 'Please increase MAX_NUM_WINDOWS'
  endif
  call MPI_WIN_ALLOCATE_SHARED(wsize, disp_unit, MPI_INFO_NULL, hostcomm, ptr, windows(num_windows), ierr)
  if (hostcomm /= 0) then
     call MPI_WIN_SHARED_QUERY(windows(num_windows), 0, wsize, disp_unit, ptr, ierr)
  endif
  call C_F_POINTER(ptr, array, arrshape)
  call MPI_WIN_FENCE(0, windows(num_windows), ierr)

  if (h_rank == 0) then
     array = 0.0_wp
  endif
  end subroutine create_shared_mem_1D
  !---------------------------
  subroutine create_shared_mem_2D(array,arrshape)
  implicit none
  real(kind=wp), pointer, intent(inout) :: array(:,:)
  integer,                intent(in)    :: arrshape(2)
  integer :: ierr

  integer(kind=MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr
  integer     :: disp_unit

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (h_rank == 0) then
     wsize = int(arrshape(1)*arrshape(2), MPI_ADDRESS_KIND)*wp_size
  else
     wsize = 0_MPI_ADDRESS_KIND
  endif

  disp_unit   = wp_size
  num_windows = num_windows + 1
  if (num_windows > MAX_NUM_WINDOWS) then
     stop 'Please increase MAX_NUM_WINDOWS'
  endif
  call MPI_WIN_ALLOCATE_SHARED(wsize, disp_unit, MPI_INFO_NULL, hostcomm, ptr, windows(num_windows), ierr)
  if (h_rank /= 0) then
     call MPI_WIN_SHARED_QUERY(windows(num_windows), 0, wsize, disp_unit, ptr, ierr)
  endif
  call C_F_POINTER(ptr, array, arrshape)
  call MPI_WIN_FENCE(0, windows(num_windows), ierr)

  if (h_rank == 0) then
     array = 0.0_wp
  endif
  end subroutine create_shared_mem_2D
  !---------------------------
  subroutine create_shared_mem_3D(array,arrshape)
  implicit none
  real(kind=wp), pointer, intent(inout) :: array(:,:,:)
  integer,                intent(in)    :: arrshape(3)
  integer :: ierr

  integer(kind=MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr
  integer     :: disp_unit

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (h_rank == 0) then
     wsize = int(arrshape(1)*arrshape(2)*arrshape(3), MPI_ADDRESS_KIND)*wp_size
  else
     wsize = 0_MPI_ADDRESS_KIND
  endif

  disp_unit   = wp_size
  num_windows = num_windows + 1
  if (num_windows > MAX_NUM_WINDOWS) then
     stop 'Please increase MAX_NUM_WINDOWS'
  endif
  call MPI_WIN_ALLOCATE_SHARED(wsize, disp_unit, MPI_INFO_NULL, hostcomm, ptr, windows(num_windows), ierr)
  if (h_rank /= 0) then
     call MPI_WIN_SHARED_QUERY(windows(num_windows), 0, wsize, disp_unit, ptr, ierr)
  endif
  call C_F_POINTER(ptr, array, arrshape)
  call MPI_WIN_FENCE(0, windows(num_windows), ierr)

  if (h_rank == 0) then
     array = 0.0_wp
  endif
  end subroutine create_shared_mem_3D
  !---------------------------
  subroutine create_shared_mem_4D(array,arrshape)
  implicit none
  real(kind=wp), pointer, intent(inout) :: array(:,:,:,:)
  integer,                intent(in)    :: arrshape(4)
  integer :: ierr

  integer(kind=MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr
  integer     :: disp_unit

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (h_rank == 0) then
     wsize = int(int(arrshape(1),int64)*arrshape(2)*arrshape(3)*arrshape(4), MPI_ADDRESS_KIND)*wp_size
  else
     wsize = 0_MPI_ADDRESS_KIND
  endif

  disp_unit   = wp_size
  num_windows = num_windows + 1
  if (num_windows > MAX_NUM_WINDOWS) then
     stop 'Please increase MAX_NUM_WINDOWS'
  endif
  call MPI_WIN_ALLOCATE_SHARED(wsize, disp_unit, MPI_INFO_NULL, hostcomm, ptr, windows(num_windows), ierr)
  if (h_rank /= 0) then
     call MPI_WIN_SHARED_QUERY(windows(num_windows), 0, wsize, disp_unit, ptr, ierr)
  endif
  call C_F_POINTER(ptr, array, arrshape)
  call MPI_WIN_FENCE(0, windows(num_windows), ierr)

  if (h_rank == 0) then
     array = 0.0_wp
  endif
  end subroutine create_shared_mem_4D
  !---------------------------
  subroutine create_shared_mem_5D(array,arrshape)
  implicit none
  real(kind=wp), pointer, intent(inout) :: array(:,:,:,:,:)
  integer,                intent(in)    :: arrshape(5)
  integer :: ierr

  integer(kind=MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr
  integer     :: disp_unit

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (h_rank == 0) then
     wsize = int(int(arrshape(1),int64)*arrshape(2)*arrshape(3)*arrshape(4)*arrshape(5), MPI_ADDRESS_KIND)*wp_size
  else
     wsize = 0_MPI_ADDRESS_KIND
  endif

  disp_unit   = wp_size
  num_windows = num_windows + 1
  if (num_windows > MAX_NUM_WINDOWS) then
     stop 'Please increase MAX_NUM_WINDOWS'
  endif
  call MPI_WIN_ALLOCATE_SHARED(wsize, disp_unit, MPI_INFO_NULL, hostcomm, ptr, windows(num_windows), ierr)
  if (h_rank /= 0) then
     call MPI_WIN_SHARED_QUERY(windows(num_windows), 0, wsize, disp_unit, ptr, ierr)
  endif
  call C_F_POINTER(ptr, array, arrshape)
  call MPI_WIN_FENCE(0, windows(num_windows), ierr)

  if (h_rank == 0) then
     array = 0.0_wp
  endif
  end subroutine create_shared_mem_5D
  !---------------------------
  subroutine create_mem_1D_int(array,arrshape,shared_memory)
  implicit none
  integer, pointer,  intent(inout) :: array(:)
  integer,           intent(in)    :: arrshape(1)
  logical, optional, intent(in)    :: shared_memory
  logical :: use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     call create_shared_mem_1D_int(array,arrshape)
  else
     if (.not. associated(array)) allocate(array(arrshape(1)))
     array = 0.0_wp
  endif
  end subroutine create_mem_1D_int
  !---------------------------
  subroutine create_mem_1D(array,arrshape,shared_memory)
  implicit none
  real(kind=wp), pointer, intent(inout) :: array(:)
  integer,                intent(in)    :: arrshape(1)
  logical, optional,      intent(in)    :: shared_memory
  logical :: use_shared_memory
 
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     call create_shared_mem_1D(array,arrshape)
  else
     if (.not. associated(array)) allocate(array(arrshape(1)))
     array = 0.0_wp
  endif
  end subroutine create_mem_1D
  !---------------------------
  subroutine create_mem_2D(array,arrshape,shared_memory)
  implicit none
  real(kind=wp), pointer, intent(inout) :: array(:,:)
  integer,                intent(in)    :: arrshape(2)
  logical, optional,      intent(in)    :: shared_memory
  logical :: use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     call create_shared_mem_2D(array,arrshape)
  else
     if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2)))
     array = 0.0_wp
  endif
  end subroutine create_mem_2D
  !---------------------------
  subroutine create_mem_3D(array,arrshape,shared_memory)
  implicit none
  real(kind=wp), pointer, intent(inout) :: array(:,:,:)
  integer,                intent(in)    :: arrshape(3)
  logical, optional,      intent(in)    :: shared_memory
  logical :: use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     call create_shared_mem_3D(array,arrshape)
  else
     if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2),arrshape(3)))
     array = 0.0_wp
  endif
  end subroutine create_mem_3D
  !---------------------------
  subroutine create_mem_4D(array,arrshape,shared_memory)
  implicit none
  real(kind=wp), pointer, intent(inout) :: array(:,:,:,:)
  integer,                intent(in)    :: arrshape(4)
  logical, optional,      intent(in)    :: shared_memory
  logical :: use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     call create_shared_mem_4D(array,arrshape)
  else
     if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2),arrshape(3),arrshape(4)))
     array = 0.0_wp
  endif
  end subroutine create_mem_4D
  !---------------------------
  subroutine create_mem_5D(array,arrshape,shared_memory)
  implicit none
  real(kind=wp), pointer, intent(inout) :: array(:,:,:,:,:)
  integer,                intent(in)    :: arrshape(5)
  logical, optional,      intent(in)    :: shared_memory
  logical :: use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     call create_shared_mem_5D(array,arrshape)
  else
     if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2),arrshape(3),arrshape(4),arrshape(5)))
     array = 0.0_wp
  endif
  end subroutine create_mem_5D
  !---------------------------
  subroutine destroy_shared_mem(win,ierr)
  implicit none
  integer, intent(in)    :: win
  integer, intent(inout) :: ierr
  call MPI_WIN_FENCE(0,win,ierr)
  call MPI_WIN_FREE(win,ierr)
  end subroutine destroy_shared_mem
  !---------------------------
  subroutine destroy_shared_mem_all()
  implicit none
  integer :: ierr, ii
  do ii=1, num_windows
     call MPI_WIN_FENCE(0,windows(ii),ierr)
     call MPI_WIN_FREE(windows(ii),ierr)
  enddo
  num_windows = 0
  end subroutine destroy_shared_mem_all
  !---------------------------
  subroutine reduce_mem_1D_int(array, shared_memory)
  implicit none
  integer, pointer,  intent(inout) :: array(:)
  logical, optional, intent(in)    :: shared_memory
  integer :: communicator, ierr, nsize
  logical :: do_reduce, use_shared_memory

  do_reduce         = .false.
  use_shared_memory = .false.
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     communicator = same_hrank_comm
     if (h_rank == 0 .and. num_hosts > 1) do_reduce = .true.
  else
     communicator = MPI_COMM_WORLD
     do_reduce = .true.
  endif
  if (do_reduce) then
     nsize = size(array)
     if (p_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,array(:),nsize,MPI_INTEGER,MPI_SUM,0,communicator,ierr)
     else
        call MPI_REDUCE(array(:),           0,nsize,MPI_INTEGER,MPI_SUM,0,communicator,ierr)
     endif
  endif
  end subroutine reduce_mem_1D_int
  !---------------------------
  subroutine reduce_mem_1D(array, shared_memory)
  implicit none
  real(kind=wp), pointer,  intent(inout) :: array(:)
  logical,       optional, intent(in)    :: shared_memory
  integer :: communicator, ierr, nsize
  logical :: do_reduce, use_shared_memory

  do_reduce         = .false.
  use_shared_memory = .false.
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     communicator = same_hrank_comm
     if (h_rank == 0 .and. num_hosts > 1) do_reduce = .true.
  else
     communicator = MPI_COMM_WORLD
     do_reduce = .true.
  endif
  if (do_reduce) then
     nsize = size(array)
     if (p_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,array(:),nsize,MPI_DOUBLE_PRECISION,MPI_SUM,0,communicator,ierr)
     else
        call MPI_REDUCE(array(:),           0,nsize,MPI_DOUBLE_PRECISION,MPI_SUM,0,communicator,ierr)
     endif
  endif
  end subroutine reduce_mem_1D
  !---------------------------
  subroutine reduce_mem_2D(array, shared_memory)
  implicit none
  real(kind=wp), pointer,  intent(inout) :: array(:,:)
  logical,       optional, intent(in)    :: shared_memory
  integer :: communicator, ierr, nsize
  logical :: do_reduce, use_shared_memory

  do_reduce         = .false.
  use_shared_memory = .false.
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     communicator = same_hrank_comm
     if (h_rank == 0 .and. num_hosts > 1) do_reduce = .true.
  else
     communicator = MPI_COMM_WORLD
     do_reduce = .true.
  endif
  if (do_reduce) then
     nsize = size(array)
     if (p_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,array(:,:),nsize,MPI_DOUBLE_PRECISION,MPI_SUM,0,communicator,ierr)
     else
        call MPI_REDUCE(array(:,:),           0,nsize,MPI_DOUBLE_PRECISION,MPI_SUM,0,communicator,ierr)
     endif
  endif
  end subroutine reduce_mem_2D
  !---------------------------
  subroutine reduce_mem_3D(array, shared_memory)
  implicit none
  real(kind=wp), pointer,  intent(inout) :: array(:,:,:)
  logical,       optional, intent(in)    :: shared_memory
  integer :: communicator, ierr, k, nsize
  logical :: do_reduce, use_shared_memory

  do_reduce         = .false.
  use_shared_memory = .false.
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     communicator = same_hrank_comm
     if (h_rank == 0 .and. num_hosts > 1) do_reduce = .true.
  else
     communicator = MPI_COMM_WORLD
     do_reduce = .true.
  endif
  if (do_reduce) then
     nsize = size(array(:,:,1))
     do k=1, size(array(1,1,:))
        if (p_rank == 0) then
           call MPI_REDUCE(MPI_IN_PLACE,array(:,:,k),nsize,MPI_DOUBLE_PRECISION,MPI_SUM,0,communicator,ierr)
        else
           call MPI_REDUCE(array(:,:,k),           0,nsize,MPI_DOUBLE_PRECISION,MPI_SUM,0,communicator,ierr)
        endif
     enddo
  endif
  end subroutine reduce_mem_3D
  !---------------------------
  subroutine reduce_mem_4D(array, shared_memory)
  implicit none
  real(kind=wp), pointer,  intent(inout) :: array(:,:,:,:)
  logical,       optional, intent(in)    :: shared_memory
  integer :: communicator, ierr, j, k, nsize
  logical :: do_reduce, use_shared_memory

  do_reduce         = .false.
  use_shared_memory = .false.
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     communicator = same_hrank_comm
     if (h_rank == 0 .and. num_hosts > 1) do_reduce = .true.
  else
     communicator = MPI_COMM_WORLD
     do_reduce = .true.
  endif
  if (do_reduce) then
     nsize = size(array(:,:,1,1))
     do k=1, size(array(1,1,1,:))
     do j=1, size(array(1,1,:,1))
        if (p_rank == 0) then
           call MPI_REDUCE(MPI_IN_PLACE,array(:,:,j,k),nsize,MPI_DOUBLE_PRECISION,MPI_SUM,0,communicator,ierr)
        else
           call MPI_REDUCE(array(:,:,j,k),           0,nsize,MPI_DOUBLE_PRECISION,MPI_SUM,0,communicator,ierr)
        endif
     enddo
     enddo
  endif
  end subroutine reduce_mem_4D
  !---------------------------
  subroutine reduce_mem_5D(array, shared_memory)
  implicit none
  real(kind=wp), pointer,  intent(inout) :: array(:,:,:,:,:)
  logical,       optional, intent(in)    :: shared_memory
  integer :: communicator, ierr, i, j, k, nsize
  logical :: do_reduce, use_shared_memory

  do_reduce         = .false.
  use_shared_memory = .false.
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     communicator = same_hrank_comm
     if (h_rank == 0 .and. num_hosts > 1) do_reduce = .true.
  else
     communicator = MPI_COMM_WORLD
     do_reduce = .true.
  endif
  if (do_reduce) then
     nsize = size(array(:,:,1,1,1))
     do k=1, size(array(1,1,1,1,:))
     do j=1, size(array(1,1,1,:,1))
     do i=1, size(array(1,1,:,1,1))
        if (p_rank == 0) then
           call MPI_REDUCE(MPI_IN_PLACE,array(:,:,i,j,k),nsize,MPI_DOUBLE_PRECISION,MPI_SUM,0,communicator,ierr)
        else
           call MPI_REDUCE(array(:,:,i,j,k),           0,nsize,MPI_DOUBLE_PRECISION,MPI_SUM,0,communicator,ierr)
        endif
     enddo
     enddo
     enddo
  endif
  end subroutine reduce_mem_5D
  !---------------------------
end module memory_mod
