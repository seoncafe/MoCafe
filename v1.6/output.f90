module output_mod
contains
!---------------------------------------------------------------------------------
  subroutine output_create(output)
  use define
  implicit none
  type(output_type), intent(inout) :: output

! The pointer should be allocated after the statement of
!    private(output)
!--- setup output array
     output%nx = par%nxim
     output%ny = par%nyim
     output%dx = par%dxim
     output%dy = par%dyim
!     output%dx = atan(grid%rmax/par%obs_distance)/(output%nx/2.0_wp) * rad2deg
!     output%dy = atan(grid%rmax/par%obs_distance)/(output%ny/2.0_wp) * rad2deg
     output%steradian_pix = (output%dx*output%dy)*(deg2rad**2)

     if (par%output_mode == 1) then
        nullify(output%tot)
        nullify(output%tot_sig)
        if (.not. associated(output%tot))     allocate(output%tot(output%nx,output%ny))
        if (.not. associated(output%tot_sig)) allocate(output%tot_sig(output%nx,output%ny))
        output%tot(:,:)     = 0.0_wp
        output%tot_sig(:,:) = 0.0_wp
     else if (par%output_mode == 2) then
        nullify(output%scatt)
        nullify(output%direc)
        nullify(output%scatt_sig)
        nullify(output%direc_sig)
        if (.not. associated(output%scatt))     allocate(output%scatt(output%nx,output%ny))
        if (.not. associated(output%direc))     allocate(output%direc(output%nx,output%ny))
        if (.not. associated(output%scatt_sig)) allocate(output%scatt_sig(output%nx,output%ny))
        if (.not. associated(output%direc_sig)) allocate(output%direc_sig(output%nx,output%ny))
        output%scatt(:,:)     = 0.0_wp
        output%direc(:,:)     = 0.0_wp
        output%scatt_sig(:,:) = 0.0_wp
        output%direc_sig(:,:) = 0.0_wp
     endif
  end subroutine output_create
!-----------------
  subroutine output_destroy(output)
  use define
  implicit none
  type(output_type), intent(inout) :: output
     if (par%output_mode == 1) then
        if (associated(output%tot))       deallocate(output%tot)
        if (associated(output%tot_sig))   deallocate(output%tot_sig)
     else if (par%output_mode == 2) then
        if (associated(output%scatt))     deallocate(output%scatt)
        if (associated(output%direc))     deallocate(output%direc)
        if (associated(output%scatt_sig)) deallocate(output%scatt_sig)
        if (associated(output%direc_sig)) deallocate(output%direc_sig)
     endif
  end subroutine output_destroy
!-----------------
  subroutine output_gather_omp(output,output_sum,my_nphotons)
  use define
  use fits_module
  use mathlib
#ifdef _OPENMP
  use omp_lib
#endif
  implicit none

  type(output_type), intent(inout) :: output
  type(output_type), intent(inout) :: output_sum
  integer(kind=i8b), intent(in)    :: my_nphotons
  ! Local variables
  real(kind=wp),       allocatable :: PSF(:,:)
  integer                          :: status

  !if (par%left_right_fold .and. par%position_angle == 0.0_wp) then
  !   if (par%output_mode == 1) then
  !       output%tot   = (output%tot   + output%tot(par%nxim:1:-1,:))/2.0_wp
  !   else if (par%output_mode == 2) then
  !       output%scatt = (output%scatt + output%scatt(par%nxim:1:-1,:))/2.0_wp
  !       output%direc = (output%direc + output%direc(par%nxim:1:-1,:))/2.0_wp
  !   endif
  !endif

  if (par%output_mode == 1) then
      output%tot   = output%tot   * (dble(par%nphotons)/dble(my_nphotons))
  else if (par%output_mode == 2) then
      output%scatt = output%scatt * (dble(par%nphotons)/dble(my_nphotons))
      output%direc = output%direc * (dble(par%nphotons)/dble(my_nphotons))
  endif

  !$OMP CRITICAL
  if (len_trim(par%psf_file) /= 0) then
     call read_2D(trim(par%psf_file),PSF,status)
#ifdef _OPENMP
     if (omp_get_thread_num() == 0) &
#endif
     write(6,*) 'Convolving images with PSF: ',trim(par%psf_file)
     if (status==0) then
        if (par%output_mode == 1) then
           call convolution_2D(output%tot,PSF)
        else if (par%output_mode == 2) then
           call convolution_2D(output%scatt,PSF)
           call convolution_2D(output%direc,PSF)
        endif
        if (allocated(PSF)) deallocate(PSF)
     endif
  endif
  !$OMP END CRITICAL

  !$OMP CRITICAL
  if (par%output_mode == 1) then
      output_sum%tot       = output_sum%tot       + output%tot
      output_sum%tot_sig   = output_sum%tot_sig   + output%tot**2
  else if (par%output_mode == 2) then
      output_sum%scatt     = output_sum%scatt     + output%scatt
      output_sum%direc     = output_sum%direc     + output%direc
      output_sum%scatt_sig = output_sum%scatt_sig + output%scatt**2
      output_sum%direc_sig = output_sum%direc_sig + output%direc**2
  endif
  !$OMP END CRITICAL
  call output_destroy(output)

  !$OMP BARRIER
  !$OMP SINGLE
  if (par%output_mode == 1) then
     output_sum%tot     = output_sum%tot/par%nthreads
     output_sum%tot_sig = output_sum%tot_sig/par%nthreads - output_sum%tot**2
     output_sum%tot_sig = sqrt(output_sum%tot_sig/(par%nthreads - 1.0_wp))
  else if (par%output_mode == 2) then
     output_sum%scatt     = output_sum%scatt/par%nthreads
     output_sum%direc     = output_sum%direc/par%nthreads
     output_sum%scatt_sig = output_sum%scatt_sig/par%nthreads - output_sum%scatt**2
     output_sum%scatt_sig = sqrt(output_sum%scatt_sig/(par%nthreads - 1.0_wp))
     output_sum%direc_sig = output_sum%direc_sig/par%nthreads - output_sum%direc**2
     output_sum%direc_sig = sqrt(output_sum%direc_sig/(par%nthreads - 1.0_wp))
  endif
  !$OMP END SINGLE

  end subroutine output_gather_omp
!-----------------
#ifdef MPI
  subroutine output_gather_mpi(output,output_sum,my_nphotons)
  use define
  use fits_module
  use mathlib
  use mpi
  implicit none

  type(output_type), intent(inout) :: output
  type(output_type), intent(inout) :: output_sum
  integer(kind=i8b), intent(in)    :: my_nphotons
  ! Local variables
  real(kind=wp),       allocatable :: PSF(:,:)
  integer                          :: status
  integer, dimension(MPI_STATUS_SIZE) :: istatus
  integer :: i,ierr,myid,nsize
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

  !if (par%left_right_fold .and. par%position_angle == 0.0_wp) then
  !   if (par%output_mode == 1) then
  !       output%tot   = (output%tot   + output%tot(par%nxim:1:-1,:))/2.0_wp
  !   else if (par%output_mode == 2) then
  !       output%scatt = (output%scatt + output%scatt(par%nxim:1:-1,:))/2.0_wp
  !       output%direc = (output%direc + output%direc(par%nxim:1:-1,:))/2.0_wp
  !   endif
  !endif

  if (par%output_mode == 1) then
      output%tot   = output%tot   * (dble(par%nphotons)/dble(my_nphotons))
  else if (par%output_mode == 2) then
      output%scatt = output%scatt * (dble(par%nphotons)/dble(my_nphotons))
      output%direc = output%direc * (dble(par%nphotons)/dble(my_nphotons))
  endif

  if (len_trim(par%psf_file) /= 0) then
     call read_2D(trim(par%psf_file),PSF,status)
     if (myid == 0) write(6,*) 'Convolving images with PSF: ',trim(par%psf_file)
     if (status==0) then
        if (par%output_mode == 1) then
            call convolution_2D(output%tot,PSF)
        else if (par%output_mode == 2) then
            call convolution_2D(output%scatt,PSF)
            call convolution_2D(output%direc,PSF)
        endif
        deallocate(PSF)
     endif
  endif

  if (par%output_mode == 1) then
     output%tot_sig   = output%tot**2
  else if (par%output_mode == 2) then
     output%scatt_sig = output%scatt**2
     output%direc_sig = output%direc**2
  endif

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  nsize = output%nx * output%ny

  if (par%iseed == 0 .or. par%iseed == -999) then
     !--- The final result will be always the same when the initial seed is the same.
     if (myid == 0) then
        select case(par%output_mode)
        case(1)
           output_sum%tot     = output%tot
           output_sum%tot_sig = output%tot_sig
           do i=1,par%nthreads - 1
              call MPI_RECV(output%tot,     nsize, MPI_DOUBLE_PRECISION, i, 1, MPI_COMM_WORLD, istatus, ierr)
              call MPI_RECV(output%tot_sig, nsize, MPI_DOUBLE_PRECISION, i, 2, MPI_COMM_WORLD, istatus, ierr)
              output_sum%tot     = output_sum%tot     + output%tot
              output_sum%tot_sig = output_sum%tot_sig + output%tot_sig
           enddo
        case(2)
           output_sum%direc     = output%direc
           output_sum%direc_sig = output%direc_sig
           do i=1,par%nthreads - 1
              call MPI_RECV(output%scatt,     nsize, MPI_DOUBLE_PRECISION, i, 1, MPI_COMM_WORLD, istatus, ierr)
              call MPI_RECV(output%scatt_sig, nsize, MPI_DOUBLE_PRECISION, i, 2, MPI_COMM_WORLD, istatus, ierr)
              call MPI_RECV(output%direc,     nsize, MPI_DOUBLE_PRECISION, i, 3, MPI_COMM_WORLD, istatus, ierr)
              call MPI_RECV(output%direc_sig, nsize, MPI_DOUBLE_PRECISION, i, 4, MPI_COMM_WORLD, istatus, ierr)
              output_sum%scatt     = output_sum%scatt     + output%scatt
              output_sum%direc     = output_sum%direc     + output%direc
              output_sum%scatt_sig = output_sum%scatt_sig + output%scatt_sig
              output_sum%direc_sig = output_sum%direc_sig + output%direc_sig
           enddo
        end select
     else
        select case(par%output_mode)
        case(1)
           call MPI_SEND(output%tot,     nsize, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, ierr)
           call MPI_SEND(output%tot_sig, nsize, MPI_DOUBLE_PRECISION, 0, 2, MPI_COMM_WORLD, ierr)
        case(2)
           call MPI_SEND(output%scatt,     nsize, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, ierr)
           call MPI_SEND(output%scatt_sig, nsize, MPI_DOUBLE_PRECISION, 0, 2, MPI_COMM_WORLD, ierr)
           call MPI_SEND(output%direc,     nsize, MPI_DOUBLE_PRECISION, 0, 3, MPI_COMM_WORLD, ierr)
           call MPI_SEND(output%direc_sig, nsize, MPI_DOUBLE_PRECISION, 0, 4, MPI_COMM_WORLD, ierr)
        end select
     endif
  else
     !--- Use this part if you don't care whether the result changes even when the initial seed is the same.
     if (par%output_mode == 1) then
        call MPI_REDUCE(output%tot,       output_sum%tot,       nsize, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(output%tot_sig,   output_sum%tot_sig,   nsize, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     else if (par%output_mode == 2) then
        call MPI_REDUCE(output%scatt,     output_sum%scatt,     nsize, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(output%direc,     output_sum%direc,     nsize, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(output%scatt_sig, output_sum%scatt_sig, nsize, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(output%direc_sig, output_sum%direc_sig, nsize, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     endif
  endif

  call output_destroy(output)

  if (myid == 0) then
     if (par%output_mode == 1) then
        output_sum%tot     = output_sum%tot/par%nthreads
        output_sum%tot_sig = output_sum%tot_sig/par%nthreads - output_sum%tot**2
        output_sum%tot_sig = sqrt(output_sum%tot_sig/(par%nthreads - 1.0_wp))
     else if (par%output_mode == 2) then
        output_sum%scatt     = output_sum%scatt/par%nthreads
        output_sum%direc     = output_sum%direc/par%nthreads
        output_sum%scatt_sig = output_sum%scatt_sig/par%nthreads - output_sum%scatt**2
        output_sum%scatt_sig = sqrt(output_sum%scatt_sig/(par%nthreads - 1.0_wp))
        output_sum%direc_sig = output_sum%direc_sig/par%nthreads - output_sum%direc**2
        output_sum%direc_sig = sqrt(output_sum%direc_sig/(par%nthreads - 1.0_wp))
     endif
  endif

  end subroutine output_gather_mpi
#endif
!-----------------
  subroutine output_gather_serial(output,output_sum,my_nphotons,do_final)
  use define
  use fits_module
  use mathlib
  implicit none
  type(output_type), intent(inout) :: output
  type(output_type), intent(inout) :: output_sum
  integer(kind=i8b), intent(in)    :: my_nphotons
  logical, optional, intent(in)    :: do_final       ! do_final is only for serial code.
  ! Local variables
  real(kind=wp),       allocatable :: PSF(:,:)
  integer                          :: status
  logical                          :: this_is_final

  ! Use do_final only in serial code.
  if (present(do_final)) then
     this_is_final = .false.
     if (do_final) this_is_final = .true.
  else
     this_is_final = .true.
  endif

  if (par%left_right_fold .and. par%position_angle == 0.0_wp) then
     if (par%output_mode == 1) then
         output%tot   = (output%tot   + output%tot(par%nxim:1:-1,:))/2.0_wp
     else if (par%output_mode == 2) then
         output%scatt = (output%scatt + output%scatt(par%nxim:1:-1,:))/2.0_wp
         output%direc = (output%direc + output%direc(par%nxim:1:-1,:))/2.0_wp
     endif
  endif

  if (par%output_mode == 1) then
      output%tot   = output%tot   * (dble(par%nphotons)/dble(my_nphotons))
  else if (par%output_mode == 2) then
      output%scatt = output%scatt * (dble(par%nphotons)/dble(my_nphotons))
      output%direc = output%direc * (dble(par%nphotons)/dble(my_nphotons))
  endif

  if (len_trim(par%psf_file) /= 0) then
     call read_2D(trim(par%psf_file),PSF,status)
     if (this_is_final) write(6,*) 'Convolving images with PSF: ',trim(par%psf_file)
     if (status==0) then
        if (par%output_mode == 1) then
            call convolution_2D(output%tot,PSF)
        else if (par%output_mode == 2) then
            call convolution_2D(output%scatt,PSF)
            call convolution_2D(output%direc,PSF)
        endif
        deallocate(PSF)
     endif
  endif

  if (par%output_mode == 1) then
      output_sum%tot       = output_sum%tot       + output%tot
      output_sum%tot_sig   = output_sum%tot_sig   + output%tot**2
  else if (par%output_mode == 2) then
      output_sum%scatt     = output_sum%scatt     + output%scatt
      output_sum%direc     = output_sum%direc     + output%direc
      output_sum%scatt_sig = output_sum%scatt_sig + output%scatt**2
      output_sum%direc_sig = output_sum%direc_sig + output%direc**2
  endif
  call output_destroy(output)

  if (this_is_final) then
     if (par%output_mode == 1) then
        output_sum%tot     = output_sum%tot/par%nthreads
        output_sum%tot_sig = output_sum%tot_sig/par%nthreads - output_sum%tot**2
        output_sum%tot_sig = sqrt(output_sum%tot_sig/(par%nthreads - 1.0_wp))
     else if (par%output_mode == 2) then
        output_sum%scatt     = output_sum%scatt/par%nthreads
        output_sum%direc     = output_sum%direc/par%nthreads
        output_sum%scatt_sig = output_sum%scatt_sig/par%nthreads - output_sum%scatt**2
        output_sum%scatt_sig = sqrt(output_sum%scatt_sig/(par%nthreads - 1.0_wp))
        output_sum%direc_sig = output_sum%direc_sig/par%nthreads - output_sum%direc**2
        output_sum%direc_sig = sqrt(output_sum%direc_sig/(par%nthreads - 1.0_wp))
     endif
  endif

  end subroutine output_gather_serial
end module output_mod
