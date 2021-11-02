module read_mod
  use define
  use fitsio_mod
contains
  !=================================================================
  ! get 3D dimeisions from density fits file.
  !subroutine get_dimension(fname,nx,ny,nz,xmax,ymax,zmax,status,reduce_factor)
  subroutine get_dimension(fname,nx,ny,nz,status,reduce_factor)
  implicit none
  character(len=*), intent(in)  :: fname
  integer,          intent(out) :: nx,ny,nz
  !real(kind=wp),    intent(out) :: xmax,ymax,zmax
  integer,           intent(inout) :: status
  integer, optional, intent(inout) :: reduce_factor

  !--- local variables
  integer :: unit

  status = 0
  call fits_open_old(unit,trim(fname),status)
  if (status == 0) then
     call fits_get_keyword(unit,'NAXIS1',nx,status)
     call fits_get_keyword(unit,'NAXIS2',ny,status)
     call fits_get_keyword(unit,'NAXIS3',nz,status)
     !call fits_get_keyword(unit,'xmax',xmax,status)
     !call fits_get_keyword(unit,'ymax',ymax,status)
     !call fits_get_keyword(unit,'zmax',zmax,status)
     if (present(reduce_factor)) then
        if (reduce_factor > 1) then
           if ((nx/reduce_factor) > 0 .and. (ny/reduce_factor) > 0 .and. (nz/reduce_factor) > 0) then
              nx = nx/reduce_factor
              ny = ny/reduce_factor
              nz = nz/reduce_factor
           else
              reduce_factor = 1
              write(*,*) 'Oops, Using the original size...'
           endif
        endif
     endif
  else
     write(*,*) 'Error in reading the data file :', trim(fname),' status = ',status
  endif
  call fits_close(unit,status)
  end subroutine
  !=================================================================
  subroutine read_3D(fname,array,reduce_factor,centering)
  implicit none
  character(len=*),    intent(in) :: fname
  real(kind=wp),    intent(inout) :: array(:,:,:)
  integer, optional,   intent(in) :: reduce_factor, centering

  !--- local variables
  integer :: unit, status = 0
  integer :: n1,n2,n3
  integer :: loc(3),n1cen,n2cen,n3cen
  integer :: i,j,k,i1,i2,j1,j2,k1,k2
  !--- note the following should be always a single precision variable.
  real(kind=wp), allocatable :: arr(:,:,:)

  call fits_open_old(unit,trim(fname),status)
  if (status == 0) then
     call fits_get_keyword(unit,'NAXIS1',n1,status)
     call fits_get_keyword(unit,'NAXIS2',n2,status)
     call fits_get_keyword(unit,'NAXIS3',n3,status)
     if (.not. allocated(arr)) allocate(arr(n1,n2,n3))
     call fits_read_image(unit,arr,status)
     call fits_close(unit,status)
  else
     write(*,*) 'Error in reading the data file :', trim(fname),' status = ',status
  endif

  if (present(reduce_factor)) then
     if (reduce_factor > 1) then
        n1 = n1/reduce_factor
        n2 = n2/reduce_factor
        n3 = n3/reduce_factor
        do k=1,n3
        do j=1,n2
        do i=1,n1
           i2 = reduce_factor*i
           i1 = i2 - (reduce_factor-1)
           j2 = reduce_factor*j
           j1 = j2 - (reduce_factor-1)
           k2 = reduce_factor*k
           k1 = k2 - (reduce_factor-1)
           array(i,j,k) = sum(arr(i1:i2,j1:j2,k1:k2))/reduce_factor**3
        enddo
        enddo
        enddo
     else
        array(:,:,:) = arr(:,:,:)
     endif
  else
     array(:,:,:) = arr(:,:,:)
  endif
  if (allocated(arr)) deallocate(arr)

  if (present(centering)) then
     if (centering == 1 .or. centering == 2) then
        allocate(arr(n1,n2,n3))
        n1cen = n1/2
        n2cen = n2/2
        n3cen = n3/2
        if (centering == 1) then
           loc = minloc(array)
        else
           loc = maxloc(array)
        endif
        if (loc(1) < n1cen) then
           arr(n1cen-loc(1)+1:n1,:,:)   = array(1:loc(1)+n1cen,:,:)
           arr(1:n1cen-loc(1),:,:)      = array(loc(1)+n1cen+1:n1,:,:)
           array(:,:,:) = arr(:,:,:)
        else if (loc(1) > n1cen) then
           arr(1:3*n1cen-loc(1),:,:)    = array(loc(1)-n1cen+1:n1,:,:)
           arr(3*n1cen-loc(1)+1:n1,:,:) = array(1:loc(1)-n1cen,:,:)
           array(:,:,:) = arr(:,:,:)
        endif
        if (loc(2) < n2cen) then
           arr(:,n2cen-loc(2)+1:n2,:)   = array(:,1:loc(2)+n2cen,:)
           arr(:,1:n2cen-loc(2),:)      = array(:,loc(2)+n2cen+1:n2,:)
           array(:,:,:) = arr(:,:,:)
        else if (loc(2) > n2cen) then
           arr(:,1:3*n2cen-loc(2),:)    = array(:,loc(2)-n2cen+1:n2,:)
           arr(:,3*n2cen-loc(2)+1:n2,:) = array(:,1:loc(2)-n2cen,:)
           array(:,:,:) = arr(:,:,:)
        endif
        if (loc(3) < n3cen) then
           arr(:,:,n3cen-loc(3)+1:n3)   = array(:,:,1:loc(3)+n3cen)
           arr(:,:,1:n3cen-loc(3))      = array(:,:,loc(3)+n2cen+1:n3)
           array(:,:,:) = arr(:,:,:)
        else if (loc(3) > n3cen) then
           arr(:,:,1:3*n3cen-loc(3))    = array(:,:,loc(3)-n3cen+1:n3)
           arr(:,:,3*n3cen-loc(3)+1:n3) = array(:,:,1:loc(3)-n3cen)
           array(:,:,:) = arr(:,:,:)
        endif
        deallocate(arr)
        if (centering == 1) then
           loc = minloc(array)
        else
           loc = maxloc(array)
        endif
        if (loc(1) /= n1cen .or. loc(2) /= n2cen .or. loc(3) /= n3cen) then
           write(6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
           write(6,*) 'Something wrong. Please check the centering algorithm.'
           write(6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        endif
     endif
  endif

  end subroutine read_3D
!--------------------------------------------------------------
end module read_mod
