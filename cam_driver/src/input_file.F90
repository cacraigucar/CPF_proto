module input_file

  use machine,only: rk => kind_phys

  use netcdf, only: nf90_open, nf90_nowrite, nf90_noerr, nf90_inq_dimid, nf90_inquire_dimension
  use netcdf, only: nf90_inq_varid, nf90_get_var, nf90_inquire_attribute, nf90_get_att
  
  implicit none

  private

  integer, public, parameter :: MAX_ATT_LEN=128
  type, public :: input_file_type

     private

     integer :: ncid

     real, allocatable :: lons(:)
     real, allocatable :: lats(:)
     real, allocatable :: levs(:)
     real, allocatable :: times(:)

     real, allocatable :: hyam(:)
     real, allocatable :: hybm(:)
     real, allocatable :: hyai(:)
     real, allocatable :: hybi(:)
     real, allocatable :: zi(:)
     real, allocatable :: lnpint(:)

     integer :: nlons
     integer :: nlats
     integer :: nlevs
     integer :: ntimes

   contains
     procedure :: open => input_file_open
     procedure :: set_slice => input_file_slice
     procedure :: extract => input_file_extract_slice
     procedure :: extract_srf => input_file_extract_slice3d
     procedure :: get_ntimes => input_file_get_ntimes
     procedure :: get_nlevels => input_file_get_nlevels
     procedure :: get_time => input_file_get_time
     procedure :: get_dtime => input_file_get_dtime
     procedure :: get_times => input_file_get_times
     procedure :: get_hyam=>input_file_get_hyam
     procedure :: get_hybm=>input_file_get_hybm
     procedure :: get_hyai=>input_file_get_hyai
     procedure :: get_hybi=>input_file_get_hybi
     procedure :: get_zi=>input_file_get_zi
     procedure :: get_lnpint=>input_file_get_lnpint
     procedure :: get_units=>input_file_get_units
  end type input_file_type

contains

  function input_file_get_hyam(this) result(x)
    class(input_file_type), intent(in) :: this
    real :: x(this%nlevs)
    x(:) = this%hyam(:)
  end function input_file_get_hyam
  function input_file_get_hybm(this) result(x)
    class(input_file_type), intent(in) :: this
    real :: x(this%nlevs)
    x(:) = this%hybm(:)
  end function input_file_get_hybm
  
  function input_file_get_hyai(this) result(x)
    class(input_file_type), intent(in) :: this
    real :: x(this%nlevs+1)
    x(:) = this%hyai(:)
  end function input_file_get_hyai
  function input_file_get_zi(this) result(x)
    class(input_file_type), intent(in) :: this
    real :: x(this%nlevs+1)
    x(:) = this%zi(:)
  end function input_file_get_zi
  function input_file_get_lnpint(this) result(x)
    class(input_file_type), intent(in) :: this
    real :: x(this%nlevs+1)
    x(:) = this%lnpint(:)
  end function input_file_get_lnpint
  function input_file_get_hybi(this) result(x)
    class(input_file_type), intent(in) :: this
    real :: x(this%nlevs+1)
    x(:) = this%hybi(:)
  end function input_file_get_hybi
  
  real function input_file_get_time(this, n)
    class(input_file_type), intent(in) :: this
    integer, intent(in) :: n
    input_file_get_time = this%times(n)
  end function input_file_get_time
  
  real function input_file_get_dtime(this)
    class(input_file_type), intent(in) :: this
    input_file_get_dtime = this%times(2) - this%times(1)
  end function input_file_get_dtime
  
  function input_file_get_times(this) result(x)
    class(input_file_type), intent(in) :: this
    real(rk) :: x(this%ntimes)
    x(:) = this%times(:)
  end function input_file_get_times
  
  integer function input_file_get_ntimes(this)
    class(input_file_type), intent(in) :: this
    input_file_get_ntimes = this%ntimes
  end function input_file_get_ntimes

  integer function input_file_get_nlevels(this)
    class(input_file_type), intent(in) :: this
    input_file_get_nlevels = this%nlevs
  end function input_file_get_nlevels

  subroutine input_file_open(this, filename )

    class(input_file_type), intent(inout) :: this

    character(len=*), intent(in) :: filename
    
    integer :: status
    
    status = nf90_open(filename, nf90_nowrite, this%ncid)
    if(status /= nf90_noerr) call handle_err(status)

    call set_coordinate(this%ncid, 'lat',  this%nlats,  this%lats )
    call set_coordinate(this%ncid, 'lon',  this%nlons,  this%lons )
    call set_coordinate(this%ncid, 'lev',  this%nlevs,  this%levs )
    call set_coordinate(this%ncid, 'time', this%ntimes, this%times )

    call set_hybrid_coord(this%ncid, 'hyam', this%nlevs, this%hyam)
    call set_hybrid_coord(this%ncid, 'hybm', this%nlevs, this%hybm)    
    call set_hybrid_coord(this%ncid, 'hyai', this%nlevs+1, this%hyai)
    call set_hybrid_coord(this%ncid, 'hybi', this%nlevs+1, this%hybi)
    
  contains

    subroutine set_hybrid_coord(fileid,coordname,ndim, coordvar)
      integer, intent(in) :: fileid
      character(len=*), intent(in) :: coordname
      integer, intent(in) :: ndim
      real, allocatable, intent(out) :: coordvar(:)
      
      integer :: varid

      allocate(coordvar(ndim))

      status = nf90_inq_varid(fileid, coordname, varid)
      if(status /= nf90_noerr) call handle_err(status)
      status = nf90_get_var(fileid, varid, coordvar)
      if(status /= nf90_noerr) call handle_err(status)
    end subroutine set_hybrid_coord
    

    subroutine set_coordinate(fileid, dimname, ndim, coordvar )

      integer, intent(in) :: fileid
      character(len=*), intent(in) :: dimname
      integer, intent(out) :: ndim
      real, allocatable, intent(out) :: coordvar(:)

      integer :: dimid, varid

      status = nf90_inq_dimid(fileid, dimname, dimid)
      if (status /= nf90_noerr) call handle_err(status)
      
      status = nf90_inquire_dimension(fileid, dimid, len = ndim)
      if (status /= nf90_noerr) call handle_err(status)

      allocate(coordvar(ndim))

      status = nf90_inq_varid(fileid, dimname, varid)
      if(status /= nf90_noerr) call handle_err(status)

      status = nf90_get_var(fileid, varid, coordvar)
      if(status /= nf90_noerr) call handle_err(status)
 
    end subroutine set_coordinate
  
  end subroutine input_file_open

  function input_file_slice(this, beglon,endlon, beglat,endlat, beglev,endlev, begtime,endtime) result( slice )
    use input_slice, only : slice_type

    class(input_file_type), intent(in) :: this
    real, optional, intent(in) :: beglon,endlon, beglat,endlat, beglev,endlev
    real, optional, intent(in) :: begtime,endtime
    
    type(slice_type) :: slice

    integer :: i

    if (present(beglon) .and. present(endlon)) then
       call get_slice( this%nlons,this%lons,  beglon,endlon, slice%beglon, slice%nlons )
    else
       slice%beglon = 1
       slice%nlons = this%nlons
    endif

    if (present(beglat) .and. present(endlat)) then
       call get_slice( this%nlats,this%lats,  beglat,endlat, slice%beglat, slice%nlats )
    else
       slice%beglat = 1
       slice%nlats = this%nlats
    endif
   
    if (present(beglev) .and. present(endlev)) then
       call get_slice( this%nlevs,this%levs,  beglev,endlev, slice%beglev, slice%nlevs )
    else
       slice%beglev = 1
       slice%nlevs = this%nlevs
    endif
    
    if (present(begtime) .and. present(endtime)) then
       call get_slice( this%ntimes,this%times,  begtime,endtime, slice%begtime, slice%ntimes )
    else
       slice%begtime = 1
       slice%ntimes = 1
    end if
    
  contains
    
    subroutine get_slice(ncoord,coords, begcoord,endcoord, begslice, nslice )

      integer, intent(in) :: ncoord
      real, intent(in) :: coords(:)
      real, intent(in) :: begcoord,endcoord
      integer, intent(out) :: begslice, nslice

      nslice = 0
      begslice = -1
      findslice: do i = 1,ncoord
         if (coords(i)>= begcoord .and. begslice<1) then
            begslice = i
            nslice = nslice+1
        elseif (coords(i)>= begcoord .and. coords(i)<=endcoord) then
            nslice = nslice+1
         endif
         if (coords(i) >= endcoord) then
            exit findslice
         endif
      end do findslice

      if (begslice<1) then
         if (begcoord<coords(1)) then
            begslice=1
            nslice=1
         elseif (begcoord>coords(ncoord)) then
            begslice=ncoord
            nslice=1
         end if
      end if

    end subroutine get_slice
  
  end function input_file_slice

  function input_file_get_units(this, varname) result(units)
    class(input_file_type), intent(inout) :: this

    character(len=*), intent(in) :: varname

    character(len=MAX_ATT_LEN) :: units

    integer :: status, varid, length
    
    status = nf90_inq_varid(this%ncid, varname, varid)
    if(status /= nf90_noerr) call handle_err(status)
    status = nf90_inquire_attribute(this%ncid, varid, "units", len = length)
    if(status /= nf90_noerr) call handle_err(status)
    if (length>MAX_ATT_LEN) then
       write(*,*) 'ERROR: input_file_get_units: units length too long on var '//trim(varname)//'... Length = ',length
       stop
    end if
    
    status = nf90_get_att(this%ncid, varid, "units", units)
    if(status /= nf90_noerr) call handle_err(status)
   
  end function input_file_get_units

  function input_file_extract_slice(this, varname, slice, default_value) result(data)
    use input_slice, only : slice_type

    class(input_file_type), intent(inout) :: this
    
    character(len=*), intent(in) :: varname
    type(slice_type), intent(in) :: slice
    real(rk), optional, intent(in) :: default_value
    
    real(rk), pointer :: data(:,:,:,:)

    integer :: varid, status

    allocate(data(slice%nlons,slice%nlats,slice%nlevs,slice%ntimes))
    data = -1.e36_rk
    
    status = nf90_inq_varid(this%ncid, varname, varid)
    
    if (status == nf90_noerr) then
       status = nf90_get_var(this%ncid, varid, data, &
            start = (/ slice%beglon, slice%beglat, slice%beglev, slice%begtime /),     &
            count = (/ slice%nlons,  slice%nlats,  slice%nlevs,  slice%ntimes /))
       if(status /= nf90_noerr) call handle_err(status)
    elseif (present(default_value)) then
       data = default_value
    else
       call handle_err(status)
    end if

  end function input_file_extract_slice

  function input_file_extract_slice3d(this, varname, slice ) result(data)
    use input_slice, only : slice_type

    class(input_file_type), intent(inout) :: this
    
    character(len=*), intent(in) :: varname
    type(slice_type), intent(in) :: slice
    
    real, pointer :: data(:,:,:)

    integer :: varid, status

    allocate(data(slice%nlons,slice%nlats,slice%ntimes))
    data = -1.e36
    
    status = nf90_inq_varid(this%ncid, varname, varid)
    if(status /= nf90_noerr) call handle_err(status)

    status = nf90_get_var(this%ncid, varid, data, &
                            start = (/ slice%beglon, slice%beglat, slice%begtime /),     &
                            count = (/ slice%nlons,  slice%nlats,  slice%ntimes /))

    if(status /= nf90_noerr) call handle_err(status)
    
  end function input_file_extract_slice3d

  subroutine handle_err(status)

    integer, intent(in) :: status

    write(*,*) 'ERROR input_file; status = ', status
    call exit(status)
    
  end subroutine handle_err
  
end module input_file
