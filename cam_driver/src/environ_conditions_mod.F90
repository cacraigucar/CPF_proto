module environ_conditions_mod

!---------------------------------------------------------------------
! This module retrieves the enviromental conditions from a named file
!---------------------------------------------------------------------

  use machine,     only: rk => kind_phys
  use input_file,  only: input_file_type
  use input_slice, only: slice_type

  implicit none

  type, public :: environ_conditions
     private
     type(input_file_type) :: inputfile
     type(slice_type) :: slice
     real(rk) :: wghts(2) = 0._rk
     integer :: num_times
     real, allocatable :: times(:)
   contains
     procedure :: getvar => environ_conditions_getvar
     procedure :: getcol => environ_conditions_getcol
     procedure :: getlatcol => environ_conditions_getlatcol
     procedure :: getsrf => environ_conditions_getsrf
     procedure :: press_mid => environ_conditions_press_mid
     procedure :: press_int => environ_conditions_press_int
     procedure :: dtime => environ_conditions_dtime
     procedure :: ntimes => environ_conditions_ntimes
     procedure :: get_times => environ_conditions_times
     procedure :: nlevels => environ_conditions_nlevels
     procedure :: levnum => environ_conditions_levnum
     procedure :: environ_conditions_update_ndx
     procedure :: environ_conditions_update_flt
     generic   :: update => environ_conditions_update_ndx, environ_conditions_update_flt
  end type environ_conditions
  
contains


  function environ_conditions_create( infilepath, lat, lon, lev ) result(env_cond)
    use input_file, only: MAX_ATT_LEN
    
    character(len=*), intent(in) :: infilepath
    real, intent(in), optional ::  lat
    real, intent(in), optional ::  lon
    real, intent(in), optional :: lev
    
    type(environ_conditions), pointer :: env_cond
    character(len=MAX_ATT_LEN) :: units
    
    allocate(env_cond)
    
    call env_cond%inputfile%open( infilepath )
    if (present(lat)) then
       if (present(lon)) then
          if (present(lev)) then
             env_cond%slice = env_cond%inputfile%set_slice( beglat=lat,endlat=lat, beglon=lon,endlon=lon, beglev=lev,endlev=lev)
          else
             env_cond%slice = env_cond%inputfile%set_slice( beglat=lat,endlat=lat, beglon=lon,endlon=lon )
          endif
       else
          if (present(lev)) then
             env_cond%slice = env_cond%inputfile%set_slice( beglat=lat,endlat=lat,  beglev=lev,endlev=lev)
          else
             env_cond%slice = env_cond%inputfile%set_slice( beglat=lat,endlat=lat )
          endif
       endif
     else
       if (present(lon)) then
          if (present(lev)) then
             env_cond%slice = env_cond%inputfile%set_slice( beglon=lon,endlon=lon, beglev=lev,endlev=lev)
          else
             env_cond%slice = env_cond%inputfile%set_slice( beglon=lon,endlon=lon )
          endif
       else
          if (present(lev)) then
             env_cond%slice = env_cond%inputfile%set_slice( beglev=lev,endlev=lev)
          else
             env_cond%slice = env_cond%inputfile%set_slice(  )
          endif
       endif
    end if
     
    env_cond%slice%ntimes = 1

    env_cond%num_times = env_cond%inputfile%get_ntimes()
    allocate(env_cond%times(env_cond%num_times))

    units = env_cond%inputfile%get_units('time')
    if (index(units,'days')>0) then
       ! convert to seconds
       env_cond%times = 24._rk*3600._rk* env_cond%inputfile%get_times()
    else
       write(*,*) 'ERROR: Do not recognize time units in file: '//trim(infilepath)
       call abort()
    end if
 
  end function environ_conditions_create

  subroutine environ_conditions_update_ndx(this, record_num)
    class(environ_conditions), intent(inout) :: this
    integer, intent(in) :: record_num
    
    this%slice%begtime = record_num
  end subroutine environ_conditions_update_ndx

  subroutine environ_conditions_update_flt(this, time)
    class(environ_conditions), intent(inout) :: this
    real(rk), intent(in) :: time

    integer :: ndx

    if (time<=this%times(1)) then
       this%slice%begtime = 1
       this%slice%ntimes = 1
       this%wghts = 0._rk
    elseif (time>=this%times(this%num_times)) then
       this%slice%begtime = this%num_times
       this%slice%ntimes = 1      
       this%wghts = 0._rk
    else
       findtime: do ndx = 1,this%num_times-1
          if (this%times(ndx)>=time) then
             exit findtime
          end if
       end do findtime
       ndx=ndx-1
       this%slice%begtime = ndx
       this%slice%ntimes = 2
       this%wghts(2) = (time - this%times(ndx))/(this%times(ndx+1)-this%times(ndx))
       this%wghts(1) = 1._rk-this%wghts(2)
    end if
    
  end subroutine environ_conditions_update_flt

  function environ_conditions_getvar(this, var, default_value) result(thevalue)
    class(environ_conditions), intent(inout) :: this
    character(len=*), intent(in) :: var
    real(rk), optional, intent(in) :: default_value

    real(rk) :: thevalue

    real(rk), pointer :: data(:,:,:,:)
    
    data => this%inputfile%extract(var, this%slice, default_value)
    if (this%slice%ntimes == 2) then
       thevalue = this%wghts(1)*data(1,1,1,1) + this%wghts(2)*data(1,1,1,2)
    else
       thevalue = data(1,1,1,1)
    end if
    

  end function environ_conditions_getvar
  
  function environ_conditions_getlatcol(this, var, nlat, nlev) result(thelatcol)
    class(environ_conditions), intent(inout) :: this
    character(len=*), intent(in) :: var
    integer,intent(in) :: nlat,nlev
    real(rk) :: thelatcol(nlat,nlev)

    real(rk), pointer :: data(:,:,:,:)
    
    data => this%inputfile%extract(var, this%slice )
    if (this%slice%ntimes == 2) then
       thelatcol(1:nlat, 1:nlev) = this%wghts(1)*data(1,1:nlat,1:nlev,1) + this%wghts(2)*data(1,1:nlat,1:nlev,2)
    else
       thelatcol(1:nlat, 1:nlev) = data(1,1:nlat,1:nlev,1)
    end if

  end function environ_conditions_getlatcol

  function environ_conditions_getcol(this, var, nlev) result(thecol)
    class(environ_conditions), intent(inout) :: this
    character(len=*), intent(in) :: var
    integer,intent(in) :: nlev
    real(rk) :: thecol(nlev)

    real(rk), pointer :: data(:,:,:,:)
    
    data => this%inputfile%extract(var, this%slice )
    if (this%slice%ntimes == 2) then
       thecol(:) = this%wghts(1)*data(1,1,:,1) + this%wghts(2)*data(1,1,:,2)
    else
       thecol(:) = data(1,1,:,1)
    end if

  end function environ_conditions_getcol
  
  function environ_conditions_getsrf(this, var ) result(theval)
    class(environ_conditions), intent(inout) :: this
    character(len=*), intent(in) :: var
    real(rk) :: theval

    real, pointer :: data(:,:,:)
    
    data => this%inputfile%extract_srf(var, this%slice )
    if (this%slice%ntimes == 2) then
       theval = this%wghts(1)*data(1,1,1) + this%wghts(2)*data(1,1,2)
    else
       theval = data(1,1,1)
    end if
 
  end function environ_conditions_getsrf
  
  function environ_conditions_press_mid(this,nlev) result(press_mid)
    class(environ_conditions), intent(inout) :: this
    integer,intent(in) :: nlev
    real(rk) :: press_mid(nlev)

    real, pointer :: data_srf(:,:,:)
    real :: ps
    
    data_srf=>this%inputfile%extract_srf('PS',this%slice)
    if (this%slice%ntimes == 2) then
       ps = this%wghts(1)*data_srf(1,1,1) + this%wghts(2)*data_srf(1,1,2)
    else
       ps = data_srf(1,1,1)
    end if

    press_mid(:) = 1.e5*this%inputfile%get_hyam() + ps*this%inputfile%get_hybm()

  end function environ_conditions_press_mid
  
  function environ_conditions_press_int(this,nlev) result(press_int)
    class(environ_conditions), intent(inout) :: this
    integer,intent(in) :: nlev
    real(rk) :: press_int(nlev)

    real, pointer :: data_srf(:,:,:)
    real :: ps
    
    data_srf=>this%inputfile%extract_srf('PS',this%slice)
    if (this%slice%ntimes == 2) then
       ps = this%wghts(1)*data_srf(1,1,1) + this%wghts(2)*data_srf(1,1,2)
    else
       ps = data_srf(1,1,1)
    end if

    press_int(:) = 1.e5*this%inputfile%get_hyai() + ps*this%inputfile%get_hybi()

  end function environ_conditions_press_int
  
  function environ_conditions_ntimes(this) result(ntimes)
    class(environ_conditions), intent(in) :: this

    integer :: ntimes
    ntimes = this%num_times
    
  end function environ_conditions_ntimes
  function environ_conditions_times(this) result(times)
    class(environ_conditions), intent(in) :: this

    real(rk) :: times(this%num_times)
    times(:) = this%times(:)
    
  end function environ_conditions_times
  
  function environ_conditions_nlevels(this) result(nlevs)
    class(environ_conditions), intent(in) :: this

    integer :: nlevs
    nlevs = this%inputfile%get_nlevels()
    
  end function environ_conditions_nlevels
  
  function environ_conditions_dtime(this) result(dtime)
    class(environ_conditions), intent(in) :: this

    real :: dtime
    dtime = real(int(86400.0*this%inputfile%get_dtime())) ! convert to seconds
  end function environ_conditions_dtime

  function environ_conditions_levnum(this) result(levnum)
    class(environ_conditions), intent(in) :: this
    integer :: levnum

    levnum = this%slice%beglev

  end function environ_conditions_levnum
  
end module environ_conditions_mod
