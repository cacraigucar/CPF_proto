module cam_kessler_main

  use machine, only: kind_phys
  use const_props_mod,        only: const_props_type
  use environ_conditions_mod, only: environ_conditions_create, environ_conditions

  implicit none

contains

  subroutine cam_kessler_main_sub()

    use ppgrid,           only: pcols, pver, pverp, pcnst
    use physics_types,    only: state, tend, physics_state, physics_type_alloc
    use CAM_ccpp_cap,     only: CAM_ccpp_physics_initialize
    use CAM_ccpp_cap,     only: CAM_ccpp_physics_timestep_initial
    use CAM_ccpp_cap,     only: CAM_ccpp_physics_run
    use CAM_ccpp_cap,     only: CAM_ccpp_physics_timestep_final
    use CAM_ccpp_cap,     only: CAM_ccpp_physics_finalize
    use ccpp_physics_api, only: ccpp_physics_suite_list
    use ccpp_physics_api, only: ccpp_physics_suite_part_list

    integer,            parameter   :: begchunk = 33 ! Not needed for CAM7
    integer,            parameter   :: endchunk = 33 ! Not needed for CAM7
    integer,            parameter   :: nlat = 192
    integer                         :: nkessler_loop  = 3

    integer                         :: i, k, rk, ind_kessler_loop
    integer                         :: ierr
    integer                         :: col_start, col_end
    integer                         :: ncol, nwrite, pver_in, nwrite_in, nstep, ncnst
    real(kind_phys)                 :: ztodt
    real(kind_phys)                 :: precl(pcols)
    real(kind_phys)                 :: scratch(pcols,pver)
    character(len=20)               :: string
    character(len=512)              :: errmsg
    character(len=128), allocatable :: part_names(:)
    integer                         :: errflg
    real(kind_phys)                 :: pmiddry_top2bot(pcols,pver)
    real(kind_phys)                 :: pint_top2bot(pcols,pverp)
    real(kind_phys)                 :: pmid_top2bot(pcols,pver)
    real(kind_phys)                 :: pdel_top2bot(pcols,pver)
    real(kind_phys)                 :: rpdel_top2bot(pcols,pver)
    real(kind_phys)                 :: pdeldry_top2bot(pcols,pver)
    real(kind_phys)                 :: zi_top2bot(pcols,pverp)
    real(kind_phys)                 :: zm_top2bot(pcols,pver)
    real(kind_phys)                 :: exner_top2bot(pcols,pver)
    real(kind_phys)                 :: t_top2bot(pcols,pver)
    real(kind_phys)                 :: s_top2bot(pcols,pver)
    real(kind_phys)                 :: q_top2bot(pcols,pver,pcnst)
    real(kind_phys)                 :: lnpint_top2bot(pcols,pverp)
    real(kind_phys)                 :: lnpmid_top2bot(pcols,pver)
    real(kind_phys)                 :: ttend_top2bot(pcols,pver)

    real(kind=kind_phys), allocatable :: wghts(:)
    integer :: nlevels, photo_lev

    character(len=255) :: model_name


    character(len=120) :: jsonfile  ! TEMPORARY VARIABLE
    integer                         :: njRxt, nkRxt, nTotRxt
    type(environ_conditions),pointer :: theEnvConds => null()
    type(environ_conditions),pointer :: colEnvConds => null()
    real(kind_phys),allocatable     :: vmr(:)
    real(kind_phys)                 :: dt 
    real(kind_phys)                 :: total_dens 
    integer                         :: file_ntimes
    real(kind_phys), allocatable    :: file_times(:)
    type(const_props_type), allocatable  :: cnst_info(:) 
    character(len=16)  :: cnst_name
  character(len=120) :: env_conds_file = '../data/env_conditions.nc'
  real, parameter :: NOT_SET = -huge(1.0)
  real :: env_lon = 20._kind_phys
  real :: env_lev = 100._kind_phys ! mbar

  real(kind_phys) :: TimeStart, TimeEnd, TimeCurrent
  integer :: n

  ! If true, run the basic kessler test which reads in fort.60 and writes out
  ! fort.61 for comparison with the corresponding CESM run 
  logical :: test_kessler_cesm = .false.


    ! Allocate the host variables
    call physics_type_alloc(state, tend, pcols)


    ! counter for output
    nstep=0

    ! Read the one time environmental conditions information and allocate arrays
    ! if using MICM along with kessler
    if (.not. test_kessler_cesm)  then
       theEnvConds => environ_conditions_create( env_conds_file, lon=10., lev=env_lev )

       file_ntimes= theEnvConds%ntimes()
       allocate(file_times(file_ntimes))
       file_times = theEnvConds%get_times()
! Not sure how to read in these variables - does not affect the run, but will
! certainly not produce correct results without these variables
!        precl(:)  = colEnvConds%getcol('PRECL', pcols)
!        state%phis(:)  = colEnvConds%getcol('PHIS', pcols)
       precl = 0._kind_phys
       state%phis = 0._kind_phys

       colEnvConds => environ_conditions_create( env_conds_file, lon=env_lon )
       nlevels = colEnvConds%nlevels()
       photo_lev = theEnvConds%levnum()

       allocate (vmr(ncnst))
    else
       allocate (vmr(1))  ! array just needs to be allocated as it is in the CCPP calling lists for now
    end if

    ! Use the suite information to setup the run
    call CAM_ccpp_physics_initialize('cam_kessler_test', precl, vmr, total_dens, errmsg, errflg)
    if (errflg /= 0) then
       write(6, *) trim(errmsg)
       stop
    end if
    call CAM_ccpp_physics_initialize('cam_micm', precl, vmr, total_dens, errmsg, errflg)
    if (errflg /= 0) then
       write(6, *) trim(errmsg)
       stop
    end if

    dt = 150.
    TimeStart = 0._kind_phys
    TimeEnd = TimeStart + dt
    TimeCurrent = TimeStart

    do ind_kessler_loop  = 1, nkessler_loop

     ! test_kessler_cesm assumes that fort.60 has been linked from the appropriate directory - as of March 2019 this is
     ! /home/cacraig/kessler00_004_data

     if (test_kessler_cesm)  then

       ncol = pcols
       read(60,fmt='(a10,i4)') string(1:8),nwrite_in
       read(60,fmt='(a20,2i4,f20.13)') string(1:19),ncol, pver_in, ztodt
       read(60,fmt='(a20,(e25.18))') string(1:20),pmiddry_top2bot(:ncol,:pver_in)
       read(60,fmt='(a20,(e25.18))') string(1:20),pint_top2bot(:ncol,:pverp)
       read(60,fmt='(a20,(e25.18))') string(1:20),pmid_top2bot(:ncol,:pver_in)
       read(60,fmt='(a20,(e25.18))') string,pdel_top2bot(:ncol,:pver_in)
       read(60,fmt='(a20,(e25.18))') string,rpdel_top2bot(:ncol,:pver_in)
       read(60,fmt='(a20,(e25.18))') string,pdeldry_top2bot(:ncol,:pver_in)
       read(60,fmt='(a20,(e25.18))') string,exner_top2bot(:ncol,:pver_in)
       read(60,fmt='(a20,(e25.18))') string,t_top2bot(:ncol,:pver_in)
       read(60,fmt='(a20,(e25.18))') string,s_top2bot(:ncol,:pver_in)
       read(60,fmt='(a20,(e25.18))') string,q_top2bot(:ncol,:pver_in,1)
       read(60,fmt='(a20,(e25.18))') string,q_top2bot(:ncol,:pver_in,2)
       read(60,fmt='(a20,(e25.18))') string,q_top2bot(:ncol,:pver_in,3)
       read(60,fmt='(a20,(e25.18))') string,zi_top2bot(:ncol,:pverp)
       read(60,fmt='(a20,(e25.18))') string,zm_top2bot(:ncol,:pver)
       read(60,fmt='(a20,(e25.18))') string,lnpint_top2bot(:ncol,:pverp)
       read(60,fmt='(a20,(e25.18))') string,lnpmid_top2bot(:ncol,:pver)
       read(60,fmt='(a20,(e25.18))') string,precl(:ncol)
       read(60,fmt='(a20,(e25.18))') string,state%phis(:ncol)
       read(60,fmt='(a20,(e25.18))') string,ttend_top2bot(:ncol,:pver_in)
     else
       ! For testing purposes (and to limit the test size) set the number of
       ! columns to the number of latitudes.  pcols is set to the number of
       ! latitudes in the ppgrid.F90-with_micm file)
       ncol = pcols

       call colEnvConds%update(TimeCurrent)
       do n = 1,ncnst
          call cnst_info(n)%print()
          cnst_name = cnst_info(n)%get_name()
          if (cnst_name == 'N2') then
             vmr(n) = theEnvConds%getvar(cnst_name,default_value=0.79_kind_phys)
          else if (cnst_name == 'O2') then
             vmr(n) = theEnvConds%getvar(cnst_name,default_value=0.21_kind_phys)
          else
             vmr(n) = theEnvConds%getvar(cnst_name,default_value=0.00_kind_phys)
          end if

          write(*,fmt="(' cnst name : ',a20,' init value : ',e13.6)") cnst_name, vmr(n)
          if (allocated(wghts) .and. cnst_name == 'CL2') then
             wghts(n) = 2._kind_phys
          end if

       enddo


! Don't know how to read the pverp variables - does not appear to be needed to
! run, but certainly won't give correct answers
!        zi_top2bot(:,:)  = colEnvConds%getlatcol('ZI',nlat, pverp)
!        lnpint_top2bot(:,:)  = colEnvConds%getlatcol('LNPINT',nlat, pverp)
!        pint_top2bot(:,:) = colEnvConds%getlatcol('PINT',nlat, pverp)
        zi_top2bot = 0._kind_phys
        lnpint_top2bot = 0._kind_phys
        pint_top2bot = 0._kind_phys

       ! To make the size more reasonable, read all latitudes at only one longitude location
        pmiddry_top2bot(:,:)=colEnvConds%getlatcol('PMIDDRY',nlat, pver)
        pmid_top2bot(:,:) = colEnvConds%getlatcol('PMID',nlat, pver)
        pdel_top2bot(:,:) = colEnvConds%getlatcol('PDEL',nlat, pver)
        rpdel_top2bot(:,:) = colEnvConds%getlatcol('RPDEL',nlat, pver)
        pdeldry_top2bot(:,:) = colEnvConds%getlatcol('PDELDRY',nlat, pver)
        exner_top2bot(:,:)  = colEnvConds%getlatcol('EXNER',nlat, pver)
        t_top2bot(:,:) = colEnvConds%getlatcol('T',nlat, pver)
        s_top2bot(:,:)  = colEnvConds%getlatcol('S',nlat, pver)
        q_top2bot(:,:,1)  = colEnvConds%getlatcol('Q',nlat, pver)
        q_top2bot(:,:,2)  = colEnvConds%getlatcol('CLDLIQ',nlat, pver)
        q_top2bot(:,:,3)  = colEnvConds%getlatcol('RAINQM',nlat, pver)
        zm_top2bot(:,:)  = colEnvConds%getlatcol('ZM',nlat, pver)
        lnpmid_top2bot(:,:)  = colEnvConds%getlatcol('LNPMID',nlat, pver)
        ttend_top2bot(:,:)  = colEnvConds%getlatcol('PTTEND',nlat, pver)

     end if


       ! Need to swap the bottom and top
       do k=1,pver
         rk= pver - k +1 
         state%pmiddry(:ncol,rk) = pmiddry_top2bot(:ncol,k)
         state%pmid(:ncol,rk)    = pmid_top2bot(:ncol,k)
         state%pdel(:ncol,rk)    = pdel_top2bot(:ncol,k)
         state%rpdel(:ncol,rk)   = rpdel_top2bot(:ncol,k)
         state%pdeldry(:ncol,rk) = pdeldry_top2bot(:ncol,k)
         state%exner(:ncol,rk)   = exner_top2bot(:ncol,k)
         state%t(:ncol,rk)       = t_top2bot(:ncol,k)
         state%s(:ncol,rk)       = s_top2bot(:ncol,k)
         state%zm(:ncol,rk)      = zm_top2bot(:ncol,k)
         state%lnpmid(:ncol,rk)  = lnpmid_top2bot(:ncol,k)
         state%q(:ncol,rk,1)     = q_top2bot(:ncol,k,1)
         state%q(:ncol,rk,2)     = q_top2bot(:ncol,k,2)
         state%q(:ncol,rk,3)     = q_top2bot(:ncol,k,3)
       end do

       do k=1,pverp
         rk= pverp - k +1 
         state%pint(:ncol,rk)    = pint_top2bot(:ncol,k)
         state%lnpint(:ncol,rk)  = lnpint_top2bot(:ncol,k)
         state%zi(:ncol,rk)      = zi_top2bot(:ncol,k)
       end do



!        ! This is a temporary read until these values are provided by the
!      ! configuratore either within metadata or via the namelist
!#include "model_name.inc"
!      jsonfile = '../../MICM_chemistry/generated/'//trim(model_name)//'/molec_info.json'
!      call json_loader_init( jsonfile, ncnst, nkrxt, nkrxt )
!      nTotrxt = nkrxt + nkrxt

      if (model_name == 'terminator') then
         allocate(wghts(ncnst))
         wghts(:) = 1._kind_phys
      endif

       ! Initialize the timestep
       call CAM_ccpp_physics_timestep_initial('cam_kessler_test', precl, vmr, total_dens, errmsg, errflg)
       if (errflg /= 0) then
          write(6, *) trim(errmsg)
          stop
       end if
        call CAM_ccpp_physics_timestep_initial('cam_micm', precl, vmr, total_dens, errmsg, errflg)
        if (errflg /= 0) then
           write(6, *) trim(errmsg)
           stop
        end if

       col_start = 1
       col_end = ncol


       ! Initialize the total tendency after the timestep initialization
       do k=1,pver
         rk= pver - k +1 
         tend%dtdt(:ncol,rk)     = ttend_top2bot(:ncol,k)
       end do

       call CAM_ccpp_physics_run('cam_kessler_test', 'physics', col_start, col_end, precl, vmr, total_dens, errmsg, errflg)
       if (errflg /= 0) then
          write(6, *) trim(errmsg)
          call ccpp_physics_suite_part_list('cam_kessler_test', part_names, errmsg, errflg)
          write(6, *) 'Available suite parts are:'
          do nwrite = 1, size(part_names)
             write(6, *) trim(part_names(nwrite))
          end do
          stop
       end if

       call CAM_ccpp_physics_run('cam_micm', 'physics', col_start, col_end, precl, vmr, total_dens, errmsg, errflg)
       if (errflg /= 0) then
          write(6, *) trim(errmsg)
          call ccpp_physics_suite_part_list('cam_kessler_test', part_names, errmsg, errflg)
          write(6, *) 'Available suite parts are:'
          do nwrite = 1, size(part_names)
             write(6, *) trim(part_names(nwrite))
          end do
          stop
       end if


       write(6,*) 'At time step', ind_kessler_loop, 'in host model Temperature =', state%T(8, :pver)


       call CAM_ccpp_physics_timestep_final('cam_kessler_test', precl, vmr, total_dens, errmsg, errflg)
       if (errflg /= 0) then
          write(6, *) trim(errmsg)
          stop
       end if

        call CAM_ccpp_physics_timestep_final('cam_micm', precl, vmr, total_dens, errmsg, errflg)
        if (errflg /= 0) then
           write(6, *) trim(errmsg)
           stop
        end if

         write(61,'(a10,i4)') 'nstep=',nstep
         write(61,'(a20,2i4,f20.13)') 'ncol, pver, ztodt=',ncol, pver, ztodt
         write(61,fmt='(a10,(e25.18))') 'pmiddry=',state%pmiddry(:ncol,pver:1:-1)
         write(61,fmt='(a10,(e25.18))') 'pint=',state%pint(:ncol,pverp:1:-1)
         write(61,fmt='(a10,(e25.18))') 'pmid=',state%pmid(:ncol,pver:1:-1)
         write(61,fmt='(a10,(e25.18))') 'pdel=',state%pdel(:ncol,pver:1:-1)
         write(61,fmt='(a10,(e25.18))') 'rpdel=',state%rpdel(:ncol,pver:1:-1)
         write(61,fmt='(a10,(e25.18))') 'pdeldry=',state%pdeldry(:ncol,pver:1:-1)
         write(61,fmt='(a10,(e25.18))') 'exner=',state%exner(:ncol,pver:1:-1)
         write(61,fmt='(a10,(e25.18))') 'state%t=',state%t(:ncol,pver:1:-1)
         write(61,fmt='(a10,(e25.18))') 'state%s=',state%s(:ncol,pver:1:-1)
         write(61,fmt='(a10,(e22.15))') 'qv=',state%q(:ncol,pver:1:-1,1)
         write(61,fmt='(a10,(e25.18))') 'qc=',state%q(:ncol,pver:1:-1,2)
         write(61,fmt='(a10,(e25.18))') 'qr=',state%q(:ncol,pver:1:-1,3)
         write(61,fmt='(a20,(e25.18))') 'zi=',state%zi(:ncol,pverp:1:-1)
         write(61,fmt='(a20,(e25.18))') 'zm=',state%zm(:ncol,pver:1:-1)
         write(61,fmt='(a20,(e25.18))') 'lnpint=',state%lnpint(:ncol,pverp:1:-1)
         write(61,fmt='(a20,(e25.18))') 'lnpmid=',state%lnpmid(:ncol,pver:1:-1)
         write(61,fmt='(a10,(e25.18))') 'precl=',precl(:ncol)
         write(61,fmt='(a10,(e25.18))') 'phis=',state%phis(:ncol)
         write(61,fmt='(a10,(e25.18))') 'tend%dtdt=',tend%dtdt(:ncol,pver:1:-1)


       ! Update the timestep

       nstep = nstep + 1
       TimeCurrent = TimeStart + ind_kessler_loop*dt


    end do


    call CAM_ccpp_physics_finalize('cam_kessler_test', precl, vmr, total_dens, errmsg, errflg)
    if (errflg /= 0) then
       write(6, *) trim(errmsg)
       write(6,'(a)') 'An error occurred in ccpp_timestep_final, Exiting...'
       stop
    end if

    call CAM_ccpp_physics_finalize('cam_micm', precl, vmr, total_dens, errmsg, errflg)
    if (errflg /= 0) then
       write(6, *) trim(errmsg)
       write(6,'(a)') 'An error occurred in ccpp_timestep_final, Exiting...'
       stop
    end if


  end subroutine cam_kessler_main_sub

end module cam_kessler_main

!> \brief Main SCM program that calls the main SCM subroutine
!!
!! The Doxygen documentation system cannot handle in-body comments in Fortran main programs, so the "main" program was put in the
!! subroutine \ref cam_kessler_main_sub above.
program cam_kessler
  use cam_kessler_main, only: cam_kessler_main_sub
  call cam_kessler_main_sub()
end program cam_kessler
