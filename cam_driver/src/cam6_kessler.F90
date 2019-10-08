module cam_kessler_main

  use ccpp_kinds, only: kind_phys

  implicit none

contains

  subroutine cam_kessler_main_sub()

    use ppgrid,           only: pcols, pver, pverp
    use physics_types,    only: state, tend, physics_state, physics_type_alloc
    use CAM_ccpp_cap,     only: CAM_ccpp_physics_initialize
    use CAM_ccpp_cap,     only: CAM_ccpp_physics_timestep_initial
    use CAM_ccpp_cap,     only: CAM_ccpp_physics_run
    use CAM_ccpp_cap,     only: CAM_ccpp_physics_timestep_final
    use CAM_ccpp_cap,     only: CAM_ccpp_physics_finalize
    use ccpp_physics_api, only: ccpp_physics_suite_list
    use ccpp_physics_api, only: ccpp_physics_suite_part_list

    integer,            parameter   :: begchunk = 33
    integer,            parameter   :: endchunk = 33
    integer,            parameter   :: ncols = 1
    integer,            parameter   :: ntimes = 11

    integer                         :: i, j
    integer                         :: ierr
    integer                         :: lchnk
    integer                         :: ncol, nwrite
    character(len=20)               :: string
    character(len=512)              :: errmsg
    character(len=128), allocatable :: part_names(:)
    integer                         :: errflg

    ! Allocate the host variables
    call physics_type_alloc(state, tend, begchunk, endchunk, pcols)

    !  open(unit=50, file='../src/fort.50')
    !  read(50, fmt='(2i4)') icols, iver

    ! Use the suite information to setup the run
    call CAM_ccpp_physics_initialize('cam_kessler_test', begchunk, endchunk, errmsg, errflg)
    if (errflg /= 0) then
       write(6, *) trim(errmsg)
       stop
    end if

    ! loop over all time steps
    nwrite=0
    do j = 1, ntimes
!      read(50,fmt='(a10,i4)') string(1:19),nwrite
!      read(50,fmt='(a20,2i4,f20.13)') string(1:19),ncol, pver_in, ztodt
!      read(50,fmt='(a10,(e22.15))') string,rho(:ncol,:)
!      read(50,fmt='(a10,(e22.15))') string,z(:ncol,:)
!      read(50,fmt='(a10,(e22.15))') string,pk(:ncol,:)
!      read(50,fmt='(a10,(e22.15))') string,th(:ncol,:)
!      read(50,fmt='(a10,(e22.15))') string,qv(:ncol,:)
!      read(50,fmt='(a10,(e22.15))') string,qc(:ncol,:)
!      read(50,fmt='(a10,(e22.15))') string,qr(:ncol,:)
!      read(50,fmt='(a10,(e22.15))') string,precl(:ncol)
       call CAM_ccpp_physics_timestep_initial('cam_kessler_test', begchunk, endchunk, errmsg, errflg)
       do lchnk = begchunk, endchunk
          ncol = pcols
          call CAM_ccpp_physics_run('cam_kessler_test', 'physics', begchunk, endchunk, lchnk, ncol, errmsg, errflg)
          if (errflg /= 0) then
             write(6, *) trim(errmsg)
             call ccpp_physics_suite_part_list('cam_kessler_test', part_names, errmsg, errflg)
             write(6, *) 'Available suite parts are:'
             do nwrite = 1, size(part_names)
                write(6, *) trim(part_names(nwrite))
             end do
             stop
          end if
          write(6,*) 'At time step', j, 'in host model Temperature =', state(lchnk)%T(1, pver)
       end do
        ! write(51,'(a10,i4)') 'nwrite=',nwrite
        ! write(51,'(a20,2i4,f20.13)') 'ncol, pver, ztodt=',ncol, pver, ztodt
        ! write(51,fmt='(a10,(e22.15))') 'rho=',rho(:ncol,:)
        ! write(51,fmt='(a10,(e22.15))') 'z=',z(:ncol,:)
        ! write(51,fmt='(a10,(e22.15))') 'pk=',pk(:ncol,:)
        ! write(51,fmt='(a10,(e22.15))') 'th=',th(:ncol,:)
        ! write(51,fmt='(a10,(e22.15))') 'qv=',qv(:ncol,:)
        ! write(51,fmt='(a10,(e22.15))') 'qc=',qc(:ncol,:)
        ! write(51,fmt='(a10,(e22.15))') 'qr=',qr(:ncol,:)
        ! write(51,fmt='(a10,(e22.15))') 'precl=',precl(:ncol)
       nwrite = nwrite + 1

      ! Decrease the temperature every time step
!      state_host%Temperature = state_host%Temperature - 100._kind_phys
       call CAM_ccpp_physics_timestep_final('cam_kessler_test', begchunk, endchunk, errmsg, errflg)

    end do


    call CAM_ccpp_physics_finalize('cam_kessler_test', begchunk, endchunk, errmsg, errflg)
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
