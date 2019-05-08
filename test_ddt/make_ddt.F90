!Hello demonstration parameterization
!

MODULE make_ddt

  USE ccpp_kinds, ONLY: kind_phys

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: make_ddt_init
  PUBLIC :: make_ddt_run
  PUBLIC :: make_ddt_finalize
  PUBLIC :: time_and_temp_type
  
  type time_and_temp_type
   real(kind_phys)  :: timestep
   REAL(kind_phys)  :: temp_layer(ncol, lev)
   REAL(kind_phys)  :: temp_level(ncol, ilev)
  end type
    

CONTAINS

!> \section arg_table_make_ddt_run  Argument Table
!! \htmlinclude arg_table_make_ddt_run.html
!!
  SUBROUTINE make_ddt_run(ncol, lev, ilev, time_and_temp, timestep, temp_level,          &
       temp_layer, errmsg, errflg)
!----------------------------------------------------------------
   IMPLICIT NONE
!----------------------------------------------------------------

   integer,            intent(in)    :: ncol, lev, ilev
   type(time_and_temp_type) , intent(inout) :: time_and_temp
   REAL(kind_phys),    intent(inout) :: temp_level(ncol, ilev)
   real(kind_phys),    intent(in)    :: timestep
   REAL(kind_phys),    INTENT(out)   :: temp_layer(ncol, lev)
   character(len=512), intent(out)   :: errmsg
   integer,            intent(out)   :: errflg
!----------------------------------------------------------------

   integer :: col_index
   integer :: lev_index

    errmsg = ''
    errflg = 0
 
    time_and_temp%temp_level(:,:) = temp_level(:,:)
    time_and_temp%timestep = timestep
    time_and_temp%temp_layer = 0._kind_phys

  END SUBROUTINE make_ddt_run

!> \section arg_table_make_ddt_init  Argument Table
!! \htmlinclude arg_table_make_ddt_init.html
!!
  subroutine make_ddt_init (errmsg, errflg)

    character(len=512),      intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

    ! This routine currently does nothing

    errmsg = ''
    errflg = 0

  end subroutine make_ddt_init

!> \section arg_table_make_ddt_finalize  Argument Table
!! \htmlinclude arg_table_make_ddt_finalize.html
!!
  subroutine make_ddt_finalize (errmsg, errflg)

    character(len=512),      intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

    ! This routine currently does nothing

    errmsg = ''
    errflg = 0

  end subroutine make_ddt_finalize

END MODULE make_ddt
