!Hello demonstration parameterization
!

MODULE kinetics_solve

  USE ccpp_kinds, ONLY: kind_phys
  use vmr_ddt_mod,   only: vmr_type

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: kinetics_solve_init
  PUBLIC :: kinetics_solve_run
  PUBLIC :: kinetics_solve_finalize
  

CONTAINS

!> \section arg_table_kinetics_solve_run  Argument Table
!! \htmlinclude arg_table_kinetics_solve_run.html
!!
  SUBROUTINE kinetics_solve_run(nbox, nlev,  vmr,        &
       errmsg, errflg)

  use kinetics_sub, only: factor

!----------------------------------------------------------------
   IMPLICIT NONE
!----------------------------------------------------------------

   integer,            intent(in)    :: nbox, nlev
   type(vmr_type),     intent(inout) :: vmr
   character(len=512), intent(out)   :: errmsg
   integer,            intent(out)   :: errflg
 
   integer :: i
!----------------------------------------------------------------
   

    errmsg = ''
    errflg = 0
 

    ! NOTE -- This is prototyping one approach to passing a large number of
    ! chemical VMR values and is the predecssor for adding in methods and maybe
    ! nesting DDTs (especially for aerosols)
    
    do i = 1,vmr%nvmr
!      vmr%vmr_array(:,:,i) = vmr%vmr_array(:,:,i)/2.0_kind_phys
      vmr%vmr_array(:,:,i) = vmr%vmr_array(:,:,i)/factor
    end do

  END SUBROUTINE kinetics_solve_run

!> \section arg_table_kinetics_solve_init  Argument Table
!! \htmlinclude arg_table_kinetics_solve_init.html
!!
  subroutine kinetics_solve_init (errmsg, errflg)

    character(len=512),      intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

    ! This routine currently does nothing

    errmsg = ''
    errflg = 0

  end subroutine kinetics_solve_init

!> \section arg_table_kinetics_solve_finalize  Argument Table
!! \htmlinclude arg_table_kinetics_solve_finalize.html
!!
  subroutine kinetics_solve_finalize (errmsg, errflg)

    character(len=512),      intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

    ! This routine currently does nothing

    errmsg = ''
    errflg = 0

  end subroutine kinetics_solve_finalize

END MODULE kinetics_solve
