!Hello demonstration parameterization
!

MODULE finish_ddt

  USE ccpp_kinds, ONLY: kind_phys
  use vmr_ddt_mod,   only: vmr_type

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: finish_ddt_init
  PUBLIC :: finish_ddt_run
  PUBLIC :: finish_ddt_finalize
  
CONTAINS

!> \section arg_table_finish_ddt_run  Argument Table
!! \htmlinclude arg_table_finish_ddt_run.html
!!
  SUBROUTINE finish_ddt_run(nbox, nlev, vmr, O3, HNO3,        &
       errmsg, errflg)
!----------------------------------------------------------------
   IMPLICIT NONE
!----------------------------------------------------------------

   integer,            intent(in)    :: nbox, nlev
   type(vmr_type),     intent(in)    :: vmr
   REAL(kind_phys),    intent(out)   :: O3(nbox,nlev)
   REAL(kind_phys),    intent(out)   :: HNO3(nbox,nlev)
   character(len=512), intent(out)   :: errmsg
   integer,            intent(out)   :: errflg
!----------------------------------------------------------------

    errmsg = ''
    errflg = 0
 

    ! NOTE -- This is prototyping one approach to passing a large number of
    ! chemical VMR values and is the predecssor for adding in methods and maybe
    ! nesting DDTs (especially for aerosols)
    O3(:,:)   = vmr%vmr_array(:,:,1)
    HNO3(:,:) = vmr%vmr_array(:,:,2)

    write(6,*) 'in finish O3(1,1)=',O3(1,1),O3(1,1)*2._kind_phys

  END SUBROUTINE finish_ddt_run

!> \section arg_table_finish_ddt_init  Argument Table
!! \htmlinclude arg_table_finish_ddt_init.html
!!
  subroutine finish_ddt_init (errmsg, errflg)

    character(len=512),      intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

    ! This routine currently does nothing

    errmsg = ''
    errflg = 0

  end subroutine finish_ddt_init

!> \section arg_table_finish_ddt_finalize  Argument Table
!! \htmlinclude arg_table_finish_ddt_finalize.html
!!
  subroutine finish_ddt_finalize (errmsg, errflg)

    character(len=512),      intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

    ! This routine currently does nothing

    errmsg = ''
    errflg = 0

  end subroutine finish_ddt_finalize

END MODULE finish_ddt
