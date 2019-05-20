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
  SUBROUTINE finish_ddt_run(nbox, nlev, vmr, temp_constant, O3, HNO3,        &
       errmsg, errflg)
!----------------------------------------------------------------
   IMPLICIT NONE
!----------------------------------------------------------------

   integer,            intent(in)    :: nbox, nlev
   type(vmr_type),     intent(in)    :: vmr
   real(kind_phys),    intent(in)    :: temp_constant(nbox,nlev)
   REAL(kind_phys),    intent(out)   :: O3(nbox,nlev)
   REAL(kind_phys),    intent(out)   :: HNO3(nbox,nlev)
   character(len=512), intent(out)   :: errmsg
   integer,            intent(out)   :: errflg
!----------------------------------------------------------------

   REAL(kind_phys)   :: sum_result(nbox,nlev)
    
    errmsg = ''
    errflg = 0
 

    ! NOTE -- This is prototyping one approach to passing a large number of
    ! chemical VMR values and is the predecssor for adding in methods and maybe
    ! nesting DDTs (especially for aerosols)
    O3(:,:)   = vmr%vmr_array(:,:,1)
    HNO3(:,:) = vmr%vmr_array(:,:,2)

    write(6,*) 'in finish O3(1,1)=',O3(1,1),O3(1,1)*2._kind_phys
    write(6,*) 'in finish HNO3(1,1)=',HNO3(1,1),HNO3(1,1)*2._kind_phys
    write(6,*) 'in finish temp_constant(1,1) =', temp_constant(1,1)

    call vmr%sum(temp_constant, sum_result)
    write(6,*) 'sum_result(1,1)=',sum_result(1,1)

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
