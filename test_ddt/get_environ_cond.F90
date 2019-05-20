!Hello demonstration parameterization
!

MODULE get_environ_cond

  USE ccpp_kinds, ONLY: kind_phys

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: get_environ_cond_init
  PUBLIC :: get_environ_cond_run
  PUBLIC :: get_environ_cond_finalize

CONTAINS

!> \section arg_table_get_environ_cond_run  Argument Table
!! \htmlinclude arg_table_get_environ_cond_run.html
!!
  subroutine get_environ_cond_run(errmsg, errflg)

    ! This routine currently does nothing -- should update values

    character(len=512),      intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

    errmsg = ''
    errflg = 0

  END SUBROUTINE get_environ_cond_run

!> \section arg_table_get_environ_cond_init  Argument Table
!! \htmlinclude arg_table_get_environ_cond_init.html
!!
  subroutine get_environ_cond_init (nbox, nlev, O3, HNO3, temp_constant, errmsg, errflg)

   integer,            intent(in)    :: nbox, nlev
   real(kind_phys),    intent(out)   :: O3(:,:)
   real(kind_phys),    intent(out)   :: HNO3(:,:)
   real(kind_phys),    intent(out)   :: temp_constant(:,:)
   character(len=512), intent(out)   :: errmsg
   integer,            intent(out)   :: errflg
!----------------------------------------------------------------

   integer :: i, j

    errmsg = ''
    errflg = 0
 
    ! This may be replaced with MusicBox json environmental conditions reader???

    do j=1,nlev
      do i=1,nbox
         O3(i,j)   = i*1.e-6
         HNO3(i,j) = j*2000*1.e-9
      end do
    end do

    temp_constant = 250._kind_phys

    write(6,*) ' O3(1,1)=',O3(1,1)

  end subroutine get_environ_cond_init

!> \section arg_table_get_environ_cond_finalize  Argument Table
!! \htmlinclude arg_table_get_environ_cond_finalize.html
!!
  subroutine get_environ_cond_finalize (errmsg, errflg)

    character(len=512),      intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

    ! This routine currently does nothing

    errmsg = ''
    errflg = 0

  end subroutine get_environ_cond_finalize

END MODULE get_environ_cond
