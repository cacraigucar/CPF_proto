MODULE vmr_ddt_mod

  USE ccpp_kinds, ONLY: kind_phys

  PUBLIC :: vmr_type
  PUBLIC :: vmr_ddt_init
  
!> \section arg_table_vmr_type  Argument Table
!! \htmlinclude arg_table_vmr_type.html
!!
  type vmr_type
   integer :: nvmr
   real(kind_phys), allocatable :: vmr_array(:,:,:)
  end type
    
END MODULE vmr_ddt_mod
