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

   contains
     procedure :: sum => sum_vmr_array

  end type

contains
   subroutine sum_vmr_array(this, result)
      implicit none
 
      class(vmr_type) :: this
      real(kind_phys) :: result (:,:)  
  
      integer :: ivmr

      result = 0._kind_phys

      do ivmr = 1, this%nvmr
        result(:,:) = this%vmr_array(:,:,ivmr) + result(:,:)
      end do
       
   end subroutine sum_vmr_array
    
END MODULE vmr_ddt_mod
