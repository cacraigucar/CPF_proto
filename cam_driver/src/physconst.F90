MODULE physconst

   implicit none
   public

   !> \section arg_table_physconst  Argument Table
   !! \htmlinclude arg_table_physconst.html
   !!

   integer, parameter :: kind_phys = selected_real_kind(12) ! 8 byte real

   !----------------------------------------------------------------------------
   ! physical constants (all data public)
   !----------------------------------------------------------------------------

   real(kind_phys),parameter :: CPDAIR    = 1.00464e3_kind_phys              ! specific heat of dry air   ~ J/kg/K
   real(kind_phys),parameter :: RDAIR     = RGAS/MWDAIR               ! Dry air gas constant     ~ J/K/kg
   real(kind_phys),parameter :: LATVAP    = 2.501e6_kind_phys                ! latent heat of evaporation ~ J/kg
   real(kind_phys),parameter :: PSTD      = 101325.0_kind_phys               ! standard pressure ~ pascals
   real(kind_phys),parameter :: RHOFW     = 1.000e3_kind_phys                ! density of fresh water     ~ kg/m^3
   real(kind_phys),parameter :: GRAVIT    = 9.80616_kind_phys                ! acceleration of gravity ~ m/s^2
   real(kind_phys),parameter :: ZVIR      = RWV/RDAIR-1.0_kind_phys          ! RWV/RDAIR - 1.0  

   ! Values without metadata (standard names)
   real(kind_phys),parameter :: PI        = 3.14159265358979323846_kind_phys ! pi
   real(kind_phys),parameter :: AVOGADRO  = 6.02214e26_kind_phys             ! Avogadro's number ~ molecules/kmole
   real(kind_phys),parameter :: BOLTZMANN = 1.38065e-23_kind_phys            ! Boltzmann's constant ~ J/K/molecule
   real(kind_phys),parameter :: RGAS      = AVOGADRO*BOLTZMANN        ! Universal gas constant ~ J/K/kmole
   real(kind_phys),parameter :: MWDAIR    = 28.966_kind_phys                 ! molecular weight dry air ~ kg/kmole
   real(kind_phys),parameter :: MWWV      = 18.016_kind_phys                 ! molecular weight water vapor 
   real(kind_phys),parameter :: RWV       = RGAS/MWWV                 ! water vapor gas constant  ~ J/K/kg

!-----------------------------------------------------------------------------

END MODULE physconst
