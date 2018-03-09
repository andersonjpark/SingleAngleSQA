!=======================================================================
!
!     module :: parameter_module
!
!=======================================================================

      module parameter_module

      use units_module

      implicit none

      integer, parameter :: r_kind4  = kind(1.e+0) ! single precision
      integer, parameter :: r_kind8  = kind(1.d+0) ! double precision
      integer, parameter :: r_kind   = r_kind8     ! standard precision

      integer, parameter :: ne       = 12          ! number of energy bins
      integer, parameter :: ne_h     = 8           ! number of energy bins (reduced for some quantities)
      integer, parameter :: nut_h    = 3           ! number of flavours
      integer, parameter :: ntau     = 2           ! number of optical depths

!.....neutrino input grid...............................................
      integer, parameter :: cyl_h_nx = 181         ! radial dimension of the neutrino input grid
      integer, parameter :: cyl_h_nz = 65          ! vertical dimension of the neutrino input grid

!.....disc input grid...................................................
      integer, parameter :: nr = 598               ! radial dimension of the disc input grid
      integer, parameter :: nz = 600               ! vertical dimension of the disc input grid

      integer, parameter :: nth      = 20          ! number of theta bins for the local decomposition
      integer, parameter :: nphi     = 30          ! number of phi bins for the local decomposition
      integer, parameter :: nphip    = 30          ! number of phi bins for the cylindrical grids

      real(r_kind), parameter :: dx = 1.00d+5      ! resolution
      real(r_kind), parameter :: fd = (units%c/dx)**2/units%G    ![g/cm^3]

      real(r_kind), parameter :: t0 =  0.010  ![s] ! initial time after SPH remapping

      logical, parameter :: dump = .true.

      end module parameter_module

!=======================================================================
