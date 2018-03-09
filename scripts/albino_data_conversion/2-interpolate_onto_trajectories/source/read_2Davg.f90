!=======================================================================
!
!     read_2Davg
!
!=======================================================================

      subroutine read_2Davg(filename,dens,temp,ye,vel)

      use parameter_module
      use units_module

      implicit none

      character(len=100)               , intent(in)  :: filename
      real(r_kind4), dimension(nr,nz)  , intent(out) :: dens
      real(r_kind4), dimension(nr,nz)  , intent(out) :: temp
      real(r_kind4), dimension(nr,nz)  , intent(out) :: ye
      real(r_kind4), dimension(3,nr,nz), intent(out) :: vel

!-----------------------------------------------------------------------
!
!    Subroutine to read the average data produced in the nu-driven wind
!    from the aftermath of a binary NS merger simulation, and store them
!    in cylindrical grids provided in the call
!
!    The 2D-cylindrical grids are assumed
!    I)   to be smaller than the data grids
!    II)  to have an even number of vertical zones, which cover both the N
!         and S emispheres
!    III) with the same resolution as the data grid (provided in the subroutine)
!    
!
!    Input:
!
!    nr       ... number of radial zones in the grid to be filled
!    nz       ... number of vertical zones in the grid to be filled
!    filename ... name of the file with the average data
!
!    Output:
!
!    dens     ... density grid
!    temp     ... temperature grid
!    ye       ... electron fraction grid
!    vel      ... grid with cylindrical components of the velocity
!             ... vel(1) radial
!             ... vel(2) phi
!             ... vel(3) z
!
!-----------------------------------------------------------------------

      character(len=8) :: date
      integer          :: num       
      integer          :: ig_cyl
      integer          :: kg_cyl
      integer          :: i,j
      integer          :: ii,jj
      integer          :: ne_chk
      integer          :: nv_cyl
      integer          :: ns_cyl
      real(r_kind8)    :: dx_in
      real(r_kind8)    :: time


      real(r_kind4),dimension(:)        ,allocatable :: rho_cyl !rho coordinates
      real(r_kind4),dimension(:)        ,allocatable :: z_cyl   !z coordinates
      real(r_kind8),dimension(:)        ,allocatable :: e,de    !energy bins and intervals
      real(r_kind4),dimension(  :,:,:)  ,allocatable :: s_cyl   !scalar data array
      real(r_kind4),dimension(:,:,:,:)  ,allocatable :: v_cyl   !vector data array

!.....check on nz
      if(mod(nz,2).ne.0) then
        write(6,*)'Error: nz must be even!'
        write(6,*)'change nz in the main file'
        stop
      end if  

!.....open file
      write(6,*)trim(filename)
      open(1,file=trim(filename),status='unknown',form='unformatted')

!.....read date, num & time
      read(1) date
      read(1) num
      read(1) time ![s]

      write(6,*)'time from SPH mapping [ms]',1.e+3*time
      write(6,*)'time in wind simulation [ms]',1.e+3*(time-t0)

!.....read spatial resolution
      read(1) dx_in
      if (dx_in.ne.dx) then
        write(6,*)'Problem with the spatial resolution dx'
        stop
      end if
      write(6,*)'Grid spatial resolution: ',dx_in/1.e+5,' [km]'

!.....conversion factor
!      fd = (c/dx_in)**2/G

!.....read data size in rho direction and tho coordinates
      read(1) ig_cyl
      allocate(rho_cyl(ig_cyl))
      read(1) rho_cyl

!.....check on radial dimension
      if(ig_cyl.lt.nr) then
        write(6,*)'Error: requested array too large in rho direction'
        write(6,*)'Reduce nr in the main file'
        stop
      end if

!.....read data size in z direction and z coordinates
      read(1) kg_cyl
      allocate(z_cyl(kg_cyl))
      read(1) z_cyl
  
!.....check on vertical dimension
      if(kg_cyl/2.lt.nz) then
        write(6,*)'Error: requested array too large in z direction'
        write(6,*)'Reduce nz in the main file'
        stop
      end if

      write(6,*)'dimensions',ig_cyl,kg_cyl

!.....read number of energy bins
      read(1) ne_chk
      if (ne.ne.ne_chk) then
        write(6,*)'ne and ne_chk do not confirm'
        stop
      end if

!.....read energy bins and energy intervals
      allocate(e(ne_chk),de(ne_chk))
      read(1) e
      read(1) de

!.....read scalar data
      read(1) ns_cyl
      allocate(s_cyl(  ns_cyl,ig_cyl,kg_cyl))
      read(1) s_cyl(:,:,:)

!.....read vector data
      read(1) nv_cyl
      allocate(v_cyl(3,nv_cyl,ig_cyl,kg_cyl))
      read(1) v_cyl(:,:,:,:)

      close(1)

!.....place the results in the output grid..............................
      do i=1,nr
        do j=1,nz
          ii = i
!          jj = kg_cyl/2 - nz/2 + j    ! all z
          jj = kg_cyl/2 + j            ! only z > 0
          dens(i,j)  = fd * s_cyl(1,ii,jj)            ![g/cm^3]
          temp(i,j)  = s_cyl(2,ii,jj)                 ![MeV]
          ye(i,j)    = s_cyl(3,ii,jj)                 ![-]
          vel(1,i,j) = v_cyl(1,1,ii,jj)*units%c       ![cm/s]
          vel(2,i,j) = v_cyl(2,1,ii,jj)*units%c       ![cm/s]
          vel(3,i,j) = v_cyl(3,1,ii,jj)*units%c       ![cm/s]
        end do
      end do
56    format(2i5,6es14.4)

!      do i=1,nr
!        write(56,*)
!        do j=1,nz
!          write(56,66)i,j,dens(i,j)
!        end do
!      end do
!66    format(2i5,6es14.4)

!.....deallocate scalar, vector data & cylindrical coordinates arrays...
      deallocate(e,de)
      deallocate(s_cyl)
      deallocate(v_cyl)
      deallocate(rho_cyl)
      deallocate(z_cyl)

      write(6,*)'Done'

      end subroutine read_2Davg

!=======================================================================
