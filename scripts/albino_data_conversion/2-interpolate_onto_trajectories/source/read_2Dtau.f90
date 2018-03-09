!=======================================================================
!
!     read_2Dtau
!
!=======================================================================

      subroutine read_2Dtau(filename,tau)

      use parameter_module
      use units_module

      implicit none

      character(len=100)                       , intent(in)  :: filename
      real(r_kind), dimension(ne_h,nut_h,nr,nz), intent(out) :: tau

!-----------------------------------------------------------------------
!
!    Subroutine to read the cylindrical optical depth in the nu-driven wind
!    from the aftermath of a binary NS merger simulation, and store them
!    in cylindrical grids provided in the call
!
!    Input:
!
!    filename ... name of the file with the average data
!
!    Output:
!
!    tau      ... optical depth grid
!
!-----------------------------------------------------------------------

      integer       :: ig_cyl
      integer       :: kg_cyl
      integer       :: i,j
      integer       :: ii,jj
      integer       :: ie,it,il
      integer       :: ne_chk
      integer       :: nut_chk
      integer       :: ntau_chk
      integer       :: ns_cyl
      real(r_kind8) :: dx_tau

      real(r_kind4),dimension(:)         , allocatable :: rho_cyl      !rho coordinates
      real(r_kind4),dimension(:)         , allocatable :: z_cyl        !z coordinates
      real(r_kind4),dimension(:,:,:)     , allocatable :: s_cyl        !scalar data array
      real(r_kind4),dimension(:,:,:,:,:) , allocatable :: tau_cyl      !tau data array

      real(r_kind8),dimension(:)         , allocatable :: e_tau         !energy bins
      real(r_kind8),dimension(:)         , allocatable :: de_tau        !energy intervals

!.....open file.........................................................
      open(1,file=trim(filename),status='old',form='unformatted')

!.....read spatial resolution...........................................
      read(1) dx_tau
      write(6,*)'Grid spatial resolution for tau: ',dx_tau/1.e+5,' [km]'

!.....read data size in rho direction and tho coordinates...............
      read(1) ig_cyl
      allocate(rho_cyl(ig_cyl))
      read(1) rho_cyl

      read(1) kg_cyl
      allocate(z_cyl(kg_cyl))
      read(1) z_cyl

!.....read number of energy bins........................................
      read(1) ne_chk

!.....check on on the number of energy bins
      if (ne.ne.ne_chk) then
        write(6,*)'Error: number of energy bins do not conform'
        stop
      end if

      allocate(e_tau(ne_chk),de_tau(ne_chk))

!.....read energy bins and energy intervals
      read(1) e_tau
      read(1) de_tau

!.....read scalar data
      read(1) ns_cyl
      allocate(s_cyl(ns_cyl,ig_cyl,kg_cyl))
      read(1) s_cyl(:,:,:)

      !read number of optical depths
      read(1) ntau_chk

!.....check on on the number of optical depths
      if (ntau.ne.ntau_chk) then
        write(6,*)'Error: number of optical depths do not conform'
        stop
      end if

      !read number of energy bins and neutrino flavors for tau
      read(1) ne_chk
      read(1) nut_chk

!.....check on on the number of optical depths
      if (nut_h.ne.nut_chk) then
        write(6,*)'Error: number of neutrino flavors do not conform'
        stop
      end if

!.....read the optical depths
      allocate(tau_cyl(ig_cyl,kg_cyl,ne_chk,nut_chk,ntau_chk))
      read(1) tau_cyl(:,:,:,:,:)

      close(1)

!.....set the output
      do it=1,nut_h
        do ie=1,ne_h
          do i=1,nr
            do j=1,nz
              if (i.le.ig_cyl.and.j.le.kg_cyl/2) then
                ii = i
                jj = kg_cyl/2 + j
                tau(ie,it,i,j)  = tau_cyl(ii,jj,ie,it,2)      ![-]
              else
                tau(ie,it,i,j)  = 1.e-4                       ![-]
              end if
            end do
          end do
        end do
      end do

!.....deallocate scalar, vector data & cylindrical coordinates arrays...
      deallocate(s_cyl)
      deallocate(tau_cyl)
      deallocate(rho_cyl)
      deallocate(z_cyl)

      write(6,*)'Optical depth read'

      end subroutine read_2Dtau

!=======================================================================
