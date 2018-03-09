!=======================================================================
!
!     read_luminosity: reading cylindrically averaged neutrino luminosity
!
!=======================================================================

      subroutine read_luminosity(filename,time,rate,direction,map,e,de)

      use parameter_module

      implicit none

      character(len=100)                                     ,intent(in)  :: filename
      real(r_kind8)                                          ,intent(out) :: time
      real(r_kind4),dimension(ne_h,nut_h,cyl_h_nx,cyl_h_nz)  ,intent(out) :: rate
      real(r_kind4),dimension(nut_h,cyl_h_nx,cyl_h_nz,2)     ,intent(out) :: direction
      logical      ,dimension(ne_h,nut_h,cyl_h_nx,cyl_h_nz,4),intent(out) :: map
      real(r_kind ),dimension(ne_h)                          ,intent(out) :: e
      real(r_kind ),dimension(ne_h)                          ,intent(out) :: de


!-----------------------------------------------------------------------

      character(6)   :: suffix
      character(8)   :: date
      integer        :: m

      integer :: num_chk
      integer :: cyl_h_nx_chk
      integer :: cyl_h_nz_chk
      integer :: ne_h_chk
      integer :: nut_h_chk

      real(r_kind8),dimension(:)        ,allocatable :: e_in,de_in    !energy bins and intervals

!.....read 2D cylindrically averaged data...............................
      open(1,file=trim(filename),status='old',form='unformatted')

      !read date, num & time
      read(1) date
      read(1) num_chk
      read(1) time    ![s]
      !read data size in rho direction
      read(1) cyl_h_nx_chk
      if (cyl_h_nx_chk.ne.cyl_h_nx) then 
        write(6,*)'Wrong cyl_h_nx!'
        write(6,*)'given:',cyl_h_nx
        write(6,*)'expected:',cyl_h_nx_chk
      end if
      !read data size in z direction
      read(1) cyl_h_nz_chk
      if (cyl_h_nz_chk.ne.cyl_h_nz) then 
        write(6,*)'Wrong cyl_h_nz!'
        write(6,*)'given:',cyl_h_nz
        write(6,*)'expected:',cyl_h_nz_chk
      end if
      !read number of energy bins
      read(1) ne_h_chk
      if (ne_h_chk.ne.ne_h) then 
        write(6,*)'Wrong ne_h!'
        write(6,*)'given:',ne_h
        write(6,*)'expected:',ne_h_chk
      end if
      !read energy bins and energy intervals
      allocate(e_in(ne_h),de_in(ne_h))
      read(1) e_in(1:ne_h)
      read(1) de_in(1:ne_h)
      e = e_in
      de = de_in
      !read number of flavour
      read(1) nut_h_chk
      if (nut_h_chk.ne.nut_h) then 
        write(6,*)'Wrong nut_h!'
        write(6,*)'given:',nut_h
        write(6,*)'expected:',nut_h_chk
      end if
      !read rate
      read(1) rate
      !read direction
      read(1) direction
      !read map
      read(1) map
      close(1)
 
!      do i=1,cyl_h_nx
!        do j=1,cyl_h_nz
!          write(61,*)direction(1,i,j,1),direction(1,i,j,2)
!          write(62,*)i,j,map(6,1,i,j,1)
!          write(63,*)i,j,map(6,1,i,j,2)
!          write(64,*)i,j,map(6,1,i,j,3)
!        end do
!      end do

      deallocate(e_in,de_in)

      end subroutine read_luminosity

!=======================================================================
