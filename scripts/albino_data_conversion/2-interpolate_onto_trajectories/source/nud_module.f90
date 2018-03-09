!=======================================================================
!
!     module :: nud_module
!
!=======================================================================

      module nud_module

      use parameter_module
      use tool_module
      use units_module

      implicit none

      type spherical_dist
        integer(2), dimension(2) :: point
        real(r_kind4) :: dist
      end type

      integer :: n_dist_red

      real(r_kind) :: domega_res
      real(r_kind),  dimension(ne_h)           :: e,de
      real(r_kind4), dimension(nth)            :: theta
      real(r_kind4), dimension(nphi)           :: phi
      real(r_kind4), dimension(nphip)          :: phip
      real(r_kind4), dimension(nphip)          :: cosphip
      real(r_kind4), dimension(nphip)          :: sinphip
      real(r_kind4), dimension(3,nth,nphi)     :: n_dir
      real(r_kind4), dimension(nth,nphi,nphip) :: nrho

      real(r_kind4), dimension(cyl_h_nx) :: rp
      real(r_kind4), dimension(-cyl_h_nz+1:cyl_h_nz) :: zp
      real(r_kind4), dimension(cyl_h_nx,cyl_h_nz) :: cyl_h_dV
      real(r_kind4), dimension(cyl_h_nx,cyl_h_nz) :: cyl_h_dV_red

      type(spherical_dist),dimension(:,:,:), allocatable :: sph_dist
      type(spherical_dist),dimension(:,:,:), allocatable :: sph_dist_red

!.....array with neutrino data..........................................
      real(r_kind4), dimension(ne_h,nut_h,cyl_h_nx,cyl_h_nz) :: ser2dc
      real(r_kind4), dimension(2,nut_h,cyl_h_nx,cyl_h_nz)    :: n_tau
      logical, dimension(ne_h,nut_h,cyl_h_nx,cyl_h_nz,4)     :: map_read
      logical, dimension(4,ne_h,nut_h,cyl_h_nx,cyl_h_nz)     :: map

!.....array with disc data..............................................
      real(r_kind8) :: time
      real(r_kind4), dimension(nr,nz)   :: dens
      real(r_kind4), dimension(nr,nz)   :: temp
      real(r_kind4), dimension(nr,nz)   :: ye
      real(r_kind4), dimension(3,nr,nz) :: vel

      contains

!=======================================================================
!
!     subroutine: direction_initialize
!
!=======================================================================

      subroutine direction_initialize

      implicit none

      integer :: i,j,k,ii,jj,kk
      real(r_kind) :: dcosth
      real(r_kind) :: dphi
      real(r_kind) :: dphip
      real(r_kind) :: tmp
      real(r_kind) :: domega_max
      real(r_kind) :: inner_dV
      real(r_kind), dimension(:), allocatable :: sinth
      real(r_kind), dimension(:), allocatable :: costh
      real(r_kind), dimension(:), allocatable :: cosphi,sinphi
      type(spherical_dist) :: sph_dist_tmp

      allocate(costh(nth),sinth(nth))
      allocate(cosphi(nphi),sinphi(nphi))
      allocate(sph_dist(nth,nphi,nth*nphi))

!.....initialize possible theta angles and their cosines................
      dcosth = 2./dble(nth)             ! [-]
      costh(1) = max(1. - dcosth/2.,-1.)          ! [-]
      sinth(1) = sqrt(1.-costh(1)**2)    ! [-]
      theta(1) = acos(costh(1))          ! [rad]
      do i=2,nth
        costh(i) = max(costh(i-1) - dcosth,-1.)   ! [-]
        sinth(i) = sqrt(1.-costh(i)**2)  ! [-]
        theta(i) = acos(costh(i))        ! [rad]
      end do

!.....initialize possible phi angles and their sines and cosines........
      dphi = (2.*units%pi)/dble(nphi)   ! [rad]
      phi(1) = dphi/2.                   ! [rad]
      cosphi(1) = cos(phi(1))            ! [-]
      sinphi(1) = sin(phi(1))            ! [-]
      do j=2,nphi
        phi(j) = phi(j-1) + dphi         ! [rad]
        cosphi(j) = cos(phi(j))          ! [-]
        sinphi(j) = sin(phi(j))          ! [-]
      end do

!.....initialize the solid angle resolution.............................
      domega_res = dcosth*dphi           ! [sterad]

!.....initialize phip angle of the integration ring.....................
      dphip = (2.*units%pi)/dble(nphip) ! [rad]
      phip(1) = dphip/2.                 ! [rad]
      cosphip(1) = cos(phip(1))          ! [-]
      sinphip(1) = sin(phip(1))          ! [-]
      do j=2,nphip
        phip(j) = phip(j-1) + dphip      ! [rad]
        cosphip(j) = cos(phip(j))        ! [-]
        sinphip(j) = sin(phip(j))        ! [-]
      end do

!.....calculate the direction unitary vectors...........................
      do i=1,nth
        do j=1,nphi
          n_dir(1,i,j) = sinth(i) * cosphi(j)  ! [-]
          n_dir(2,i,j) = sinth(i) * sinphi(j)  ! [-]
          n_dir(3,i,j) = costh(i)              ! [-]
        end do
      end do

      write(6,*)'Angles and unitary directions have been initialized'

!.....caluclate spherical distances for the direction vectors...........
      do i=1,nth
        do j=1,nphi
          k = 1
          do ii=1,nth
            do jj=1,nphi
              tmp = min(max(                                            &
     &            sum(dble(n_dir(:,i,j)*n_dir(:,ii,jj))),-1.),1.)   ![-]
              sph_dist(i,j,k)%point = (/ ii , jj /)    ![-]
              sph_dist(i,j,k)%dist = acos(tmp)         ![-]
              k = k + 1
            end do
          end do
        end do
      end do

!.....order the distances, from the smallest to the largest.............
      do i=1,nth
        do j=1,nphi
           
          do k = 1,nth*nphi-1
            do kk = 1,nth*nphi-k
              if(sph_dist(i,j,kk)%dist.gt.sph_dist(i,j,kk+1)%dist)then
                sph_dist_tmp = sph_dist(i,j,kk)
                sph_dist(i,j,kk) = sph_dist(i,j,kk+1)
                sph_dist(i,j,kk+1) = sph_dist_tmp
              else
                continue
              end if
            end do
          end do
        
        end do
      end do

!.....reduce the distances vector to the first neighbours...............
      domega_max = 0.9
      n_dist_red = floor(domega_max/domega_res) + 2
      allocate(sph_dist_red(nth,nphi,n_dist_red))

      do i=1,nth
        do j=1,nphi
          do k=1,n_dist_red
            sph_dist_red(i,j,k)%dist = sph_dist(i,j,k)%dist
            sph_dist_red(i,j,k)%point(1:2) = sph_dist(i,j,k)%point(1:2)
          end do
        end do
      end do
     
      write(6,*)'Spherical distances have been calculated and ordered'

!.....radial cylindrical coordinate: high resolution part...............
      rp(1) = 0.5
      do i=2,cyl_h_nx
        rp(i) = rp(i-1) + 1.   ![c.u.]
      end do

!.....vertical cylindrical coordinate...................................
      zp(0) = - 0.5
      zp(1) = 0.5
      do i=2,cyl_h_nz
        zp(i) = zp(i-1) + 1.   ![c.u.]
        zp(-i+1) = - zp(i)
      end do

      write(6,*)'heating cylindrical (centred) coordinates initialized'

!.....initialize the heating cylindrical volumes........................
      cyl_h_dV(1,:) = units%pi*(rp(1) + 0.5 )**2
      inner_dV = cyl_h_dV(1,1)
      do i=2,cyl_h_nx
        cyl_h_dV(i,:) = units%pi*(rp(i) + 0.5)**2
        cyl_h_dV(i,:) = cyl_h_dV(i,:) - inner_dV
        inner_dV = inner_dV + cyl_h_dV(i,1)
      end do

!.....initialize the reduced heating cylindrical volumes................
      cyl_h_dV_red = cyl_h_dV/dble(nphip)

      write(6,*)'heating cylindrical volumes initialized'

!.....calculate the radial (cylindrical) components of the unitary ..... to be checked!!!
!     vectors, for each possible phi' value............................. 
      do k=1,nphip
        do j=1,nphi
          do i=1,nth
            nrho(i,j,k) = n_dir(1,i,j)*cosphip(k) +                     &
     &               n_dir(2,i,j)*sinphip(k)           ![-]
          end do
        end do
      end do
      write(6,*)'radial components of the directions initialized'

      deallocate(costh,sinth)
      deallocate(sinphi,cosphi)
      deallocate(sph_dist)

      end subroutine direction_initialize

!=======================================================================

!=======================================================================
!
!     subroutine: pref_dir
!
!=======================================================================

      subroutine pref_dir(r,ith_c,iphi_c)

      implicit none

      real(r_kind), dimension(3), intent(in)  :: r
      integer                   , intent(out) :: ith_c
      integer                   , intent(out) :: iphi_c

!-----------------------------------------------------------------------
!
!     Input:
!     r      ... spherical coordinates of the point A [c.u.]
!
!     Output:
!     ith_c  ... index for the central theta angle bin [-]
!     iphi_c ... index for the central phi angle bin [-]
!
!-----------------------------------------------------------------------

      integer :: delta

!.....search theta angle closest bin....................................
      call bisec_search(nth,theta,r(2),ith_c,delta)

!.....search phi angle closest bin......................................
      call bisec_search(nphi,phi,r(3),iphi_c,delta)

      end subroutine pref_dir

!=======================================================================

!=======================================================================
!
!     subroutine: relative_v
!
!=======================================================================

      subroutine relative_v(p,x,ic,kc,iphip,r_ca)

      implicit none

      logical                   ,intent(in)  :: p
      real(r_kind), dimension(3),intent(in)  :: x
      integer                   ,intent(in)  :: ic
      integer                   ,intent(in)  :: kc
      integer                   ,intent(in)  :: iphip
      real(r_kind), dimension(3),intent(out) :: r_ca

      r_ca(1) = x(1)-rp(ic)*cosphip(iphip)   ![c.u.]
      r_ca(2) = x(2)-rp(ic)*sinphip(iphip)   ![c.u.]
      if (p) then
        r_ca(3) = x(3)-zp(kc)                ![c.u.]
      else
        r_ca(3) = x(3)-zp(-kc+1)             ![c.u.]
      end if

      end subroutine relative_v

!=======================================================================

!=======================================================================
!
!     subroutine: pot_calc
!
!=======================================================================

      subroutine pot_calc(x,qx,ic,kc,iphip,sr,nu_mp_loc,nu_mp,          &
     &                    tau_fin,tau_source,dn,dl)
      implicit none

      real(r_kind), dimension(3)               ,intent(in)    :: x
      real(r_kind), dimension(3)               ,intent(in)    :: qx
      integer                                  ,intent(in)    :: ic
      integer                                  ,intent(in)    :: kc
      integer                                  ,intent(in)    :: iphip
      real(r_kind4), dimension(ne_h,nut_h)     ,intent(in)    :: sr
      logical,dimension(4,ne_h,nut_h)          ,intent(in)    :: nu_mp_loc
      logical,dimension(4,ne_h,nut_h)          ,intent(in)    :: nu_mp
      real(r_kind), dimension(ne_h,nut_h)      ,intent(in)    :: tau_fin
      real(r_kind), dimension(ne_h,nut_h)      ,intent(in)    :: tau_source
      real(r_kind), dimension(ne_h,nut_h)      ,intent(inout) :: dn
      real(r_kind), dimension(ne_h,nut_h)      ,intent(inout) :: dl

      logical :: p
      integer :: ith,iphi,ip
      integer :: ith_c,iphi_c
      integer :: ie,it,kk
      integer :: num_omega
      real(r_kind) :: domega
      real(r_kind) :: mu
      real(r_kind) :: dist_fac
      real(r_kind) :: tmp,tmp1,tmp2
      real(r_kind) :: mu_coeff
      real(r_kind) :: qx_mod
      real(r_kind) :: costh_qqp
      real(r_kind), parameter :: pref1 = dx / units%c * fd  ![cm/lenght[c.u.] * s/cm * g/cm^3/density[c.u.] ]
      real(r_kind), dimension(3) :: r_ca,r_sp,r_cy
      real(r_kind), dimension(ne_h,nut_h) :: dmp_tau

      real(r_kind), parameter :: dist_norm2 = 1./units%pi        !gamma = 1.
!      real, parameter :: dist_norm2 = 3./4./units%pi     !gamma = 0.5
      real(r_kind), parameter :: dist_norm3 = 1./2./units%pi     !gamma = 0.
   
      character(len=100) :: fileNumber
      COMMON /YZHU14/ fileNumber
  
!      open(66,file='./output/xcosf'//trim(fileNumber)//'.txt',status='unknown')
!      open(88,file='./output/kkitie'//trim(fileNumber)//'.txt',status='unknown')
      do ip=1,2

        if (ip.eq.1) then
          p=.true.               ! z > 0
          mu_coeff = +1. 
        else if (ip.eq.2) then
          p=.false.              ! z < 0
          mu_coeff = -1. 
        end if      

!.....compute dump factors based on the optical depth...................
        do ie=1,ne_h 
          do it=1,nut_h
            tmp1 = min(tau_source(ie,it),2./3.)
            tmp2 = min(tau_fin(ie,it),tmp1)
            dmp_tau(ie,it) = exp( tmp2 - tmp1 )
          end do
        end do 

!.....relative position vector in Cartesian coordinates.................
        call relative_v(p,x,ic,kc,iphip,r_ca)

!.....relative position vector in spherical coordinates.................
        call cartesian2spherical(r_ca,r_sp)

!.....relative position vector in cylindrical coordinates...............
        call cartesian2cylindrical(r_ca,r_cy)

        if (r_sp(1).lt.1.) return

        call pref_dir(r_sp,ith_c,iphi_c)

!.....compute cosine of the angle formed with the test neutrino.........
        qx_mod = sqrt(sum(qx*qx))
        costh_qqp = sum(r_ca*qx)/r_sp(1)/qx_mod
!        write(66,*)costh_qqp
        domega = 10.**(-1.9765*log10(r_sp(1))-0.069649)      ![sterad]
        num_omega = floor(domega/domega_res) + 1             ![-]
!        if (iphip.eq.1.and.ic.eq.10.and.kc.eq.10) then
!          write(*,*)ip,r_sp(1)!
!        end if
        tmp = (r_sp(1)*r_sp(1))*dble(num_omega)
        tmp = pref1 * cyl_h_dV_red(ic,kc) / tmp

!.....sum over the relevant angles......................................
        do kk=1,num_omega

          ith  = sph_dist_red(ith_c,iphi_c,kk)%point(1)     ![-]
          iphi = sph_dist_red(ith_c,iphi_c,kk)%point(2)     ![-]

          do it=1,nut_h

              mu = (n_tau(1,it,ic,kc)*nrho(ith,iphi,iphip) +            &
     &              mu_coeff*n_tau(2,it,ic,kc)*n_dir(3,ith,iphi))      ![-]

            if (mu.le.0.) cycle

            do ie=1,ne_h
!               write(88,*)ic,kc,ip,kk,ie,it,&
!     &           costh_qqp
             ! if (iphip.eq.1.AND.ic.eq.1.AND.kc.eq.1) then
!             if (kk.eq.2.or.ip.eq.2) then
!                open(66,file='./output/xcosf111.txt',status='unknown')
!                write(66,*)kk,ie,it,ip,costh_qqp,num_omega,domega
!             end if
!.....don't compute the density when the point where the external    ...
!     point is inside the trapped regime (optically thick +          ...
!.....neutrino surface).................................................
              if (nu_mp_loc(1,ie,it).or.nu_mp_loc(2,ie,it)) then
                continue
              else

                if (dump) then
                  tmp1 = tmp * dmp_tau(ie,it)     ! with dump
                else
                  tmp1 = tmp                      ! without dump
                end if

                if (nu_mp(2,ie,it)) then
!               dist_fac  = 1.
                  dl(ie,it) = dl(ie,it) + tmp1*mu*sr(ie,it)*dist_norm2* &
     &                         (1.-costh_qqp)                          ![particles/cm^3]
                  dn(ie,it) = dn(ie,it) + tmp1*mu*sr(ie,it)*dist_norm2 ![particles/cm^3]
!                  write(66,*)'nu_mp(2)',iphip,ic,kc,kk,ie,it,dl(ie,it), &
!     &          sr(ie,it),tmp1*mu*sr(ie,it)*dist_norm2,costh_qqp

                else if (nu_mp(3,ie,it)) then
!               dist_fac  = 0.
!               dist_norm = 1./(2.*units%pi) !* domega_res
                  dl(ie,it) = dl(ie,it) + tmp1*sr(ie,it)*dist_norm3 *    &
     &                         (1.-costh_qqp)                         ![particles/cm^3]
                  dn(ie,it) = dn(ie,it) + tmp1*sr(ie,it)*dist_norm3   ![particles/cm^3]
!                  write(66,*)'nu_mp(3)',iphip,ic,kc,kk,ie,it,dl(ie,it),  &
!     &         sr(ie,it),tmp1*sr(ie,it)*dist_norm3,costh_qqp
                 else
                  continue
                end if
              end if

            end do
          end do
        end do

      end do

      end subroutine pot_calc

!=======================================================================

!=======================================================================
!
!     subroutine: nu_potential
!
!=======================================================================

      subroutine nu_potential(num,x,qx,d,t,y,d_nu,v_nu)

      implicit none

!      include 'omp_lib.h'

      integer                            , intent(in)  :: num
      real(r_kind), dimension(3)         , intent(in)  :: x
      real(r_kind), dimension(3)         , intent(in)  :: qx
      real(r_kind)                       , intent(out) :: d
      real(r_kind)                       , intent(out) :: t
      real(r_kind)                       , intent(out) :: y
      real(r_kind), dimension(ne_h,nut_h), intent(out) :: d_nu
      real(r_kind), dimension(ne_h,nut_h), intent(out) :: v_nu


      integer :: j
      integer :: ie,it
      integer :: ic,kc
      integer :: ir,iz
      integer :: ic1,kc1,iphi1
      integer :: ic_s,kc_s

      logical, dimension(4,ne_h,nut_h) :: map_loc

      real(r_kind)                              :: rr,rz
      real(r_kind)                              :: lumin
      real(r_kind), dimension(3)                :: rc
      real(r_kind), dimension(ne_h,nut_h,nr,nz) :: tau

      character(len=100) :: filename
      character(len=6)   :: suffix

      logical, save :: dir_init=.true.
      logical, save :: read_file=.true.

!.....initialize coordinates, directions and volumes....................
      if (dir_init) then
        call direction_initialize
        dir_init = .false.
      end if

!.....read-in input files...............................................
      if (read_file) then

!.....generate file suffix..............................................
        write(suffix,88) num
88      format('_',i5)
        do j=1,6
          if(suffix(j:j).eq.' ')suffix(j:j)='0'
        end do

!.....generate file name for the neutrino quantities....................
        filename = '../input_data/data'//suffix//'.lum'
        write(6,*)'I am reading the neutrino input file'
        call read_luminosity(filename,time,ser2dc,n_tau,map_read,e,de)

!.....change map disposition............................................
        do kc1=1,cyl_h_nz
          do ic1=1,cyl_h_nx
            do it=1,nut_h
              do ie=1,ne_h
                map(1:4,ie,it,ic1,kc1) =  map_read(ie,it,ic1,kc1,1:4)
              end do
            end do
          end do
        end do

!.....generate file name for the disc quantities........................
        filename = '../input_data/data'//suffix//'.cyl'
        write(6,*)'I am reading the disc profile input file'
        call read_2Davg(filename,dens,temp,ye,vel)

!.....generate file name for the neutrino optical depth.................
        filename = '../input_data/opdep_2D'//suffix//'.cyl'
        write(6,*)'I am reading the neutrino optical depths'
        call read_2Dtau(filename,tau)

        read_file = .false.

      end if

!.....convert from Cartesian to cylindrical coordinates.................
      call cartesian2cylindrical(x,rc)

!.....find local integers and displacements in the cylindrical grid.....
      ic = floor(1.e+5*rc(1)/dx)+1
      rr = 1.e+5*rc(1)/dx - float(ic-1)
      kc = floor(abs(1.e+5*rc(3)/dx))+1
      rz = 1.e+5*rc(3)/dx - float(kc-1)

!.....calculate the neutrino potentials.................................
!$OMP parallel &
!$OMP default(none) private(kc1,ic1,iphi1,map_loc) &
!$OMP shared(ic,kc,x,ser2dc,map,dnu)

!$OMP do collapse(3) schedule(guided) reduction(+:dnu)

      v_nu = 0.
      d_nu = 0.

      do iphi1=1,nphip
        do kc1=1,cyl_h_nz
          do ic1=1,cyl_h_nx

            if(ic.lt.cyl_h_nx.and.kc.lt.cyl_h_nz) then
              map_loc(:,:,:) = map(:,:,:,ic,kc)
            else
              map_loc(1,:,:) = .false.
              map_loc(2,:,:) = .false.
              map_loc(3,:,:) = .true.
              map_loc(4,:,:) = .false.
            end if

            ic_s = ic1 + nint(abs(n_tau(1,1,ic1,kc1)))
            kc_s = kc1 + nint(abs(n_tau(2,1,ic1,kc1)))
!            ic_s = ic1 + 1
!            kc_s = kc1 + 1

            call pot_calc(x,qx,ic1,kc1,iphi1,ser2dc(:,:,ic1,kc1),       &
     &                map_loc,map(:,:,:,ic1,kc1),                       &
     &                tau(:,:,ic,kc),tau(:,:,ic_s,kc_s),                &
     &                d_nu,v_nu)
            !open(66,file='./output/xcosf1.txt',status='unknown')
            !write(66,*)iphi1,kc1,ic1,d_nu,v_nu!,costh_qqp
          end do
        end do
      end do
!$OMP end do
!$OMP end parallel


!.....interpolate the density (log)..................................... 
      d = (1.-rr)*((1.-rz) * log10(dens(ic  ,kc  ))                     &
     &  +              rz  * log10(dens(ic  ,kc+1)))                    &
     &  +     rr *((1.-rz) * log10(dens(ic+1,kc  ))                     &
     &  +              rz  * log10(dens(ic+1,kc+1)))
      d = 10.**d

!.....interpolate the temperature (log)................................. 
      t = (1.-rr)*((1.-rz) * log10(temp(ic  ,kc  ))                     &
     &  +              rz  * log10(temp(ic  ,kc+1)))                    &
     &  +     rr *((1.-rz) * log10(temp(ic+1,kc  ))                     &
     &  +              rz  * log10(temp(ic+1,kc+1)))
      t = 10.**t

!.....interpolate the electron fraction (lin)...........................
      y = (1.-rr)*((1.-rz) * (ye(ic  ,kc  ))                            &
     &  +              rz  * (ye(ic  ,kc+1)))                           &
     &  +     rr *((1.-rz) * (ye(ic+1,kc  ))                            &
     &  +              rz  * (ye(ic+1,kc+1)))
 
      end subroutine nu_potential

!=======================================================================

      end module nud_module

!=======================================================================
