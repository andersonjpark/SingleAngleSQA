!***********************************************************************
!
!     module: tool_module
!
!***********************************************************************

      module tool_module

      use parameter_module

      implicit none

      contains

!=======================================================================
!
!     subroutine: bisec_search
!
!=======================================================================

      subroutine bisec_search(n,array,val,ind,delta)

      implicit none

      integer                    ,intent(in) :: n
      real(r_kind4), dimension(n),intent(in) :: array
      real(r_kind)               ,intent(in) :: val
      integer                    ,intent(out) :: ind
      integer                    ,intent(out) :: delta

      integer :: i,il,iu
      real(r_kind), dimension(2) :: eps

      il = 0
      iu = n + 1
      do while(iu-il.gt.1)
        i = (iu+il)/2
        if(array(i).gt.val) then
          iu = i
        else
          il = i
        end if      
      end do

      if(il.eq.0) then
        ind = 1
        delta = -1
      else if(iu.eq.n+1) then
        ind = n
        delta = +1
      else
        eps(1) =  abs(array(il)-val)
        eps(2) =  abs(array(iu)-val)
        if(eps(1).lt.eps(2)) then
          ind = il
          delta = -1
        else
          ind = iu
          delta = +1
        end if
      end if

      end subroutine bisec_search

!=======================================================================

!=======================================================================
!
!     subroutine cartesian2cylindrical
!
!=======================================================================

      subroutine cartesian2cylindrical(x,rc)

      use units_module

      implicit none

      real(r_kind), dimension(3), intent(in)  :: x
      real(r_kind), dimension(3), intent(out) :: rc

      real(r_kind), parameter :: halfpi      = 0.5*units%pi
      real(r_kind), parameter :: pi          = units%pi
      real(r_kind), parameter :: threehalfpi = 1.5*units%pi
      real(r_kind), parameter :: twopi       = 2.0*units%pi

!.....radial coordinate.................................................
      rc(1) = sqrt(x(1)*x(1) + x(2)*x(2))

!.....phi coordinate....................................................
      if (x(2).gt.0.) then
        if (x(1).gt.0.) then
          rc(2) = atan(x(2)/x(1))
        else if(x(1).eq.0.) then
          rc(2) = halfpi
        else if (x(1).lt.0.) then
          rc(2) = atan(x(2)/x(1)) + pi
        endif
      elseif (x(2).eq.0.) then
        if(x(1).ge.0.) then
          rc(2) = 0.
        else if (x(1).lt.0.) then
          rc(2) = pi
        endif
      else if (x(2).lt.0.) then
        if (x(1).gt.0.) then
          rc(2) = twopi + atan(x(2)/x(1)) 
        else if(x(1).eq.0.) then
          rc(2) = threehalfpi
        else if (x(1).lt.0.) then
          rc(2) = pi + atan(x(2)/x(1)) 
        endif
      endif

!.....zeta coordinate...................................................
      rc(3) = x(3)

      end subroutine cartesian2cylindrical

!=======================================================================

!=======================================================================
!
!     subroutine cartesian2spherical
!
!=======================================================================

      subroutine cartesian2spherical(x,r)

      use units_module

      implicit none

      real(r_kind), dimension(3), intent(in)  :: x
      real(r_kind), dimension(3), intent(out) :: r

      real(r_kind), parameter :: halfpi      = 0.5*units%pi
      real(r_kind), parameter :: pi          = units%pi
      real(r_kind), parameter :: threehalfpi = 1.5*units%pi
      real(r_kind), parameter :: twopi       = 2.0*units%pi

!.....radial coordinate.................................................
      r(1) = sqrt(x(1)*x(1) + x(2)*x(2) + x(3)*x(3))

!.....theta coordinate..................................................
      if(r(1).gt.0.) then
        r(2) = ACOS(x(3)/r(1))
      else
        r(2) = 0.
      endif

!.....phi coordinate....................................................
      if (x(2).gt.0.) then
        if (x(1).gt.0.) then
          r(3) = atan(x(2)/x(1))
        else if(x(1).eq.0.) then
          r(3) = halfpi
        else if (x(1).lt.0.) then
          r(3) = atan(x(2)/x(1)) + pi
        endif
      elseif (x(2).eq.0.) then
        if(x(1).ge.0.) then
          r(3) = 0.
        else if (x(1).lt.0.) then
          r(3) = pi
        endif
      else if (x(2).lt.0.) then
        if (x(1).gt.0.) then
          r(3) = twopi + atan(x(2)/x(1)) 
        else if(x(1).eq.0.) then
          r(3) = threehalfpi
        else if (x(1).lt.0.) then
          r(3) = pi + atan(x(2)/x(1)) 
        endif
      endif

      end subroutine cartesian2spherical

!=======================================================================

      end module tool_module

!***********************************************************************
