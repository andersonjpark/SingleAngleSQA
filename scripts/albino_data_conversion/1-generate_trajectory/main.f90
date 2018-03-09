      program main

      implicit none

      integer, parameter :: r_kind = 8
      integer :: i

!.....number of points..................................................
      integer, parameter :: ntr = 800

!.....step along the trajectory.........................................
      real(r_kind), parameter :: deltax = 0.2  ![km]

!.....step of starting point, for calculation ending at the same
!point,studying the effect of 1-cosTheta
      real(r_kind), parameter :: deltay = 10.  ![km] chaning along y 

!.....Cartesian direction vector (also not normalized).................
!      real(r_kind), dimension(3), parameter :: u = (/ 0.,0.,1./)

!.....Cartesian starting point..........................................
      real(r_kind), dimension(3) :: u0 =(/0.,30.05,10.41/)![km] real 
!      real(r_kind), dimension(3) :: u00 = (/0.,0.,45./) !pass this point

!.....Cartesian ending point...for calculation ending at the same
!point,studying the effect of 1-cosTheta
      real(r_kind), dimension(3), parameter :: uf = (/0.,0.,44.4/) ![km]
      real(r_kind) :: nNorm

!.....output filename...................................................

      real(r_kind), dimension(3) :: un
      real(r_kind), parameter :: pi=3.1415926
      real(r_kind), parameter :: thetaMin=0
      real(r_kind) :: thetaMax
      real(r_kind) :: deltaTheta
      real(r_kind) :: phi
      integer :: j,k,p,pp,fN
      integer, parameter :: nTheta=50
      integer, parameter :: nStart=0
      character(len=100) :: format_str
      character(len=100) :: fileNumber
      character(len=100) :: format_num 
      thetaMax=pi/2
      deltaTheta=(thetaMax-thetaMin)/nTheta
      phi=pi/2
!.....normalize the direction vector....................................
!        un = u/sqrt(sum(u*u))
      do p=0,1
        if (p.eq.0) then
          pp=-1
        else
          pp=1
        endif
        do j=0,nStart
!          u0(2)=0+j*deltay  
!          un(1)=(uf(1)-u0(1))/sqrt((uf(1)-u0(1))**2+(uf(2)-u0(2))**2+(uf(3)-u0(3))**2)
!          un(2)=(uf(2)-u0(2))/sqrt((uf(1)-u0(1))**2+(uf(2)-u0(2))**2+(uf(3)-u0(3))**2)
!          un(3)=(uf(3)-u0(3))/sqrt((uf(1)-u0(1))**2+(uf(2)-u0(2))**2+(uf(3)-u0(3))**2)

!  ..........go over range of angle from a fixed starting point
          do k=0,nTheta
            un(1)=SIN(thetaMin+k*deltaTheta)*COS(pp*phi)
            un(2)=SIN(thetaMin+k*deltaTheta)*SIN(pp*phi)
            un(3)=COS(thetaMin+k*deltaTheta)
!            u0(1)=u00(1)+un(1)*20
!            u0(2)=u00(2)+un(2)*20
!            u0(3)=u00(3)+un(3)*20
!            un(1)=-un(1)
!            un(2)=-un(2)
!            un(3)=-un(3)
            fN= k+j*(nTheta+1)+p*(1+nStart)*(1+nTheta)
            if (fN < 10) then
                format_str = "(I1.1)"
            else if (fN<100) then
                format_str = "(I2.2)"
            else
                format_str = "(I3.3)"
            endif
            write (fileNumber,format_str)fN
            write (*,*) fN,pp,j,k,u0(1),u0(2),u0(3),un(1),un(2),un(3)
!  .....read-in the points of the trajectory..............................
            open(8,file='output/trajectory'//trim(fileNumber)//'.txt',status='unknown')
            do i=1,ntr
              if (u0(1)+un(1)*(i-1)*deltax>0 .AND. u0(1)+un(1)*(i-1)*deltax <10) then
                format_num = "(3X,F17.15)"
              else if (u0(1)+un(1)*(i-1)*deltax>0 .AND. u0(1)+un(1)*(i-1)*deltax<100) then
                format_num = "(3X,F17.14)"
              else if (u0(1)+un(1)*(i-1)*deltax>0) then
                format_num = "(3X,F17.13)"
              else if (u0(1)+un(1)*(i-1)*deltax<0 .AND. u0(1)+un(1)*(i-1)*deltax > -10) then
                format_num = "(3X,F17.14)"
              else if (u0(1)+un(1)*(i-1)*deltax<-10 .AND. u0(1)+un(1)*(i-1)*deltax > -100) then
                format_num = "(3X,F17.13)"
              else if (u0(1)+un(1)*(i-1)*deltax<0) then
                format_num = "(3X,F17.12)"
              else
                format_num = "(3X,F17.15)"
              endif
              write(8,format_num)u0(1)+un(1)*(i-1)*deltax
    
              if (u0(2)+un(2)*(i-1)*deltax>0 .AND. u0(2)+un(2)*(i-1)*deltax < 10) then
                format_num = "(3X,F17.15)"
              else if (u0(2)+un(2)*(i-1)*deltax>0 .AND. u0(2)+un(2)*(i-1)*deltax<100) then
                format_num = "(3X,F17.14)"
              else if (u0(2)+un(2)*(i-1)*deltax>0) then
                format_num = "(3X,F17.13)"
              else if (u0(2)+un(2)*(i-1)*deltax<0 .AND. u0(2)+un(2)*(i-1)*deltax >-10) then
                format_num = "(3X,F17.14)"
              else if (u0(2)+un(2)*(i-1)*deltax<-10 .AND. u0(2)+un(2)*(i-1)*deltax > -100) then
                format_num = "(3X,F17.13)"
              else if (u0(2)+un(2)*(i-1)*deltax<0) then
                format_num = "(3X,F17.12)"
              else
                format_num = "(3X,F17.15)"
              endif
              write(8,format_num)u0(2)+un(2)*(i-1)*deltax
    
              if (u0(3)+un(3)*(i-1)*deltax>0 .AND. u0(3)+un(3)*(i-1)*deltax < 10) then
                format_num = "(3X,F17.15)"
              else if (u0(3)+un(3)*(i-1)*deltax>0 .AND. u0(3)+un(3)*(i-1)*deltax<100) then
                format_num = "(3X,F17.14)"
              else if (u0(3)+un(3)*(i-1)*deltax>0) then
                format_num = "(3X,F17.13)"
              else if (u0(3)+un(3)*(i-1)*deltax<0 .AND. u0(3)+un(3)*(i-1)*deltax > -10) then
                format_num = "(3X,F17.14)"
              else if (u0(3)+un(3)*(i-1)*deltax<-10 .AND. u0(3)+un(3)*(i-1)*deltax > -100) then
                format_num = "(3X,F17.13)"
              else if (u0(3)+un(3)*(i-1)*deltax<0) then
                format_num = "(3X,F17.12)"
              else
                format_num = "(3X,F17.15)"
              endif
              write(8,format_num)u0(3)+un(3)*(i-1)*deltax
    
            end do
            close(8)
          end do
        end do
      end do
      end program main
