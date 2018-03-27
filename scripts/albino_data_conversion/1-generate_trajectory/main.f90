program main
  
  implicit none
  
  integer, parameter :: r_kind = 8
  integer :: i
  
  !.....number of points..................................................
  integer, parameter :: ntr = 800 ! number of points
  real(r_kind), parameter :: deltax = 0.2  ![km]
  real(r_kind), dimension(3) :: u0 =(/0.,30.05,10.41/)![km] real 
  real(r_kind), dimension(3) :: un ! direction vector
  real(r_kind), parameter :: pi=3.1415926
  real(r_kind), parameter :: phi = pi/2.0
  real(r_kind), parameter :: theta = pi/4.0
  character(len=100) :: format_num 

  !.....normalize the direction vector....................................
  un(1)=SIN(theta)*COS(phi)
  un(2)=SIN(theta)*SIN(phi)
  un(3)=COS(theta)
  write (*,*) u0(1),u0(2),u0(3),un(1),un(2),un(3)
  !  .....read-in the points of the trajectory..............................
  open(8,file='output/trajectory.txt',status='unknown')
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
end program main
