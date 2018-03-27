program main
  
  implicit none
  
  integer, parameter :: r_kind = 8
  integer :: i
  
  !.....number of points..................................................
  integer, parameter :: ntr = 800 ! number of points
  real(r_kind), parameter :: deltax = 0.2  ![km]
  real(r_kind), dimension(3) :: u0 =(/0.,0.,18./)![km] real 
  real(r_kind), dimension(3) :: un =(/0.,sqrt(2.),3./)! direction vector
  character(len=100) :: format_num 

  !.....normalize the direction vector....................................
  un = un/sqrt(sum(un*un))
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
