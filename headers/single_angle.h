double CSI(double r);
double dCSIdr(double r);

// *********************************************************************

// geometric functions for single angle self interaction- first is from Duan et al, second from Dighe et al

//double CSI(double r){ return 0.5*pow(1.-sqrt(1.-Rnu/r)*sqrt(1.+Rnu/r),2.);}
//double dCSIdr(double r){ return -pow(Rnu/r,2.)/r*( 1./sqrt(1.-Rnu/r)/sqrt(1.+Rnu/r) -1.);}

double CSI(double r){ double y=r/Rnu; return 2.*pow(y-sqrt(y-1.)*sqrt(y+1.),2.)-0.5/(y*y);}
double dCSIdr(double r){ double y=r/Rnu, dydr=1./Rnu; return dydr*0.5*(-8.*pow(y-sqrt(y-1.)*sqrt(y+1.),2.)/sqrt(y-1.)/sqrt(y+1.)+2./(y*y*y));}
