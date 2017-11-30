double Ve(double rho, double Ye);
double dVedr(double rho, double drhodr, double Ye, double dYedr);

double Vmu(double rho, double Ye);
double dVmudr(double rho, double drhodr, double Ye, double dYedr);

// **************************************************************************

double Ve(double rho, double Ye){ return (M_SQRT2*cgs::constants::GF/cgs::constants::Mp)*rho*Ye;}

double dVedr(double rho, double drhodr, double Ye, double dYedr){ return (M_SQRT2*cgs::constants::GF/cgs::constants::Mp) * (drhodr*Ye + rho*dYedr );}

double Vmu(double rho, double Ye){ return 0.;}

double dVmudr(double rho, double drhodr, double Ye, double dYedr){ return 0.;}



