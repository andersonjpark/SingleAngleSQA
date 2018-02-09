double deltaV(const double E){ // 1/cm
  return abs(dm21)*cgs::constants::c4 / (2.*E) / (cgs::constants::hbarc*2.*M_PI);
}

double Ve(double rho, double Ye);
double dVedr(double rho, double drhodr, double Ye, double dYedr);

double Vmu(double rho, double Ye);
double dVmudr(double rho, double drhodr, double Ye, double dYedr);

// **************************************************************************

double Ve(double rho, double Ye){
  return 1000*deltaV(E[0])*cgs::constants::hbarc*2.*M_PI;//
  //UNCOMMENT return (M_SQRT2*cgs::constants::GF/cgs::constants::Mp)*rho*Ye;
}

double dVedr(double rho, double drhodr, double Ye, double dYedr){
  return 0;//
  //UNCOMMENT return (M_SQRT2*cgs::constants::GF/cgs::constants::Mp) * (drhodr*Ye + rho*dYedr );
}

double Vmu(double rho, double Ye){ return 0.;}

double dVmudr(double rho, double drhodr, double Ye, double dYedr){ return 0.;}



