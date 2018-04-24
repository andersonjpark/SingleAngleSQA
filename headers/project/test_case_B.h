#include "../nulib_interface.h"
#include "../isospin.h"

void initialize(vector<vector<MATRIX<complex<double>,NF,NF> > >& fmatrixf,
		double r, double rho, double T, double Ye){
  for(int i=0; i<NE; i++){
    for(state m=matter; m<=antimatter; m++)
      for(flavour f1=e; f1<=mu; f1++)
	for(flavour f2=e; f2<=mu; f2++) 
	  fmatrixf[m][i][f1][f2] = 0;

    fmatrixf[    matter][i][e ][e ] = 1;//eas.Bnu(0,i);
    fmatrixf[    matter][i][mu][mu] = 0;//eas.Bnu(2,i);
    fmatrixf[antimatter][i][e ][e ] = 4./3.;//eas.Bnu(1,i);
    fmatrixf[antimatter][i][mu][mu] = 0;//eas.Bnu(2,i);
  }
}

double get_rho(const double r){
  return -1;
}
double get_drhodr(const double rrho, const double r){
  return -1;
}
double get_Ye(const double r){
  return -1;
}
double get_dYedr(const double r){
  return -1;
}

double deltaV(const double E){ // erg
  return abs(dm21)*cgs::constants::c4 / (2.*E);
}
//=======//
// Ebins //
//=======//
void set_Ebins(vector<double>& E){
  const double NEP=8;
  cout<<endl << "NE="<<NE << " NEP="<<NEP << endl;
  E.resize(NE);
  for(int i=0;i<NE;i++){
    unsigned ind = i;
    if(NE==1||NE==2||NE==4) ind = i*NEP/NE+NEP/NE/2;

    E[i] = 1.e6*cgs::units::eV * exp(log(2*1000000) + ind*(log(37.48*1000000) - log(2*1000000))/(NE-1))/1000000;

    cout<<E[i]<<" "<<deltaV(E[i]) << " " << (cgs::constants::hbarc * 2.*M_PI)/deltaV(E[i]) << endl;
    cout.flush();
  }
  cout.flush();
}

//===================//
// Vacuum Potentials //
//===================//

void set_kV(vector<vector<double> >& kV){
  assert(NF==2);
  for(int i=0;i<NE;i++){
    kV[i][0] = m1*m1 * cgs::constants::c4 /2./E[0];
    kV[i][1] = kV[i][0] + deltaV(E[0]);
  }
}

//===================//
// Matter Potentials //
//===================//
double Ve(double rho, double Ye){
  return 1000*deltaV(E[0]);
}

double dVedr(double rho, double drhodr, double Ye, double dYedr){
  return 0;
}

double Vmu(double rho, double Ye){ return 0.;}

double dVmudr(double rho, double drhodr, double Ye, double dYedr){ return 0.;}


//=============================//
// Self-Interaction Potentials //
//=============================//
double MU(const double r, const double E){ // erg
  double dV = deltaV(E);
  double r_dimless = r*dV / (cgs::constants::hbarc * 2.*M_PI);
  return 1e4 * dV * exp(-r_dimless / 10.)  / (double)NE;
}
void getPunosc(const double r, const state m, const unsigned ig,
	       MATRIX<complex<double>,NF,NF>& p_unosc){
  
  // construct potential from distribution function 
  for(flavour f=e;f<=mu;f++)
    for(flavour fp=e; fp<=mu; fp++)
      p_unosc[f][fp] = 0;
  p_unosc[e][e] = MU(r, E[0]) * (m==matter ? 1. : 4./3.);
}

void my_interact(vector<vector<MATRIX<complex<double>,NF,NF> > >& fmatrixf,
	      double rho, double T, double Ye, double r, double dr){
  double kappa = deltaV(E[0]) / (cgs::constants::hbarc * M_2PI)/100.;
  double tmp = 0;
  for(int i=0; i<NE; i++){
    fmatrixf[    matter][i][e][e ] *= exp(-kappa    * dr);
    fmatrixf[    matter][i][e][mu] *= exp(-kappa/2. * dr);
    fmatrixf[    matter][i][mu][e] = conj(fmatrixf[matter][i][e][mu]);
  }
}
