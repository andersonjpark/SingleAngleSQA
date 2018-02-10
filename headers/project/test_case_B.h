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
double ALPHA(const double r){
  return 4./3.;
}
void getP(vector<vector<MATRIX<complex<double>,NF,NF> > >& pmatrixf0, const double r)
{
  for(int i=0;i<=NE-1;i++){
    for(flavour f=e;f<=mu;f++)
      for(flavour fp=e; fp<=mu; fp++){
	pmatrixf0[matter    ][i][f][fp] = 0;
	pmatrixf0[antimatter][i][f][fp] = 0;
      }
    pmatrixf0[matter    ][i][e ][e ]=MU(r, E[0]);
    pmatrixf0[antimatter][i][e ][e ]=MU(r, E[0]) * ALPHA(r);
    pmatrixf0[antimatter][i][mu][mu]=0;
    pmatrixf0[matter    ][i][mu][mu]=0;
  }
}
