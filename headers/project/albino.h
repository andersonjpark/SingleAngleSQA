double get_rho(const double r){
  return exp(lnrho(log(r)));
}
double get_drhodr(const double rrho, const double r){
  return rrho*lnrho.Derivative(log(r))/r;
}
double get_Ye(const double r){
  return Ye(r);
}
double get_dYedr(const double r){
  return Ye.Derivative(r);
}

//=======//
// Ebins //
//=======//
void set_Ebins(vector<double>& E){
  const double NEP=8;
  cout<<"NE="<<NE << " NEP="<<NEP;
  E.resize(NE);
  for(int i=0;i<NE;i++){
    unsigned ind = i;
    if(NE==1||NE==2||NE==4) ind = i*NEP/NE+NEP/NE/2;

    E[i] = 1.e6*cgs::units::eV * exp(log(2*1000000) + ind*(log(37.48*1000000) - log(2*1000000))/(NE-1))/1000000;

    cout<<E[i]<<endl;
    cout.flush();
  }
  cout.flush();
}

//===================//
// Vacuum Potentials //
//===================//
double deltaV(const double E){ // erg
  return abs(dm21)*cgs::constants::c4 / (2.*E);
}

void set_kV(vector<vector<double> >& kV){
  assert(NF==2);
  for(int i=0;i<NE;i++){
    kV[i][0] = m1*m1 * cgs::constants::c4 /2./E[i];
    kV[i][1] = kV[i][0] + deltaV(E[i]);
  }
}

//===================//
// Matter Potentials //
//===================//
double Ve(double rho, double Ye){
  return (M_SQRT2*cgs::constants::GF/cgs::constants::Mp)*rho*Ye;
}

double dVedr(double rho, double drhodr, double Ye, double dYedr){
  return (M_SQRT2*cgs::constants::GF/cgs::constants::Mp) * (drhodr*Ye + rho*dYedr );
}

double Vmu(double rho, double Ye){ return 0.;}

double dVmudr(double rho, double drhodr, double Ye, double dYedr){ return 0.;}


//=============================//
// Self-Interaction Potentials //
//=============================//
void getP(vector<vector<MATRIX<complex<double>,NF,NF> > >& pmatrixf0, const double r)
{
  for(int i=0;i<=NE-1;i++){
    for(flavour f=e;f<=mu;f++)
      for(flavour fp=e; fp<=mu; fp++){
	pmatrixf0[matter    ][i][f][fp] = 0;
	pmatrixf0[antimatter][i][f][fp] = 0;
      }
    pmatrixf0[matter    ][i][e ][e ]=eP[i](r);
    pmatrixf0[antimatter][i][e ][e ]=eBarP[i](r);
    pmatrixf0[antimatter][i][mu][mu]=xP[i](r);
    pmatrixf0[matter    ][i][mu][mu]=xP[i](r);
  }
}
