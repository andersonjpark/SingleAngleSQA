#include "../nulib_interface.h"
#include "../isospin.h"

void initialize(vector<vector<MATRIX<complex<double>,NF,NF> > >& fmatrixf,
		double r, double rho, double T, double Ye){
  // T should be MeV
  double T_tmp = 10.0;
  nulibtable_range_species_range_energy_(&rho, &T_tmp, &Ye, &eas.storage.front(),
  					 &__nulibtable_MOD_nulibtable_number_species,
  					 &__nulibtable_MOD_nulibtable_number_groups,
  					 &__nulibtable_MOD_nulibtable_number_easvariables);
  eas.fix_units();
  
  for(int i=0; i<NE; i++){
    for(state m=matter; m<=antimatter; m++)
      for(flavour f1=e; f1<=mu; f1++)
	for(flavour f2=e; f2<=mu; f2++) 
	  fmatrixf[m][i][f1][f2] = 0;

    double Elow = i>0    ? E[i-1] : 0;
    double Ehi  = i<NE-1 ? E[i+1] : E[i] + 0.5*(E[i]-E[i-1]);
    double phaseSpaceVol = 4./3.*M_PI * (pow(Ehi,3) - pow(Elow,3)) / pow(2.*M_PI*cgs::constants::hbarc,3);
    fmatrixf[    matter][i][e ][e ] = eD[i](r)    / phaseSpaceVol; //eas.emis(0,i) / eas.abs(0,i);
    fmatrixf[    matter][i][mu][mu] = xD[i](r)    / phaseSpaceVol; //eas.emis(2,i) / eas.abs(2,i);
    fmatrixf[antimatter][i][e ][e ] = eBarD[i](r) / phaseSpaceVol; //eas.emis(1,i) / eas.abs(1,i);
    fmatrixf[antimatter][i][mu][mu] = xD[i](r)    / phaseSpaceVol; //eas.emis(2,i) / eas.abs(2,i);
  }  
}

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
void getP(const double r,
	  const vector<vector<MATRIX<complex<double>,NF,NF> > > U0, 
	  const vector<vector<MATRIX<complex<double>,NF,NF> > > fmatrixf, 
	  vector<vector<MATRIX<complex<double>,NF,NF> > >& pmatrixf0,
	  vector<vector<MATRIX<complex<double>,NF,NF> > >& pmatrixm0){

  double hf[4]; // [state][pauli index] coefficients of Pauli matrices for f
  double hp[4]; // pauli matrix coefficients for p
  double hf_norm[3]; // components of isospin spatial unit vector
  double hp_unosc[4]; // unoscillated potential has only diagonal elements
  MATRIX<complex<double>,NF,NF> p_unosc;
  
  for(int i=0;i<=NE-1;i++){
    for(state m=matter; m<=antimatter; m++){

      // decompose (oscillated) distribution function
      pauli_decompose(fmatrixf[m][i], hf);
      double hf_length = sqrt(hf[0]*hf[0] + hf[1]*hf[1] + hf[2]*hf[2]);
      for(unsigned k=0; k<3; k++) hf_norm[k] = hf[k] / hf_length;

      // decompose unoscillated potential
      double P0 = (m==matter ? eP[i](r) : eBarP[i](r));
      double P1 = xP[i](r);
      p_unosc[0][0] = complex<double>(P0,0);
      p_unosc[1][0] = complex<double>(0,0);
      p_unosc[0][1] = complex<double>(0,0);
      p_unosc[1][1] = complex<double>(P1,0);
      pauli_decompose(p_unosc, hp_unosc);

      // re-distribute p Pauli coefficients to have same flavor angle as f
      assert(hp_unosc[0] == 0);
      assert(hp_unosc[1] == 0);
      for(unsigned k=0; k<3; k++)
	hp[k] = hf_norm[k] * abs(hp_unosc[2]);
      hp[3] = hp_unosc[3];

      // reconstruct the potential matrix
      pauli_reconstruct(hp, pmatrixf0[m][i]);

      // put in mass basis
      pmatrixm0[m][i] = Adjoint(U0[m][i])
	* pmatrixf0[m][i]
	* U0[m][i];
    }
  }
}

void interact(vector<vector<MATRIX<complex<double>,NF,NF> > >& fmatrixf,
	      double rho, double T, double Ye, double dr){
  //return;
  // don't do anything if too sparse
  //if(log10(rho) <= __nulibtable_MOD_nulibtable_logrho_min)
  //  return;

  // T should be MeV
  double T_tmp = 10.0;
  //nulibtable_range_species_range_energy_(&rho, &T_tmp, &Ye, &eas.storage.front(),
  //					 &__nulibtable_MOD_nulibtable_number_species,
  //					 &__nulibtable_MOD_nulibtable_number_groups,
  //					 &__nulibtable_MOD_nulibtable_number_easvariables);
  //eas.fix_units();


  double absopac = 1./(100.*1e5);
  double tmp = 0;
  for(int i=0; i<NE; i++){
    // scale the diagonal components
    fmatrixf[    matter][i][e ][e ] += (/*eas.emis(0,i)*/ - absopac/*eas.abs(0,i)*/*fmatrixf[    matter][i][e ][e ]) * dr;
    fmatrixf[    matter][i][mu][mu] += (/*eas.emis(2,i)*/ - absopac/*eas.abs(2,i)*/*fmatrixf[    matter][i][mu][mu]) * dr;
    fmatrixf[antimatter][i][e ][e ] += (/*eas.emis(1,i)*/ - absopac/*eas.abs(1,i)*/*fmatrixf[antimatter][i][e ][e ]) * dr;
    fmatrixf[antimatter][i][mu][mu] += (/*eas.emis(2,i)*/ - absopac/*eas.abs(2,1)*/*fmatrixf[antimatter][i][mu][mu]) * dr;

    // scale the off-diagonal components
    double kappa_e  = absopac; //eas.abs(0,i);
    double kappa_mu = absopac; //eas.abs(2,i);
    double kappa_avg = 0.5 * (kappa_e+kappa_mu);
    //fmatrixf[    matter][i][e ][e ] *= exp(kappa_e   * dr);
    //fmatrixf[    matter][i][mu][mu] *= exp(kappa_mu  * dr);
    fmatrixf[    matter][i][e ][mu] -= kappa_avg * dr; //*= exp(kappa_avg * dr);
    fmatrixf[    matter][i][mu][e ] = conj(fmatrixf[matter][i][e][mu]);
    fmatrixf[antimatter][i][e ][mu] -= kappa_avg * dr; //*= exp(kappa_avg * dr);
    fmatrixf[antimatter][i][mu][e ] = conj(fmatrixf[matter][i][e][mu]);
  }
}
