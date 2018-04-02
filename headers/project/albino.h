#include "../nulib_interface.h"
#include "../isospin.h"


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
double dE(const unsigned i){
  double dlogE = (log(E[NE-1]) - log(E[0])) / (NE-1.);
  double Elow = exp(log(E[0]) + (i-0.5)*dlogE);
  double Ehi  = exp(log(E[0]) + (i+0.5)*dlogE);
  return Ehi - Elow;
}
double phaseVolDensity(const double density, const unsigned i){
  double phaseSpaceVol = 4.*M_PI * E[i]*E[i]*dE(i) / pow(2.*M_PI*cgs::constants::hbarc,3);
  return density / phaseSpaceVol;
}

//============//
// Initialize //
//============//
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

    fmatrixf[    matter][i][e ][e ] = phaseVolDensity(eD[i](r)   , i); //eas.emis(0,i) / eas.abs(0,i);
    fmatrixf[    matter][i][mu][mu] = phaseVolDensity(xD[i](r)   , i); //eas.emis(2,i) / eas.abs(2,i);
    fmatrixf[antimatter][i][e ][e ] = phaseVolDensity(eBarD[i](r), i); //eas.emis(1,i) / eas.abs(1,i);
    fmatrixf[antimatter][i][mu][mu] = phaseVolDensity(xD[i](r)   , i); //eas.emis(2,i) / eas.abs(2,i);
  }
  
  for(state m=matter; m<=antimatter; m++)
    for(int i=0; i<NE; i++)
      for(flavour f1=e; f1<=mu; f1++)
	for(flavour f2=e; f2<=mu; f2++)
	  assert(fmatrixf[m][i][f1][f2] == fmatrixf[m][i][f1][f2]);
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
	      double rho, double T, double Ye, double r, double dr){
  //return;

  // set up rate matrix
  MATRIX<complex<double>,NF,NF> dfdr, dfbardr;
  
  // don't do anything if too sparse
  if(log10(rho) <= __nulibtable_MOD_nulibtable_logrho_min)
    return;

  // T should be MeV
  double T_tmp = 10.0;
  nulibtable_range_species_range_energy_(&rho, &T_tmp, &Ye, &eas.storage.front(),
  					 &__nulibtable_MOD_nulibtable_number_species,
  					 &__nulibtable_MOD_nulibtable_number_groups,
  					 &__nulibtable_MOD_nulibtable_number_easvariables);
  eas.fix_units();


  double tmp = 0;
  for(int i=0; i<NE; i++){
    // reset dfdr
    for(flavour f1=e; f1<=mu; f1++)
      for(flavour f2=e; f2<=mu; f2++){
	dfdr   [f1][f2] = 0;
	dfbardr[f1][f2] = 0;
      }

    // get interaction rates
    /* double absopac = eas.abs(1./(500.*1e5); */
    /* double scatopac = 1./(100.*1e5); */
    /* double rmin = 50e5; */
    /* double emis_e = eas.absopac;//phaseVolDensity(eD[i](rmin)   , i) * absopac; */
    /* double emis_a = phaseVolDensity(eBarD[i](rmin), i) * absopac; */
    /* double emis_x = phaseVolDensity(xD[i](rmin)   , i) * absopac; */

    // absorption and out-scattering
    double kappa_e    = eas.abs(0,i) + eas.scat(0,i); //absopac + scatopac; //
    double kappa_ebar = eas.abs(1,i) + eas.scat(1,i); //absopac + scatopac; //
    double kappa_mu   = eas.abs(2,i) + eas.scat(2,i); //absopac; //
    double kappa_avg    = 0.5*(kappa_e   +kappa_mu);
    double kappa_avgbar = 0.5*(kappa_ebar+kappa_mu);
    dfdr   [e ][e ] -= kappa_e      * fmatrixf[    matter][i][e ][e ];
    dfbardr[e ][e ] -= kappa_ebar   * fmatrixf[antimatter][i][e ][e ];
    dfdr   [mu][mu] -= kappa_mu     * fmatrixf[    matter][i][mu][mu];
    dfbardr[mu][mu] -= kappa_mu     * fmatrixf[antimatter][i][mu][mu];
    dfdr   [e ][mu] -= kappa_avg    * fmatrixf[    matter][i][e ][mu];
    dfbardr[e ][mu] -= kappa_avgbar * fmatrixf[antimatter][i][e ][mu];

    // emission
    dfdr   [e ][e ] += eas.emis(0,i);//emis_e;
    dfbardr[e ][e ] += eas.emis(1,i);//emis_a;
    dfdr   [mu][mu] += eas.emis(2,i);//emis_x;
    dfbardr[mu][mu] += eas.emis(2,i);//emis_x;
    
    // in-scattering (currently assumes charged-current scattering)
    dfdr   [e ][e ] += phaseVolDensity(eD[i](r)   , i) * eas.scat(0,i);//scatopac;
    dfbardr[e ][e ] += phaseVolDensity(eBarD[i](r), i) * eas.scat(1,i);//scatopac;

    // Make sure dfdr is Hermitian
    dfdr   [mu][e ] = conj(dfdr   [e][mu]);
    dfbardr[mu][e ] = conj(dfbardr[e][mu]);

    // update fmatrixf
    for(flavour f1=e; f1<=mu; f1++)
      for(flavour f2=e; f2<=mu; f2++){
	fmatrixf[    matter][i][f1][f2] += dfdr[f1][f2]    * dr;
	fmatrixf[antimatter][i][f1][f2] += dfbardr[f1][f2] * dr;
      }

    // check that everything makes sense
    for(state s=matter; s<=antimatter; s++)
      for(flavour f1=e; f1<=mu; f1++)
	for(flavour f2=e; f2<=mu; f2++){
	  assert(fmatrixf[s][i][f1][f2] == fmatrixf[s][i][f1][f2]);
      }
  }
}
