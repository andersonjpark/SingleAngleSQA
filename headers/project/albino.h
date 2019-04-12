#ifndef ALBINO_H
#define ALBINO_H

#include "../nulib_interface.h"
#include "../isospin.h"

//=======//
// Ebins //
//=======//
double dE3(const array<double,NE>& E, const unsigned i){
  double dlogE = (log(E[NE-1]) - log(E[0])) / (NE-1.);
  double Elow = exp(log(E[0]) + (i-0.5)*dlogE);
  double Ehi  = exp(log(E[0]) + (i+0.5)*dlogE);
  return pow(Ehi,3) - pow(Elow,3);
}
array<double,NE> set_Ebins(){
  array<double,NE> E;

  const double NEP=8;
  cout << endl;
  cout<<"NE="<<NE << " NEP="<<NEP << endl;
  for(int i=0;i<NE;i++){
    unsigned ind = i;
    if(NE==1||NE==2||NE==4) ind = i*NEP/NE+NEP/NE/2;

    double lEtop = log(37.48 * 1e6*cgs::units::eV) ; //erg
    double lEbottom = log(2. * 1e6*cgs::units::eV) ; //erg
    double dlE = (lEtop-lEbottom)/(NE-1);
    E[i] =  exp(lEbottom + ind*dlE);

    cout << E[i]/(1.e6*cgs::units::eV) << " ";
    cout << exp(lEbottom + (ind-0.5)*dlE)/(1.e6*cgs::units::eV) << " ";
    cout << exp(lEbottom + (ind+0.5)*dlE)/(1.e6*cgs::units::eV) << endl;
  }
  cout.flush();

  return E;
}
template<typename T>
T phaseVolDensity(const array<double,NE>& E, const T density, const unsigned i){
  double phaseSpaceVol = 4.*M_PI * dE3(E,i)/3. / pow(2.*M_PI*cgs::constants::hbarc,3);
  return density / phaseSpaceVol;
}

//============//
// Initialize //
//============//
void initialize(State& s,
		double r,
		array<double,NE> E,
		const array<array<array<DISCONTINUOUS,NF>,NE>,NM>& D_unosc){
  // T should be MeV
  cout << "Setting initial data." << endl;
  cout << "rho = " << s.rho << " g/ccm" << endl;
  cout << "T = " << s.T << " MeV" << endl;
  cout << "Ye = " << s.Ye << endl;
  s.Ye = max(s.Ye,__nulibtable_MOD_nulibtable_ye_min);
  /* nulibtable_range_species_range_energy_(&s.rho, &s.T, &s.Ye, &eas.eas.front(), */
  /* 					 &__nulibtable_MOD_nulibtable_number_species, */
  /* 					 &__nulibtable_MOD_nulibtable_number_groups, */
  /* 					 &__nulibtable_MOD_nulibtable_number_easvariables); */
  
  for(int i=0; i<NE; i++){
    for(state m=matter; m<=antimatter; m++){
      s.fmatrixf[m][i] = MATRIX<complex<double>,NF,NF>();
      for(flavour f=e; f<=mu; f++)
	s.fmatrixf[m][i][f][f] = phaseVolDensity(E, D_unosc[m][i][f](s.r), i);
    }
    
    cout << "GROUP " << i << endl;
    /* cout << "\teas.emis = {" << eas.emis(0,i) << ", " << eas.emis(1,i) << ", " << eas.emis(2,i) << "}" << endl; */
    /* cout << "\teas.abs = {" << eas.abs(0,i) << ", " << eas.abs(1,i) << ", " << eas.abs(2,i) << "}" << endl; */
    /* cout << "\tBB = {" << eas.emis(0,i)/eas.abs(0,i) << ", " << eas.emis(1,i)/eas.abs(1,i) << ", " << eas.emis(2,i)/eas.abs(2,i) << "}" << endl; */
    /* cout << "\teas.scat = {" << eas.scat(0,i) << ", " << eas.scat(1,i) << ", " << eas.scat(2,i) << "}" << endl; */

    cout << "\tf = {" << real(s.fmatrixf[matter][i][e][e]) << ", " << real(s.fmatrixf[antimatter][i][e][e]) << ", " << real(s.fmatrixf[matter][i][mu][mu]) << ", " << real(s.fmatrixf[antimatter][i][mu][mu]) << "}" << endl;
  }
  
  for(state m=matter; m<=antimatter; m++)
    for(int i=0; i<NE; i++)
      for(flavour f1=e; f1<=mu; f1++)
	for(flavour f2=e; f2<=mu; f2++)
	  assert(s.fmatrixf[m][i][f1][f2] == s.fmatrixf[m][i][f1][f2]);
}


//===================//
// Vacuum Potentials //
//===================//
double deltaV(const double E){ // erg
  return abs(dm21)*cgs::constants::c4 / (2.*E);
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


array<array<MATRIX<complex<double>,NF,NF>,NE>,NM>
  Kinteract(const State& s, const array<array<array<DISCONTINUOUS,NF>,NE>,NM>& D_unosc){

  // set up the array to be returned
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> dfdr;

  // don't do anything if too sparse
  if(log10(s.rho) <= __nulibtable_MOD_nulibtable_logrho_min)
    return dfdr;

  // T should be MeV
  nulibtable_range_species_range_energy_(&s.rho, &s.T, &s.Ye, &eas.eas.front(),
  					 &__nulibtable_MOD_nulibtable_number_species,
  					 &__nulibtable_MOD_nulibtable_number_groups,
  					 &__nulibtable_MOD_nulibtable_number_easvariables);

  // get oscillated background density
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> DBackground;
  for(state m=matter; m<=antimatter; m++){
    for(int ig=0; ig<NE; ig++){
      for(flavour f=e; f<=mu; f++)
	DBackground[m][ig][f][f] = D_unosc[m][ig][f](s.r);
      DBackground[m][ig] = s.Sf[m][ig] * DBackground[m][ig] * Adjoint(s.Sf[m][ig]);
    }
  }

#pragma omp parallel for collapse(2)
  for(int m=matter; m<=antimatter; m++){
    for(int i=0; i<NE; i++){

      // get nulib species indices
      const int se = (m==matter ? 0 : 1);
      const int sx = (m==matter ? 2 : 3);
      const state mbar = (m==matter ? antimatter : matter);
      
      // intermediate variables
      complex<double> unblock_in, unblock_out;
      double kappa_e, kappa_mu, kappa_avg;
      double Phi0e, Phi0x, Phi0_avg, Phi0_tilde, Phi0_emu;
      MATRIX<complex<double>,NF,NF> block, Pi_plus, Pi_minus;
      MATRIX<double,NF,NF> Phi0, Phi0avg,  Phi0tilde;      

      // emission
      dfdr[m][i][e ][e ] += eas.emis(se,i);
      dfdr[m][i][mu][mu] += eas.emis(sx,i);
      
      // absorption kappa_abs is <kappa> for absorption
      MATRIX<double,NF,NF> kappa_abs = eas.avg_matrix(eas.abs(se,i), eas.abs(sx,i));
      for(flavour f1=e; f1<=mu; f1++)
	for(flavour f2=e; f2<=mu; f2++)
	  dfdr[m][i][f1][f2] -= Phi0avg[f1][f2] * s.fmatrixf[m][i][f1][f2];
    
      // scattering and pair annihilation
      // no factor of 1/2 in front of Phi0 because
      // integrated over outgoing theta, assumed isotropic
      for(int j=0; j<NE; j++){

	// in-scattering from j to i. D*Phi0 give density scattered
	// divide by phase space vol in i to get f
	Phi0avg      = eas.avg_matrix(  eas.Phi0(0,j,i), eas.Phi0(2,j,i));
	Phi0tilde    = eas.tilde_matrix(eas.Phi0(0,j,i), eas.Phi0(2,j,i));
	Phi0     = Phi0avg    - Phi0tilde;
	block    = eas.blocking_term0(Phi0, s.fmatrixf[m][i], DBackground[m][j]);
	for(flavour f1=e; f1<=mu; f1++)
	  for(flavour f2=e; f2<=mu; f2++)
	    dfdr[m][i][f1][f2] += phaseVolDensity(s.E, DBackground[m][j][f1][f2]*Phi0[f1][f2] - block[f1][f2], i);

	// out-scattering from i to j. for blocking, get phase space vol from D[j] in j
	Phi0avg      = eas.avg_matrix(  eas.Phi0(0,i,j), eas.Phi0(2,i,j));
	Phi0tilde    = eas.tilde_matrix(eas.Phi0(0,i,j), eas.Phi0(2,i,j));
	Phi0    = Phi0avg    - Phi0tilde;
	block    = eas.blocking_term0(Phi0, s.fmatrixf[m][i], DBackground[m][j] );
	for(flavour f1=e; f1<=mu; f1++)
	  for(flavour f2=e; f2<=mu; f2++)
	    dfdr[m][i][f1][f2] += s.fmatrixf[m][i][f1][f2]*Phi0avg[f1][f2] - phaseVolDensity(s.E, block[f1][f2], j);

	// Make sure dfdr is Hermitian
	dfdr[m][i][mu][e ] = conj(dfdr[m][i][e][mu]);
	
      } // other group
    } // group
  } // state
  
  return dfdr;
}


#endif
