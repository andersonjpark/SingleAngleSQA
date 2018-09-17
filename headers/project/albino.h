#include "../nulib_interface.h"
#include "../isospin.h"


//=======//
// Ebins //
//=======//
double dE3(const unsigned i){
  double dlogE = (log(E[NE-1]) - log(E[0])) / (NE-1.);
  double Elow = exp(log(E[0]) + (i-0.5)*dlogE);
  double Ehi  = exp(log(E[0]) + (i+0.5)*dlogE);
  return pow(Ehi,3) - pow(Elow,3);
}
void set_Ebins(vector<double>& E){
  const double NEP=8;
  E.resize(NE);
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
}
template<typename T>
T phaseVolDensity(const T density, const unsigned i){
  double phaseSpaceVol = 4.*M_PI * dE3(i)/3. / pow(2.*M_PI*cgs::constants::hbarc,3);
  return density / phaseSpaceVol;
}

//============//
// Initialize //
//============//
void initialize(vector<vector<MATRIX<complex<double>,NF,NF> > >& fmatrixf,
		double r, double rho, double T, double Ye){
  // T should be MeV
  cout << "Setting initial data." << endl;
  cout << "rho = " << rho << " g/ccm" << endl;
  cout << "T = " << T << " MeV" << endl;
  cout << "Ye = " << Ye << endl;
  double Ye_in = max(Ye,__nulibtable_MOD_nulibtable_ye_min);
  nulibtable_range_species_range_energy_(&rho, &T, &Ye_in, &eas.eas.front(),
  					 &__nulibtable_MOD_nulibtable_number_species,
  					 &__nulibtable_MOD_nulibtable_number_groups,
  					 &__nulibtable_MOD_nulibtable_number_easvariables);
  
  for(int i=0; i<NE; i++){
    for(state m=matter; m<=antimatter; m++)
      for(flavour f1=e; f1<=mu; f1++)
	for(flavour f2=e; f2<=mu; f2++) 
	  fmatrixf[m][i][f1][f2] = 0;
    double De = phaseVolDensity(eD[i](r)   , i); //eas.emis(0,i) / eas.abs(0,i); //
    double Da = phaseVolDensity(eBarD[i](r), i); //eas.emis(1,i) / eas.abs(1,i); //
    double Dx = phaseVolDensity(xD[i](r)   , i); //eas.emis(2,i) / eas.abs(2,i); //
    
    fmatrixf[    matter][i][e ][e ] = De; 
    fmatrixf[    matter][i][mu][mu] = Dx;
    fmatrixf[antimatter][i][e ][e ] = Da;
    fmatrixf[antimatter][i][mu][mu] = Dx;
      
    cout << "GROUP " << i << endl;
    cout << "\teas.emis = {" << eas.emis(0,i) << ", " << eas.emis(1,i) << ", " << eas.emis(2,i) << "}" << endl;
    cout << "\teas.abs = {" << eas.abs(0,i) << ", " << eas.abs(1,i) << ", " << eas.abs(2,i) << "}" << endl;
    cout << "\tBB = {" << eas.emis(0,i)/eas.abs(0,i) << ", " << eas.emis(1,i)/eas.abs(1,i) << ", " << eas.emis(2,i)/eas.abs(2,i) << "}" << endl;
    cout << "\teas.scat = {" << eas.scat(0,i) << ", " << eas.scat(1,i) << ", " << eas.scat(2,i) << "}" << endl;

    cout << "\tf = {" << real(fmatrixf[matter][i][e][e]) << ", " << real(fmatrixf[antimatter][i][e][e]) << ", " << real(fmatrixf[matter][i][mu][mu]) << ", " << real(fmatrixf[antimatter][i][mu][mu]) << "}" << endl;
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
void getPunosc(const double r, const state m, const unsigned ig,
	       MATRIX<complex<double>,NF,NF>& p_unosc){
  
  // decompose unoscillated potential
  double P0 = (m==matter ? eP[ig](r) : eBarP[ig](r));
  double P1 = xP[ig](r);
  p_unosc[e ][e ] = complex<double>(P0,0);
  p_unosc[mu][e ] = complex<double>(0,0);
  p_unosc[e ][mu] = complex<double>(0,0);
  p_unosc[mu][mu] = complex<double>(P1,0);
}


void my_interact(vector<vector<MATRIX<complex<double>,NF,NF> > >& fmatrixf,
		 vector<vector<MATRIX<complex<double>,NF,NF> > >& Scumulative,
		 double rho, double T, double Ye, double r, double dr){
  // don't do anything if too sparse
  if(log10(rho) <= __nulibtable_MOD_nulibtable_logrho_min)
    return;
  
  // set up rate matrix
  MATRIX<complex<double>,NF,NF> dfdr, dfbardr, block, blockbar;
  MATRIX<double,NF,NF> Phi0avg, Phi0tilde, Phi0avgbar, Phi0tildebar, Phi0, Phi0bar;
  vector<MATRIX<complex<double>,NF,NF>> DBackground, DbarBackground;
  DBackground.resize(eas.ng);
  DbarBackground.resize(eas.ng);

  // T should be MeV
  nulibtable_range_species_range_energy_(&rho, &T, &Ye, &eas.eas.front(),
  					 &__nulibtable_MOD_nulibtable_number_species,
  					 &__nulibtable_MOD_nulibtable_number_groups,
  					 &__nulibtable_MOD_nulibtable_number_easvariables);

  // get background density
  for(int ig=0; ig<NE; ig++){
    DBackground[ig][e ][e ]    =    eD[ig](r);
    DBackground[ig][mu][mu]    =    xD[ig](r);
    DbarBackground[ig][e ][e ] = eBarD[ig](r);
    DbarBackground[ig][mu][mu] =    xD[ig](r);
    DBackground[ig][e][mu]     = 0;
    DBackground[ig][mu][e]     = 0;
    DbarBackground[ig][e][mu]  = 0;
    DbarBackground[ig][mu][e]  = 0;
    DBackground[ig] = U0[matter][ig]
      * Scumulative[matter][ig]
      * Adjoint(U0[matter][ig])
      * DBackground[ig]
      * U0[matter][ig]
      * Adjoint(Scumulative[matter][ig])
      * Adjoint(U0[matter][ig]);
    DbarBackground[ig] = U0[antimatter][ig]
      * Scumulative[antimatter][ig]
      * Adjoint(U0[antimatter][ig])
      * DbarBackground[ig]
      * U0[antimatter][ig]
      * Adjoint(Scumulative[antimatter][ig])
      * Adjoint(U0[antimatter][ig]);
  }

  double kappa_e, kappa_ebar, kappa_mu, kappa_avg, kappa_avgbar;
  double Phi0_e, Phi0_ebar, Phi0_mu, Phi0_avg, Phi0_avgbar, Phi0_tilde, Phi0_tildebar, Phi0_emu, Phi0_emubar;
  for(int i=0; i<NE; i++){
    // reset dfdr
    for(flavour f1=e; f1<=mu; f1++)
      for(flavour f2=e; f2<=mu; f2++){
	dfdr   [f1][f2] = 0;
	dfbardr[f1][f2] = 0;
      }
    
    // emission
    dfdr   [e ][e ] += eas.emis(0,i);
    dfbardr[e ][e ] += eas.emis(1,i);
    dfdr   [mu][mu] += eas.emis(2,i);
    dfbardr[mu][mu] += eas.emis(2,i);

    // absorption
    Phi0avg    = eas.avg_matrix(eas.abs(0,i), eas.abs(2,i));
    Phi0avgbar = eas.avg_matrix(eas.abs(1,i), eas.abs(2,i));
    for(flavour f1=e; f1<=mu; f1++)
      for(flavour f2=e; f2<=mu; f2++){
	dfdr   [f1][f2] -= Phi0avg   [f1][f2] * fmatrixf[    matter][i][f1][f2];
	dfbardr[f1][f2] -= Phi0avgbar[f1][f2] * fmatrixf[antimatter][i][f1][f2];
      }

    // scattering
    // no factor of 1/2 in front of Phi0 because integrated over outgoing theta, assumed isotropic
    for(int j=0; j<NE; j++){
      // in-scattering from j to i. D*Phi0 give density scattered, divide by phase space vol in i to get f
      Phi0avg      = eas.avg_matrix(  eas.Phi0(0,j,i), eas.Phi0(2,j,i));
      Phi0avgbar   = eas.avg_matrix(  eas.Phi0(1,j,i), eas.Phi0(2,j,i));
      Phi0tilde    = eas.tilde_matrix(eas.Phi0(0,j,i), eas.Phi0(2,j,i));
      Phi0tildebar = eas.tilde_matrix(eas.Phi0(1,j,i), eas.Phi0(2,j,i));
      Phi0     = Phi0avg    - Phi0tilde;
      Phi0bar  = Phi0avgbar - Phi0tildebar;
      block    = eas.blocking_term0(Phi0,    fmatrixf[    matter][i],    DBackground[j]);
      blockbar = eas.blocking_term0(Phi0bar, fmatrixf[antimatter][i], DbarBackground[j]);
      
      for(flavour f1=e; f1<=mu; f1++)
	for(flavour f2=e; f2<=mu; f2++){
	  dfdr   [f1][f2] += phaseVolDensity(   DBackground[j][f1][f2]*Phi0   [f1][f2] - block[f1][f2], i);
	  dfbardr[f1][f2] += phaseVolDensity(DbarBackground[j][f1][f2]*Phi0bar[f1][f2] - block[f1][f2], i);
	}

      // out-scattering from i to j. for blocking, get phase space vol from D[j] in j
      Phi0avg      = eas.avg_matrix(eas.Phi0(0,i,j), eas.Phi0(2,i,j));
      Phi0avgbar   = eas.avg_matrix(eas.Phi0(1,i,j), eas.Phi0(2,i,j));
      Phi0tilde    = eas.avg_matrix(eas.Phi0(0,i,j), eas.Phi0(2,i,j));
      Phi0tildebar = eas.avg_matrix(eas.Phi0(1,i,j), eas.Phi0(2,i,j));
      Phi0    = Phi0avg    - Phi0tilde;
      Phi0bar = Phi0avgbar - Phi0tildebar;
      block    = eas.blocking_term0(Phi0   , fmatrixf[    matter][i],    DBackground[j] );
      blockbar = eas.blocking_term0(Phi0bar, fmatrixf[antimatter][i], DbarBackground[j] );
      for(flavour f1=e; f1<=mu; f1++)
	for(flavour f2=e; f2<=mu; f2++){
	  dfdr   [f1][f2] += fmatrixf[    matter][i][f1][f2]*Phi0avg   [f1][f2] - phaseVolDensity(block   [f1][f2], j);
	  dfbardr[f1][f2] += fmatrixf[antimatter][i][f1][f2]*Phi0avgbar[f1][f2] - phaseVolDensity(blockbar[f1][f2], j);
	}
    }
    
    // Make sure dfdr is Hermitian
    dfdr   [mu][e ] = conj(dfdr   [e][mu]);
    dfbardr[mu][e ] = conj(dfbardr[e][mu]);

    // update fmatrixf
    for(flavour f1=e; f1<=mu; f1++)
      for(flavour f2=e; f2<=mu; f2++){
    	fmatrixf[    matter][i][f1][f2] +=    dfdr[f1][f2] * dr;
    	fmatrixf[antimatter][i][f1][f2] += dfbardr[f1][f2] * dr;
      }

    // check that everything makes sense
    for(state s=matter; s<=antimatter; s++){
      for(flavour f1=e; f1<=mu; f1++)
	for(flavour f2=e; f2<=mu; f2++){
	  assert(fmatrixf[s][i][f1][f2] == fmatrixf[s][i][f1][f2]);
      }
      assert(abs(fmatrixf[s][i][e ][e ]) < 1.);
      assert(abs(fmatrixf[s][i][mu][mu]) < 1.);
    }
  }
}
