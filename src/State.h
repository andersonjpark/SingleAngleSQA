#ifndef STATE_H
#define STATE_H

#include "misc.h"
#include "mixing_angles.h"
#include "isospin.h"
#include "profile.h"

class State{
 public:
  double Ecom_Elab;
  double r;
  double rho, T, Ye;

  // energy grid
  array<double,NE> E, Vphase;
  
  // distribution function in the direction of the trajectory
  // value at the last reset
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> fmatrixf;

  // Evolution variables for neutrino oscillation
  // Y Describes oscillation since Scumulative was last updated
  // Scumulative is evolution matrix from initial mass basis (s0.UU) to mass basis
  // at the time of the last update. S = WBWB(Y) * Scumulative to current mass basis.
  array<array<array<array<double,NY>,NS>,NE>,NM> Y;
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> Scumulative;
  
  // Intermediate quantities used in calculating potentials and K
  array<array<array<double,NF>,NE>,NM> kk;
  array<array<array<double,NF-1>,NE>,NM> dkk;
  array<array<array<array<double,NF>,NF>,NE>,NM> AA;
  array<MATRIX<complex<double>,NF,NF>,NM> VfSI, VfMSW, dVfMSWdr;
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> Sf, SThisStep, Hf, UU;
  array<array<array<MATRIX<complex<double>,NF,NF>,NF>,NE>,NM> CC; 
  array<array<array<MATRIX<complex<double>,NF,NF>,NS>,NE>,NM> BB,WW;

  // other matrices
  array<array<double,NM>,NE> dphi_dr_interact, dtheta_dr_interact;
  array<array<double,NM>,NE> dphi_dr_osc,      dtheta_dr_osc;

  State(const array<double,NE>& E){
    this->E = E;
    
    // set Scumulative to identity
    for(int m=0; m<NM; m++){
      for(int ig=0; ig<NE; ig++){
	for(int f1=0; f1<NF; f1++){
	  for(int f2=0; f2<NF; f2++)
	    Scumulative[m][ig][f1][f2] = 0.;
	  Scumulative[m][ig][f1][f1] = 1.;
	}
	Y[m][ig] = YIdentity;
      }
    }

    // set Vphase
    for(int i=0; i<NE; i++){
      double dlogE = (log(E[NE-1]) - log(E[0])) / (NE-1.);
      double Elow = exp(log(E[0]) + (i-0.5)*dlogE);
      double Ehi  = exp(log(E[0]) + (i+0.5)*dlogE);
      double dE3 =  pow(Ehi,3) - pow(Elow,3);
      Vphase[i] = 4.*M_PI * dE3/3. / pow(2.*M_PI*cgs::constants::hbarc,3);
    }
  }


  void update_potential(const Profile& profile,
			const array<array<MATRIX<complex<double>,NF,NF>,NE>,NM>& HfV,
			const State& s0){
    // fluid background
    rho = exp(profile.lnrho(r));
    T = profile.temperature(r);
    Ye = profile.Ye(r);
    Ecom_Elab = profile.Ecom_Elab(r);

    // Matter Potential
    double matter_potential = M_SQRT2*cgs::constants::GF/cgs::constants::Mp*rho*Ye*Ecom_Elab;
    VfMSW[matter][e ][e ] = matter_potential;
    VfMSW[matter][mu][mu] = 0;
    VfMSW[matter][e ][mu] = 0;
    VfMSW[matter][mu][e ] = 0;
    VfMSW[antimatter]=-Conjugate(VfMSW[matter]);

    double dlogrhodr=profile.lnrho.Derivative(r);
    double dYedr=profile.Ye.Derivative(r);
    dVfMSWdr[matter] = VfMSW[matter] * (dlogrhodr + dYedr/Ye);
    dVfMSWdr[antimatter]=-Conjugate(dVfMSWdr[matter]);
  
    // SI potential
    VfSI[matter] = MATRIX<complex<double>,NF,NF>();
    VfSI[antimatter] = MATRIX<complex<double>,NF,NF>();
    #pragma omp parallel for collapse(2)
    for(int m=matter; m<=antimatter; m++){
      for(int i=0;i<=NE-1;i++){
	// decompose unoscillated potential
	MATRIX<complex<double>,NF,NF> pmatrixf0;
	pmatrixf0[e ][e ] = sqrt(2.)*cgs::constants::GF
	  * complex<double>(profile.Dens_unosc[m][i][e](r) - 
			    profile.Flux_unosc[m][i][e](r),0);
	pmatrixf0[mu][e ] = complex<double>(0,0);
	pmatrixf0[e ][mu] = complex<double>(0,0);
	pmatrixf0[mu][mu] = sqrt(2.)*cgs::constants::GF
	  * complex<double>(profile.Dens_unosc[m][i][mu](r) -
			    profile.Flux_unosc[m][i][mu](r),0);

	// stuff that used to be in K()
	Hf[m][i] = HfV[m][i]+VfMSW[m];
	kk[m][i] = k(Hf[m][i]);
	dkk[m][i] = deltak(Hf[m][i]);
	CC[m][i]  = CofactorMatrices(Hf[m][i],kk[m][i]);
	AA[m][i] = MixingMatrixFactors(CC[m][i],s0.CC[m][i],s0.AA[m][i]);
	UU[m][i] = U(dkk[m][i],CC[m][i],AA[m][i]);
	BB[m][i][msw] = B(Y[m][i][msw]);
	BB[m][i][si ] = B(Y[m][i][si ]);
	WW[m][i][msw] = W(Y[m][i][msw]);
	WW[m][i][si ] = W(Y[m][i][si ]);
	SThisStep[m][i] = WW[m][i][msw] * BB[m][i][msw] * WW[m][i][si] * BB[m][i][si];
	Sf[m][i] = UU[m][i] * SThisStep[m][i] * Scumulative[m][i] * Adjoint(s0.UU[m][i]);

	// contribution to the self-interaction potential from this energy
	MATRIX<complex<double>,NF,NF> VfSIE = Sf[m][i] * pmatrixf0 * Adjoint(Sf[m][i]);
	if(m==antimatter) VfSIE = -Conjugate(VfSIE);
        #pragma omp critical
	VfSI[matter] += VfSIE;
      }
    }
    VfSI[antimatter]=-Conjugate(VfSI[matter]);
  }

  void accumulate_S(double dr, const State& sReset){
    #pragma omp parallel for collapse(2)
    for(int m=matter;m<=antimatter;m++){
      for(int i=0;i<=NE-1;i++){
	Scumulative[m][i] = SThisStep[m][i] * Scumulative[m][i];
	
	// convert fmatrix from flavor basis to (reset-point) mass basis
	// evolve fmatrix from reset-point to current-point mass basis
	// convert fmatrix from (current-point) mass basis to flavor basis
	MATRIX<complex<double>,NF,NF> SfThisStep =
	  UU[m][i]
	  * SThisStep[m][i]
	  * Adjoint(sReset.UU[m][i]);
	fmatrixf[m][i] = SfThisStep * fmatrixf[m][i] * Adjoint(SfThisStep);
	    
	// reset the evolution matrix to identity
	Y[m][i] = YIdentity;

	// get rate of change of fmatrix from oscillation
	array<double,4> hold = pauli_decompose(sReset.fmatrixf[m][i]);
	array<double,4> hnew = pauli_decompose(       fmatrixf[m][i]);
	double oldmag   = sqrt(hold[0]*hold[0] + hold[1]*hold[1] + hold[2]*hold[2]);
	double newmag   = sqrt(hnew[0]*hnew[0] + hnew[1]*hnew[1] + hnew[2]*hnew[2]);
	double costheta = (hold[0]*hnew[0] + hold[1]*hnew[1] + hold[2]*hnew[2]) / (newmag*oldmag);
	assert(costheta-1. < 1e-10);
	costheta = min(1.,costheta);
	dtheta_dr_osc[i][m] = (acos(hnew[2]/newmag) - acos(hold[2]/oldmag)) / dr;
	dphi_dr_osc[i][m] = (atan2(hnew[1],hnew[0]) - atan2(hold[1],hold[0])) / dr;
      }
    }
  }

  void assert_noNaN(double accuracy){
    for(state m=matter; m<=antimatter; m++){
      for(int i=0; i<NE; i++){
	
	for(flavour f1=e; f1<=mu; f1++){
	  assert(real(fmatrixf[m][i][f1][f1]) <= 1.);
	  assert(real(fmatrixf[m][i][f1][f1]) >= 0.);
	  assert(abs(imag(fmatrixf[m][i][f1][f1])) < accuracy);
	  for(flavour f2=e; f2<=mu; f2++)
	    assert(fmatrixf[m][i][f1][f2] == fmatrixf[m][i][f1][f2]);
	}

	array<double,4> isospin = pauli_decompose(fmatrixf[m][i]);
	double fperp2 = isospin[0]*isospin[0] + isospin[1]*isospin[1];
	assert(fabs(isospin[0]) <= isospin[3]);
	assert(fperp2 <= pow(min(isospin[3], 1.-isospin[3]),2) - isospin[2]*isospin[2]);

	for(solution x=msw;x<=si;x++)
	    for(int j=0;j<=NY-1;j++)
	      assert(Y[m][i][x][j] == Y[m][i][x][j]);
      }
    }
  }


  //============//
  // Initialize //
  //============//
  void initialize(const array<array<array<DISCONTINUOUS,NF>,NE>,NM>& D_unosc){
    // T should be MeV
    cout << "Setting initial data." << endl;
    cout << "rho = " << rho << " g/ccm" << endl;
    cout << "T = " << T << " MeV" << endl;
    cout << "Ye = " << Ye << endl;
    //Ye = max(Ye,__nulibtable_MOD_nulibtable_ye_min);
    /* nulibtable_range_species_range_energy_(&rho, &T, &Ye, &eas.eas.front(), */
    /* 					 &__nulibtable_MOD_nulibtable_number_species, */
    /* 					 &__nulibtable_MOD_nulibtable_number_groups, */
    /* 					 &__nulibtable_MOD_nulibtable_number_easvariables); */
  
    for(int i=0; i<NE; i++){
      for(state m=matter; m<=antimatter; m++){
	fmatrixf[m][i] = MATRIX<complex<double>,NF,NF>();
	for(flavour f=e; f<=mu; f++)
	  fmatrixf[m][i][f][f] = D_unosc[m][i][f](r) / Vphase[i];
      }
    
      cout << "GROUP " << i << endl;
      cout << "\tf = {" << real(fmatrixf[matter][i][e][e]) << ", " << real(fmatrixf[antimatter][i][e][e]) << ", " << real(fmatrixf[matter][i][mu][mu]) <<"}" << endl;
      /* cout << "\teas.emis = {" << eas.emis(0,i) << ", " << eas.emis(1,i) << ", " << eas.emis(2,i) << "}" << endl; */
      /* cout << "\teas.abs = {" << eas.abs(0,i) << ", " << eas.abs(1,i) << ", " << eas.abs(2,i) << "}" << endl; */
      /* cout << "\tBB = {" << eas.emis(0,i)/eas.abs(0,i) << ", " << eas.emis(1,i)/eas.abs(1,i) << ", " << eas.emis(2,i)/eas.abs(2,i) << "}" << endl; */
      /* cout << "\teas.scat = {" << eas.scat(0,i) << ", " << eas.scat(1,i) << ", " << eas.scat(2,i) << "}" << endl; */
    }
  }


};



#endif
