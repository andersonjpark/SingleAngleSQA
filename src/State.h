#ifndef STATE_H
#define STATE_H

#include "misc.h"
#include "mixing_angles.h"
#include "isospin.h"
#include "profile.h"

class State{
 public:
  double Ecom_Elab, Elab_Elab0;
  double r;
  double rho, T, Ye;

  // energy grid
  array<double,NE> E, Etop;
  
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
  array<MATRIX<complex<double>,NF,NF>,NM> VfSI;
  array<array<array<array<double,NF>,NF>,NE>,NM> AA;
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> Sf, SThisStep, VfMSW, dVfMSWdr, UU;
  array<array<array<MATRIX<complex<double>,NF,NF>,NF>,NE>,NM> CC; 
  array<array<array<MATRIX<complex<double>,NF,NF>,NS>,NE>,NM> BB,WW;

  // other matrices
  array<array<double,NM>,NE> dphi_dr_interact, dtheta_dr_interact;
  array<array<double,NM>,NE> dphi_dr_osc,      dtheta_dr_osc;

  State(const array<double,NE>& E, const array<double,NE>& Etop){
    this->E    = E;
    this->Etop = Etop;
    
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
  }

  static double Vphase(double Elow, double Ehi){
    assert(Ehi>Elow);
    assert(Elow>=0);
    double dE3 = pow(Ehi,3) - pow(Elow,3);
    double result = 4.*M_PI * dE3/3. / pow(2.*M_PI*cgs::constants::hbarc,3);
    return result;
  }

  static double Vphase_overlap(double Elow1, double Ehi1, double Elow2, double Ehi2){
    assert(Elow1<Ehi1);
    assert(Elow2<Ehi2);
    assert(Elow1>=0);
    assert(Elow2>=0);
    double Elow = max(Elow1,Elow2);
    double Ehi  = min(Ehi1, Ehi2);
    if(Ehi<=Elow) return 0;
    else return Vphase(Elow, Ehi);
  }

  static double Vphase_overlap_comoving(int i0, const array<double,NE>& Etop0,
					int ilab, const array<double,NE>& Etop_lab,
					double Ecom_Elab){
    assert(i0>=0);
    assert(i0<NE);
    assert(ilab>=0);
    assert(ilab<NE);

    // compute top and bottom energies in the comoving frame
    double Etop_fromLab = Etop_lab[ilab] * Ecom_Elab;
    double Ebottom_fromLab = Ebottom(ilab, Etop_lab) * Ecom_Elab;

    // grab the top and bottom energies from E0
    double Etop_fromCom = Etop0[i0];
    double Ebottom_fromCom = Ebottom(i0, Etop0);

    // extend fromLab grid out to end of comoving grid
    // so S maps onto whole comoving bin
    if(ilab==NE-1) Etop_fromLab = max(Etop_fromLab, Etop_fromCom);

    // get the phase space overlap
    double V_overlap = Vphase_overlap(Ebottom_fromLab, Etop_fromLab, Ebottom_fromCom, Etop_fromCom);
    assert(V_overlap <= Vphase(i0,Etop0));
    return V_overlap;
  }

  static double Ebottom(int i, const array<double,NE>& Etop){
    assert(i>=0);
    assert(i<NE);
    double result = (i>0 ? Etop[i-1] : 0);
    assert(result < Etop[i]);
    return result;
  }
  
  static double Vphase(int i, const array<double,NE>& Etop){
    assert(i>=0);
    assert(i<NE);
    return Vphase(Ebottom(i,Etop), Etop[i]);
  }

  void update_background(const Profile& profile, const State& s0){
    rho = exp(profile.lnrho(r));
    T = profile.temperature(r);
    Ye = profile.Ye(r);
    Ecom_Elab = profile.Ecom_Elab(r);
    Elab_Elab0 = profile.Elab_Elab0(r);

    for(int i=0; i<NE; i++){
      E[i]    = s0.E[i]    * Elab_Elab0;
      Etop[i] = s0.Etop[i] * Elab_Elab0;
    }
  }
  
  void update_potential(const Profile& profile, const State& s0){
    // fluid background
    update_background(profile, s0);

    // vacuum potential
    array<array<double,NF>,NE> kV = set_kV(E);
    array<MATRIX<complex<double>,NF,NF>,NM> UV = Evaluate_UV();
    array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> VfVac = Evaluate_VfVac(kV,UV);

    // derivative of vacuum potential (which is proportional to 1/E)
    double VfVac_derivative_fac = -profile.Elab_Elab0.Derivative(r) / Elab_Elab0;

    // Matter Potential
    array<MATRIX<complex<double>,NF,NF>,NM> VfMatter, dVfMatterdr;
    double matter_potential=M_SQRT2*cgs::constants::GF/cgs::constants::Mp*rho*Ye*Ecom_Elab;
    VfMatter[matter][e ][e ] = matter_potential;
    VfMatter[matter][mu][mu] = 0;
    VfMatter[matter][e ][mu] = 0;
    VfMatter[matter][mu][e ] = 0;
    VfMatter[antimatter]=-Conjugate(VfMatter[matter]);

    // derivative of matter potential
    double dlogrhodr=profile.lnrho.Derivative(r);
    double dYedr=profile.Ye.Derivative(r);
    dVfMatterdr[matter] = VfMatter[matter] * (dlogrhodr + dYedr/Ye);
    dVfMatterdr[antimatter]=-Conjugate(dVfMatterdr[matter]);
  
    // SI potential
    VfSI[matter] = MATRIX<complex<double>,NF,NF>();
    VfSI[antimatter] = MATRIX<complex<double>,NF,NF>();
    #pragma omp parallel for collapse(2)
    for(int m=matter; m<=antimatter; m++){
      for(int i=0;i<=NE-1;i++){

	// stuff that used to be in K()
	MATRIX<complex<double>,NF,NF> dVfVacdr = VfVac[m][i] * VfVac_derivative_fac;
	VfMSW[m][i] = VfVac[m][i]+VfMatter[m];
	dVfMSWdr[m][i] = dVfMatterdr[m] + dVfVacdr;
	kk[m][i] = k(VfMSW[m][i]);
	dkk[m][i] = deltak(VfMSW[m][i]);
	CC[m][i]  = CofactorMatrices(VfMSW[m][i],kk[m][i]);
	AA[m][i] = MixingMatrixFactors(CC[m][i],s0.CC[m][i],s0.AA[m][i]);
	UU[m][i] = U(dkk[m][i],CC[m][i],AA[m][i]);
	BB[m][i][msw] = B(Y[m][i][msw]);
	BB[m][i][si ] = B(Y[m][i][si ]);
	WW[m][i][msw] = W(Y[m][i][msw]);
	WW[m][i][si ] = W(Y[m][i][si ]);
	SThisStep[m][i] = WW[m][i][msw] * BB[m][i][msw] * WW[m][i][si] * BB[m][i][si];
	Sf[m][i] = UU[m][i] * SThisStep[m][i] * Scumulative[m][i] * Adjoint(s0.UU[m][i]);

      }
    }

    // calculate the self-interaction potential
    #pragma omp parallel for collapse(2)
    for(int m=matter; m<=antimatter; m++){
      for(int i0=0; i0<NE; i0++){

	// decompose unoscillated potential
	// remains in comoving frame for the time being
	MATRIX<complex<double>,NF,NF> pmatrixf0;
	pmatrixf0[e ][e ] = sqrt(2.)*cgs::constants::GF
	  * complex<double>(profile.Dens_unosc[m][i0][e](r) - 
			    profile.Flux_unosc[m][i0][e](r),0);
	pmatrixf0[mu][e ] = complex<double>(0,0);
	pmatrixf0[e ][mu] = complex<double>(0,0);
	pmatrixf0[mu][mu] = sqrt(2.)*cgs::constants::GF
	  * complex<double>(profile.Dens_unosc[m][i0][mu](r) -
			    profile.Flux_unosc[m][i0][mu](r),0);

	MATRIX<complex<double>,NF,NF> VfSIE;
	double V0 = Vphase(i0, s0.Etop);
	double total_overlap_fraction = 0; // should end up being 1
	for(int ilab=0; ilab<NE; ilab++){

	  // calculate fraction of bin i0 that overlaps with bin ilab
	  double V_overlap = Vphase_overlap_comoving(i0, s0.Etop, ilab, Etop, Ecom_Elab);
	  double overlap_fraction = V_overlap / V0;
	  total_overlap_fraction += overlap_fraction;
	  
	  // calculate contribution to the potential due to this overlapping
	  // segment of bin i0, oscillating it with Sf from bin ilab
	  VfSIE += Sf[m][ilab]
	    * (pmatrixf0 * overlap_fraction)
	    * Adjoint(Sf[m][ilab]);
	}
	assert(fabs(total_overlap_fraction-1.)<1e-6);

	// accumulate the contribution from bin i0 onto the com-frame potential
        #pragma omp critical
	VfSI[matter] += (m==matter ? VfSIE : -Conjugate(VfSIE));
      }
    }

    // convert comoving-frame potential to lab-frame potential
    VfSI[matter] *= Ecom_Elab;
    
    // set antimatter potential
    VfSI[antimatter]=-Conjugate(VfSI[matter]);
  }


  //====================//
  // OSCILLATED MOMENTS //
  //====================//
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> oscillated_moments(const array<array<array<DISCONTINUOUS,NF>,NE>,NM>& D_unosc) const{
    
    array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> DBackground;
    for(state m=matter; m<=antimatter; m++){
      for(int ig=0; ig<NE; ig++){
	for(flavour f=e; f<=mu; f++)
	  DBackground[m][ig][f][f] = D_unosc[m][ig][f](r);
	DBackground[m][ig] = Sf[m][ig] * DBackground[m][ig] * Adjoint(Sf[m][ig]);
      }
    }
    return DBackground;
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
    array<array<double,NF>,NE> kV = set_kV(E);
    array<MATRIX<complex<double>,NF,NF>,NM> UV = Evaluate_UV();
    array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> VfVac = Evaluate_VfVac(kV,UV);
    array<array<MATRIX<complex<double>,NF,NF>,NF>,NE> CV = Evaluate_CV(kV, VfVac);
    array<array<array<double,NF>,NF>,NE> AV = Evaluate_AV(kV,VfVac,UV);
    for(state m=matter; m<=antimatter; m++){
      for(int i=0;i<=NE-1;i++){
	for(int j=0;j<=NF-1;j++){
	  if(real(CC[m][i][j][mu][e]*CV[i][j][mu][e]) < 0.)
	    AA[m][i][j][e]=-AV[i][j][e];
	  else AA[m][i][j][e]=AV[i][j][e];
	  AA[m][i][j][mu]=AV[i][j][mu];
	}
	UU[m][i]=U(dkk[m][i],CC[m][i],AA[m][i]);
      }
    }

    // determine eigenvalue ordering
    if(kV[0][1]>kV[0][0])
      cout<<"\n\nNormal hierarchy" << endl;
    else{
      if(kV[0][1]<kV[0][0])
	cout<<"\n\nInverted hierarchy" << endl;
      else{
	cout<<endl<<endl<<"Neither normal or Inverted"<<endl;
	abort();
      }
    }
  

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
	  fmatrixf[m][i][f][f] = D_unosc[m][i][f](r) / Vphase(i,Etop);
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
