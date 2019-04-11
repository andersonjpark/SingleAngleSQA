#ifndef STATE_H
#define STATE_H

#include "misc.h"
#include "mixing_angles.h"
#include "isospin.h"

class State{
 public:
  double r;
  double rho, T, Ye;
  double drhodr, dYedr;

  // energy grid
  array<double,NE> E;
  
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
  array<MATRIX<complex<double>,NF,NF>,NM> VfSI;
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
  }


  void update_potential(const DISCONTINUOUS& lnrho,
			const DISCONTINUOUS& temperature,
			const DISCONTINUOUS& electronfraction,
			const array<array<array<DISCONTINUOUS,NF>,NE>,NM>& P_unosc,
			const array<array<MATRIX<complex<double>,NF,NF>,NE>,NM>& HfV,
			const State& s0){
    // fluid background
    rho = exp(lnrho(r));
    T = temperature(r);
    Ye = electronfraction(r);
    drhodr=rho*lnrho.Derivative(r);
    dYedr=electronfraction.Derivative(r);

    // Matter Potential
    array<MATRIX<complex<double>,NF,NF>,NM> VfMSW;
    VfMSW[matter][e][e]=Ve(rho,Ye);
    VfMSW[matter][mu][mu]=Vmu(rho,Ye);
    VfMSW[matter][e][mu] = 0;
    VfMSW[matter][mu][e] = 0;
    VfMSW[antimatter]=-Conjugate(VfMSW[matter]);
    
    // SI potential
    VfSI[matter] = MATRIX<complex<double>,NF,NF>();
    VfSI[antimatter] = MATRIX<complex<double>,NF,NF>();
    #pragma omp parallel for collapse(2)
    for(int m=matter; m<=antimatter; m++){
      for(int i=0;i<=NE-1;i++){
	// decompose unoscillated potential
	MATRIX<complex<double>,NF,NF> pmatrixf0;
	pmatrixf0[e ][e ] = complex<double>(P_unosc[m][i][e](r),0);
	pmatrixf0[mu][e ] = complex<double>(0,0);
	pmatrixf0[e ][mu] = complex<double>(0,0);
	pmatrixf0[mu][mu] = complex<double>(P_unosc[m][i][mu](r),0);

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
    array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> old_fmatrixf = fmatrixf;
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
	double hold[4], hnew[4];
	pauli_decompose(old_fmatrixf[m][i], hold);
	pauli_decompose(    fmatrixf[m][i], hnew);
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

  void assert_noNaN(){
    for(state m=matter; m<=antimatter; m++){
      for(int i=0; i<NE; i++){

	for(flavour f1=e; f1<=mu; f1++)
	  for(flavour f2=e; f2<=mu; f2++)
	    assert(fmatrixf[m][i][f1][f2] == fmatrixf[m][i][f1][f2]);

	for(solution x=msw;x<=si;x++)
	    for(int j=0;j<=NY-1;j++)
	      assert(Y[m][i][x][j] == Y[m][i][x][j]);
      }
    }
  }
};



#endif
