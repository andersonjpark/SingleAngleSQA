#ifndef STATE_H
#define STATE_H

#include "misc.h"
#include "mixing_angles.h"

class State{
 public:
  double r;
  double rho, T, Ye;
  double drhodr, dYedr;
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> fmatrixf;
  array<array<array<array<double,NY>,NS>,NE>,NM> Y;
  
  // mixing angles to MSW basis at initial point
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> U0;
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> UWBW;
  
  // potentials and potential derivatives
  array<MATRIX<complex<double>,NF,NF>,NM> VfSI;

  // evolution matrices
  // for fm0 in initial matter basis and fm in current matter basis,
  // fm = Scumulative fm0 Scumulative^+
  // Sf goes from flavor basis to flavor basis
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> Scumulative;
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> Sf, SThisStep;

  // other matrices
  array<array<double,NM>,NE> dphi_dr_interact, dtheta_dr_interact;
  array<array<double,NM>,NE> dphi_dr_osc,      dtheta_dr_osc;

  // stuff that used to be in K()
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> Hf;
  array<array<array<double,NF>,NE>,NM> kk;
  array<array<array<double,NF-1>,NE>,NM> dkk;
  array<array<array<MATRIX<complex<double>,NF,NF>,NF>,NE>,NM> CC; 
  array<array<array<array<double,NF>,NF>,NE>,NM> AA;
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> UU,BB;
  array<array<array<MATRIX<complex<double>,NF,NF>,NS>,NE>,NM> Sa;
  
  State(const vector<double>& E){
    
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
	BB[m][i] = B(Y[m][i][msw]);
	UWBW[m][i] = UU[m][i] * W(Y[m][i][msw]) * BB[m][i] * W(Y[m][i][si]);
	Sa[m][i][si] = B(Y[m][i][si]);
	SThisStep[m][i] = W(Y[m][i][msw]) * BB[m][i] * W(Y[m][i][si]) * Sa[m][i][si];
	Sf[m][i] = UWBW[m][i] * Sa[m][i][si] * Scumulative[m][i] * Adjoint(U0[m][i]);

	// contribution to the self-interaction potential from this energy
	MATRIX<complex<double>,NF,NF> VfSIE = Sf[m][i] * pmatrixf0 * Adjoint(Sf[m][i]);
	if(m==antimatter) VfSIE = -Conjugate(VfSIE);
        #pragma omp critical
	VfSI[matter] += VfSIE;
      }
    }
    VfSI[antimatter]=-Conjugate(VfSI[matter]);


  }
};



#endif
