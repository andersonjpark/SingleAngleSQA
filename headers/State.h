#ifndef STATE_H
#define STATE_H

#include "misc.h"
#include "mixing_angles.h"

class State{
 public:
  double rho, T, Ye;
  double drhodr, dYedr;
  double r, dr_block, dr_osc, dr_int;
  int counter;
  //EAS eas;
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> fmatrixf, fmatrixm;
  hid_t dset_f, dset_r, dset_dr_osc, dset_dr_int, dset_dr_block;
  array<array<array<array<double,NY>,NS>,NE>,NM> Y;
  
  // vacuum matrices set at initial conditions
  array<array<double,NF>,NE> kV;
  array<MATRIX<complex<double>,NF,NF>,NM> UV;
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> HfV, UWBW;
  array<array<MATRIX<complex<double>,NF,NF>,NF>,NE> CV;
  array<array<array<double,NF>,NF>,NE> AV;

  // mixing angles to MSW basis at initial point
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> U0;

  // other matrices...
  array<array<array<MATRIX<complex<double>,NF,NF>,NF>,NE>,NM> C;
  array<array<array<array<double,NF>,NF>,NE>,NM> A;
  
  // potentials and potential derivatives
  array<MATRIX<complex<double>,NF,NF>,NM> VfMSW, dVfMSWdr, VfSI;
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> pmatrixf0, pmatrixm0, Scumulative;

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
  
  State(/*string nulibfilename, string eosfilename, double rho_in, double Ye_in, double T_in, double dr0, double mixing, bool do_interact*/const vector<double>& E){
  /*   r=0; */
  /*   rho = rho_in; */
  /*   T = T_in; */
  /*   Ye = Ye_in; */
  /*   eas = EAS(nulibfilename, eosfilename); */
  /*   initialize(fmatrixf,eas,rho,T,Ye, mixing, do_interact); */
  /*   dr_block = dr0; */
  /*   dr_osc = dr0; */
  /*   dr_int = dr0; */
  /*   counter = 0; */
  
    // vectors of energies and vacuum eigenvalues
    kV  = set_kV(E);
    UV  = Evaluate_UV();
    HfV = Evaluate_HfV(kV,UV);
    CV  = Evaluate_CV(kV, HfV);
    AV  = Evaluate_AV(kV,HfV,UV);
    
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


  void update_background(const DISCONTINUOUS& lnrho,
			 const DISCONTINUOUS& temperature,
			 const DISCONTINUOUS& electronfraction,
			 const array<DISCONTINUOUS,NE>& eD,
			 const array<DISCONTINUOUS,NE>& eBarD,
			 const array<DISCONTINUOUS,NE>& xD,
			 const array<DISCONTINUOUS,NE>& eP,
			 const array<DISCONTINUOUS,NE>& eBarP,
			 const array<DISCONTINUOUS,NE>& xP){
    // fluid background
    rho = exp(lnrho(r));
    T = temperature(r);
    Ye = electronfraction(r);
    drhodr=rho*lnrho.Derivative(r);
    dYedr=electronfraction.Derivative(r);

    // Matter Potential
    VfMSW[matter][e][e]=Ve(rho,Ye);
    VfMSW[matter][mu][mu]=Vmu(rho,Ye);
    VfMSW[matter][e][mu] = 0;
    VfMSW[matter][mu][e] = 0;
    VfMSW[antimatter]=-Conjugate(VfMSW[matter]);
    
    dVfMSWdr[matter][e][e]=dVedr(rho,drhodr,Ye,dYedr);
    dVfMSWdr[matter][mu][mu]=dVmudr(rho,drhodr,Ye,dYedr);
    dVfMSWdr[matter][e][mu]=0;
    dVfMSWdr[matter][mu][e]=0;
    dVfMSWdr[antimatter]=-Conjugate(dVfMSWdr[matter]);

    // SI potential
    VfSI[matter] = MATRIX<complex<double>,NF,NF>();
    VfSI[antimatter] = MATRIX<complex<double>,NF,NF>();
    #pragma omp parallel for collapse(2)
    for(int m=matter; m<=antimatter; m++){
      for(int i=0;i<=NE-1;i++){
	// decompose unoscillated potential
	double P0 = (m==matter ? eP[i](r) : eBarP[i](r));
	double P1 = xP[i](r);
	pmatrixf0[m][i][e ][e ] = complex<double>(P0,0);
	pmatrixf0[m][i][mu][e ] = complex<double>(0,0);
	pmatrixf0[m][i][e ][mu] = complex<double>(0,0);
	pmatrixf0[m][i][mu][mu] = complex<double>(P1,0);
	
	// oscillate the potential and put into the mass basis
	pmatrixm0[m][i] = Scumulative[m][i]
	  * Adjoint(U0[m][i])
	  * pmatrixf0[m][i]
	  * U0[m][i]
	  * Adjoint(Scumulative[m][i]);
	pmatrixf0[m][i] =  U0[m][i] * pmatrixm0[m][i] * Adjoint(U0[m][i]);

	// stuff that used to be in K()
	Hf[m][i] = HfV[m][i]+VfMSW[m];
	kk[m][i] = k(Hf[m][i]);
	dkk[m][i] = deltak(Hf[m][i]);
	CC[m][i]  = CofactorMatrices(Hf[m][i],kk[m][i]);
	AA[m][i] = MixingMatrixFactors(CC[m][i],C[m][i],A[m][i]);
	UU[m][i] = U(dkk[m][i],CC[m][i],AA[m][i]);
	BB[m][i] = B(Y[m][i][msw]);
	UWBW[m][i] = UU[m][i] * W(Y[m][i][msw]) * BB[m][i] * W(Y[m][i][si]);
	Sa[m][i][si] = B(Y[m][i][si]);

	// contribution to the self-interaction potential from this energy
	MATRIX<complex<double>,NF,NF> Sfm    = UWBW[m][i]*Sa[m][i][si];
	MATRIX<complex<double>,NF,NF> VfSIE = Sfm * pmatrixm0[m][i] * Adjoint(Sfm);
	if(m==antimatter) VfSIE = -Conjugate(VfSIE);
        #pragma omp critical
	VfSI[matter] += VfSIE;
      }
    }
    VfSI[antimatter]=-Conjugate(VfSI[matter]);


  }
};



#endif
