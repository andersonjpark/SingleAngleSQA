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

  // temporaries
  array<array<double,NF>,NE> kV;
  array<MATRIX<complex<double>,NF,NF>,NM> UV;
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> HfV;
  array<array<MATRIX<complex<double>,NF,NF>,NF>,NE> CV;
  array<array<array<double,NF>,NF>,NE> AV;
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> U0; // mixing angles to MSW basis at initial point
  array<MATRIX<complex<double>,NF,NF>,NM> VfMSW;
  array<array< array<double,NF>,NE>,NM> k0;
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> UWBW;
  array<array<array<MATRIX<complex<double>,NF,NF>,NS>,NE>,NM> Sa;
  array<array<double,NM>,NE> dphi_dr_interact, dtheta_dr_interact;
  array<array<double,NM>,NE> dphi_dr_osc,      dtheta_dr_osc;
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> pmatrixf0, pmatrixm0, Scumulative;

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
    
  /*   // MSW potential matrix */
  /*   VfMSW[matter][e][e]=Ve(rho,Ye); */
  /*   VfMSW[matter][mu][mu]=Vmu(rho,Ye); */
  /*   VfMSW[antimatter]=-Conjugate(VfMSW[matter]); */
    
  /*   // other matrices */
  /*   for(int m=matter; m<=antimatter; m++){ */
  /*     for(int i=0;i<NE;i++){ */
  /* 	MATRIX<complex<double>,NF,NF> Hf0=HfV[m][i]+ VfMSW[m]; */
  /* 	k0[m][i] = (m==matter? k(Hf0) : kbar(Hf0) ); */
  /* 	array<double,1> deltak0 = (m==matter ? deltak(Hf0) : deltakbar(Hf0) ); */
  /* 	array<MATRIX<complex<double>,NF,NF>,NF> C0 = CofactorMatrices(Hf0,k0[m][i]); */
  /* 	array<array<double,NF>,NF> A0; */
  /* 	for(int j=0;j<=NF-1;j++){ */
  /* 	  if(real(C0[j][mu][e]*CV[i][j][mu][e]) < 0.) */
  /* 	    A0[j][e]=-AV[i][j][e]; */
  /* 	  else A0[j][e]=AV[i][j][e]; */
  /* 	  A0[j][mu]=AV[i][j][mu]; */
  /* 	} */
  /* 	U0[m][i]=U(deltak0,C0,A0); */
  /* 	if(m==antimatter) U0[m][i] = Conjugate(U0[m][i]); */
  /*     } */
  /*   } */
    
    // set Scumulative to identity
    for(int m=0; m<NM; m++)
      for(int ig=0; ig<NE; ig++)
	for(int f1=0; f1<NF; f1++){
	  for(int f2=0; f2<NF; f2++)
	    Scumulative[m][ig][f1][f2] = 0.;
	  Scumulative[m][ig][f1][f1] = 1.;
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

    // potential
    for(state m=matter; m<=antimatter; m++){
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
      }
    }

  }
};



#endif
