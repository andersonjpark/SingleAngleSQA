#ifndef STATE_H
#define STATE_H

#include "misc.h"
#include "mixing_angles.h"

class State{
 public:
  double rho, T, Ye;
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
  }

};



#endif
