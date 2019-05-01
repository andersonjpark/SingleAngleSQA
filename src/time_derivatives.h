#ifndef TIMEDERIVATIVES_H
#define TIMEDERIVATIVES_H

#include "nulib_interface.h"

//============//
// Koscillate //
//============//
array<array<array<array<double,NY>,NS>,NE>,NM> Koscillate(const State& s){

  array<array<array<array<double,NY>,NS>,NE>,NM> K;

#pragma omp parallel for collapse(2)
  for(int m=matter; m<=antimatter; m++){
    for(int i=0;i<=NE-1;i++){

      array<double,NF-1> phase;
      MATRIX<complex<double>,NF,NF> Ha,HB;
      phase[0] = M_2PI*(s.Y[m][i][msw][4]-s.Y[m][i][msw][5]);
      Ha[0][1]=0.;
      for(int j=0;j<=NF-2;j++)
	for(int k=j+1;k<=NF-1;k++)
	  for(flavour f=e;f<=mu;f++)
	    Ha[j][k]+= conj(s.UU[m][i][f][j])*s.dVfMSWdr[m][i][f][f]*s.UU[m][i][f][k];
    
      Ha[0][1] *= I*cgs::constants::hbarc/s.dkk[m][i][0]*exp(I*phase[0]);
      Ha[1][0] = conj(Ha[0][1]);
    
      // HB = -I/cgs::constants::hbarc*Ha*BB;
      HB[0][0]=-I/cgs::constants::hbarc*( Ha[0][1]*s.BB[m][i][msw][1][0] );
      HB[0][1]=-I/cgs::constants::hbarc*( Ha[0][1]*s.BB[m][i][msw][1][1] );

      array<double,4> dvdr;
      dvdr[0]=real(HB[0][1]);
      dvdr[1]=imag(HB[0][1]);
      dvdr[2]=real(HB[0][0]);
      dvdr[3]=imag(HB[0][0]);

      MATRIX<double,3,4> JI = JInverse(s.Y[m][i][msw]);
      
      array<double,NF> dkkdr = dkdr(s.UU[m][i],s.dVfMSWdr[m][i]);
      array<MATRIX<complex<double>,NF,NF>,NF> dCCdr = CofactorMatricesDerivatives(s.dVfMSWdr[m][i],dkkdr);
      array<double,NF> QQ =  Q(s.UU[m][i],s.dkk[m][i],s.CC[m][i],dCCdr);

      for(int j=0;j<=2;j++){
	K[m][i][msw][j]=0.;
	for(int k=j;k<=3;k++)
	  K[m][i][msw][j] += JI[j][k]*dvdr[k];
      }
      K[m][i][msw][3] = 0.;
      K[m][i][msw][4] = (s.kk[m][i][0]+QQ[0])/M_2PI/cgs::constants::hbarc;
      K[m][i][msw][5] = (s.kk[m][i][1]+QQ[1])/M_2PI/cgs::constants::hbarc;

      
      // *********************
      // SI part of solution *
      // *********************
      MATRIX<complex<double>,NF,NF> UWBW =
	s.UU[m][i]
	* s.WW[m][i][msw]
	* s.BB[m][i][msw]
	* s.WW[m][i][si];
      Ha = Adjoint(UWBW)*s.VfSI[m]*UWBW;

      K[m][i][si][4]=real(Ha[0][0])/(M_2PI*cgs::constants::hbarc);
      K[m][i][si][5]=real(Ha[1][1])/(M_2PI*cgs::constants::hbarc);
    
      HB[0][0]=-I/cgs::constants::hbarc*( Ha[0][1]*s.BB[m][i][si][1][0] );
      HB[0][1]=-I/cgs::constants::hbarc*( Ha[0][1]*s.BB[m][i][si][1][1] );
    
      //array<double,4> dvdr;
      dvdr[0]=real(HB[0][1]);
      dvdr[1]=imag(HB[0][1]);
      dvdr[2]=real(HB[0][0]);
      dvdr[3]=imag(HB[0][0]);
    
      //MATRIX<double,3,4>
      JI = JInverse(s.Y[m][i][si]);
    
      for(int j=0;j<=2;j++){
	K[m][i][si][j]=0.;
	for(int k=j;k<=3;k++) K[m][i][si][j]+=JI[j][k]*dvdr[k];
      }
    
      K[m][i][si][3]=0.;

      for(solution x=msw;x<=si;x++)
	for(int j=0;j<=NY-1;j++)
	  assert(K[m][i][x][j] == K[m][i][x][j]);
    }
  }
  
  return K;
}// end of K function


//==================//
// HELPER FUNCTIONS //
//==================//
MATRIX<double,2,2> avg_matrix(const double eval, const double muval){
  MATRIX<double,2,2> result;
  result[e ][e ] = eval;
  result[mu][mu] = muval;
  result[e ][mu] = result[mu][e] = (eval + muval) / 2.;
  return result;
}
MATRIX<double,2,2> tilde_matrix(const double eval, const double muval){
  MATRIX<double,2,2> result;
  result[e ][e ] = 0;
  result[mu][mu] = 0;
  result[e ][mu] = result[mu][e] = (eval - muval) / (4.*sin2thetaW);
  return result;
}
MATRIX<complex<double>,2,2> blocking_term0(const MATRIX<double,2,2>& Phi0matrix,
					   const MATRIX<complex<double>,2,2>& f,
					   const MATRIX<complex<double>,2,2>& fp){
  MATRIX<complex<double>,2,2> result;
  for(flavour fa=e; fa<=mu; fa++)
    for(flavour fb=e; fb<=mu; fb++){
      result[fa][fb] = 0;
      for(flavour fc=e; fc<=mu; fc++){
	result[fa][fb] += 0.25 * (Phi0matrix[fc][fb]*f[fa][fc]*fp[fc][fb] + Phi0matrix[fa][fc]*fp[fa][fc]*f[fc][fb]);
      }
    }
  return result;
}


//===========//
// Kinteract //
//===========//
array<array<MATRIX<complex<double>,NF,NF>,NE>,NM>
  Kinteract(const State& s, const State& s0, const Profile& profile){

  // set up the array to be returned
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> dfdr;

  // get oscillated background density
  array<array<array<MATRIX<complex<double>,NF,NF>,NMOMENTS>,NE>,NM> MBackground = s.oscillated_moments(profile, s0);

#pragma omp parallel for collapse(2)
  for(int m=matter; m<=antimatter; m++){
    for(int i=0; i<NE; i++){

      // get nulib species indices
      const int se = (m==matter ? 0 : 1);
      const int sx = (m==matter ? 2 : 3);
      //const state mbar = (m==matter ? antimatter : matter);
      
      // intermediate variables
      complex<double> unblock_in, unblock_out;
      MATRIX<complex<double>,NF,NF> block, Pi_plus, Pi_minus;

      // absorption and emission
      MATRIX<double,NF,NF> kappa_abs_avg = avg_matrix(eas.abs(se,s.E[i]), eas.abs(sx,s.E[i]));
      double E_kT = s.E[i]/(1e6*eV) / s.T;
      dfdr[m][i][e ][e ] += kappa_abs_avg[e ][e ] * eas.fermidirac(se, E_kT);
      dfdr[m][i][mu][mu] += kappa_abs_avg[mu][mu] * eas.fermidirac(sx, E_kT);
      for(flavour f1=e; f1<=mu; f1++)
	for(flavour f2=e; f2<=mu; f2++)
	  dfdr[m][i][f1][f2] -= kappa_abs_avg[f1][f2] * s.fmatrixf[m][i][f1][f2];
    
      // scattering and pair annihilation
      // no factor of 1/2 in front of Phi0 because
      // integrated over outgoing theta, assumed isotropic
      double Vi = Vphase(i, s.Etop);
      for(int j=0; j<NE; j++){
	double Vj = Vphase(j, s.Etop);
	array<MATRIX<double,NF,NF>,2> Phi, PhiAvg,  PhiTilde;      
	array<double,2> Phi_scat_e, Phi_scat_x;
	
	// in-scattering from j to i. D*Phi0 give density scattered
	// divide by phase space vol in i to get f
	Phi_scat_e = eas.Phi_scat(se, s.E[i], s.E[j], Vi, Vj);
	Phi_scat_x = eas.Phi_scat(sx, s.E[i], s.E[j], Vi, Vj);
	for(int mom=0; mom<1; mom++){
	  PhiAvg[mom]   = avg_matrix(  Phi_scat_e[mom], Phi_scat_x[mom]);
	  PhiTilde[mom] = tilde_matrix(Phi_scat_e[mom], Phi_scat_x[mom]);
	  Phi[mom] = PhiAvg[mom] - PhiTilde[mom];
	}
	block    = blocking_term0(Phi[0], s.fmatrixf[m][i], MBackground[m][j][0]);
	for(flavour f1=e; f1<=mu; f1++)
	  for(flavour f2=e; f2<=mu; f2++)
	    dfdr[m][i][f1][f2] += (MBackground[m][j][0][f1][f2]*Phi[0][f1][f2] - block[f1][f2]) / Vi;

	//out-scattering from i to j. for blocking, get phase space vol from D[j] in j
	Phi_scat_e = eas.Phi_scat(se, s.E[i], s.E[j], Vi, Vj);
	Phi_scat_x = eas.Phi_scat(sx, s.E[i], s.E[j], Vi, Vj);
	for(int mom=0; mom<1; mom++){
	  PhiAvg[mom]   = avg_matrix(  Phi_scat_e[mom], Phi_scat_x[mom]);
	  PhiTilde[mom] = tilde_matrix(Phi_scat_e[mom], Phi_scat_x[mom]);
	  Phi[mom] = PhiAvg[mom] - PhiTilde[mom];
	}
	block    = blocking_term0(Phi[0], s.fmatrixf[m][i], MBackground[m][j][0] );
	for(flavour f1=e; f1<=mu; f1++)
	  for(flavour f2=e; f2<=mu; f2++)
	    dfdr[m][i][f1][f2] += s.fmatrixf[m][i][f1][f2]*PhiAvg[0][f1][f2] - block[f1][f2]/Vj;
	
	
      } // other group

      // Make sure dfdr is Hermitian
      dfdr[m][i][mu][e ] = conj(dfdr[m][i][e][mu]);
    } // group
  } // state
  
  return dfdr;
}



#endif
