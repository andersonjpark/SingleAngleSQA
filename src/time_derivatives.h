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


//===========//
// Kinteract //
//===========//
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
      //const state mbar = (m==matter ? antimatter : matter);
      
      // intermediate variables
      complex<double> unblock_in, unblock_out;
      MATRIX<complex<double>,NF,NF> block, Pi_plus, Pi_minus;
      MATRIX<double,NF,NF> Phi0, Phi0avg,  Phi0tilde;      

      // emission
      dfdr[m][i][e ][e ] += eas.emis(se,i);
      dfdr[m][i][mu][mu] += eas.emis(sx,i);
      
      // absorption kappa_abs is <kappa> for absorption
      Phi0avg = eas.avg_matrix(eas.abs(se,i), eas.abs(sx,i));
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
	    dfdr[m][i][f1][f2] += (DBackground[m][j][f1][f2]*Phi0[f1][f2] - block[f1][f2]) / s.Vphase(i);

	// out-scattering from i to j. for blocking, get phase space vol from D[j] in j
	Phi0avg      = eas.avg_matrix(  eas.Phi0(0,i,j), eas.Phi0(2,i,j));
	Phi0tilde    = eas.tilde_matrix(eas.Phi0(0,i,j), eas.Phi0(2,i,j));
	Phi0    = Phi0avg    - Phi0tilde;
	block    = eas.blocking_term0(Phi0, s.fmatrixf[m][i], DBackground[m][j] );
	for(flavour f1=e; f1<=mu; f1++)
	  for(flavour f2=e; f2<=mu; f2++)
	    dfdr[m][i][f1][f2] += s.fmatrixf[m][i][f1][f2]*Phi0avg[f1][f2] - block[f1][f2]/s.Vphase(j);

	// Make sure dfdr is Hermitian
	dfdr[m][i][mu][e ] = conj(dfdr[m][i][e][mu]);
	
      } // other group
    } // group
  } // state
  
  return dfdr;
}



#endif
