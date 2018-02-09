#include <cmath>
#include <complex>
using std::complex;
using std::polar;
using std::abs;
using std::arg;
using std::real;
using std::imag;
#include <cstdarg>
using std::va_list;
#include <cstdlib>

#include<iostream>
using std::cout;
#include<ostream>
using std::ostream;
using std::endl;
using std::flush;
#include<fstream>
using std::ifstream;
using std::ofstream;
#include<sstream>
using std::stringstream;

#include<algorithm>
using std::min;
using std::max;
using std::swap;
using std::lower_bound;
#include<string>
using std::string;
#include <utility>
using std::pair;

#include<functional>
#include<limits>
using std::numeric_limits;

#include<vector>
using std::vector;

#include "mstl.h"
using namespace physics;
using namespace prefixes;
using interpolation::DISCONTINUOUS;

// global variables
DISCONTINUOUS rho, lnrho, Ye, temperature; // rho is the mass density
double NSI;
int NEP(8);

// headers
#include "headers/parameters.h"
#include "headers/potentials.h"
#include "headers/single angle.h"
#include "headers/flavour basis.h"
#include "headers/eigenvalues.h"
#include "headers/mixing angles.h"
#include "headers/adiabatic basis.h"
#include "headers/flux.h"
#include "headers/jacobians.h"
#include "headers/multiEnergy.h"
#include "headers/MNR.h"
#include "headers/nulib_interface.h"

//vector<vector<MATRIX<complex<double>,NF,NF> > > rhomatrixf0(NM), rhomatrixm0(NM);
vector<vector<MATRIX<complex<double>,NF,NF> > > pmatrixf0(NM), pmatrixm0(NM);
vector<vector<MATRIX<complex<double>,NF,NF> > > fmatrixf(NM), fmatrixm(NM);
vector<DISCONTINUOUS> eP,eBarP,xP;

//===//
// B //
//===//
MATRIX<complex<double>,NF,NF> B(vector<double> y){
  MATRIX<complex<double>,NF,NF> s;
  double cPsi1=cos(y[0]),sPsi1=sin(y[0]), cPsi2=cos(y[1]),sPsi2=sin(y[1]), cPsi3=cos(y[2]),sPsi3=sin(y[2]);
  
  s[0][1] = cPsi1 + I*sPsi1*cPsi2;
  sPsi1 *= sPsi2;
  s[0][0] = sPsi1 * (cPsi3 + I*sPsi3);

  s[1][0] = -y[3]*conj(s[0][1]);
  s[1][1] =  y[3]*conj(s[0][0]);

  return s;
}


//======//
// getP //
//======//
double MU(const double r, const double E){ // erg
  double dV = deltaV(E);
  return 1e5 * dV*cgs::constants::hbarc*2.*M_PI * exp(-r*dV / 10.)  / (double)NE;
}
double ALPHA(const double r){
  return 4./3.;
}
void getP(const double r)
{
  for(int i=0;i<=NE-1;i++){
    for(flavour f=e;f<=mu;f++) for(flavour fp=e; fp<=mu; fp++){
	pmatrixf0[matter    ][i][f][fp] = 0;
	pmatrixf0[antimatter][i][f][fp] = 0;
      }
    pmatrixf0[matter    ][i][e ][e ]=MU(r, E[0]);           //UNCOMMENT eP[i](r); //
    pmatrixf0[antimatter][i][e ][e ]=MU(r, E[0]) * ALPHA(r);//UNCOMMENT eBarP[i](r); //
    pmatrixf0[antimatter][i][mu][mu]=0;//UNCOMMENT xP[i](r);
    pmatrixf0[matter    ][i][mu][mu]=0;//UNCOMMENT xP[i](r);
  }
}


//===//
// K //
//===//
void K(double r,
       double dr,
       vector<vector<vector<vector<double> > > > &Y,
       vector<vector<vector<MATRIX<complex<double>,NF,NF> > > > &C0,
       vector<vector<vector<vector<double> > > > &A0,
       vector<vector<vector<vector<double> > > > &K){

  MATRIX<complex<double>,NF,NF> VfSI,VfSIbar;  // self-interaction potential
  vector<MATRIX<complex<double>,NF,NF> > VfSIE(NE); // contribution to self-interaction potential from each energy
  MATRIX<complex<double>,NF,NF> VfMSW,VfMSWbar;
  MATRIX<complex<double>,NF,NF> dVfMSWdr,dVfMSWbardr;
  MATRIX<complex<double>,NF,NF> Hf,Hfbar, UU,UUbar;
  vector<double> kk,kkbar, dkk,dkkbar, dkkdr,dkkbardr, QQ,QQbar;
  vector<MATRIX<complex<double>,NF,NF> > CC,dCCdr;
  vector<vector<double> > AA;
  MATRIX<complex<double>,NF,NF> BB,BBbar;
  MATRIX<complex<double>,NF,NF> Sfm,Sfmbar;
  vector<vector<MATRIX<complex<double>,NF,NF> > > 
    Sa(NE,vector<MATRIX<complex<double>,NF,NF> >(NS)), 
    Sabar(NE,vector<MATRIX<complex<double>,NF,NF> >(NS));
  vector<MATRIX<complex<double>,NF,NF> > UWBW(NE);
  vector<MATRIX<complex<double>,NF,NF> > UWBWbar(NE);
  double rrho,drrhodr, YYe,dYYedr;
  MATRIX<double,3,4> JI;
  int i;
  MATRIX<complex<double>,NF,NF> Ha;
  MATRIX<complex<double>,NF,NF> HB;
  vector<double> phase(1);
  vector<double> dvdr(4);
  // *************
  rrho=-1;//UNCOMMENT exp(lnrho(log(r)));
  drrhodr=-1;//UNCOMMENT rrho*lnrho.Derivative(log(r))/r;
  YYe=-1; //UNCOMMENT Ye(r);
  dYYedr=-1; //UNCOMMENT Ye.Derivative(r);
  VfMSW[e][e]=Ve(rrho,YYe);
  VfMSW[mu][mu]=Vmu(rrho,YYe);
  VfMSWbar=-Conjugate(VfMSW);
  dVfMSWdr[e][e]=dVedr(rrho,drrhodr,YYe,dYYedr);
  dVfMSWdr[mu][mu]=dVmudr(rrho,drrhodr,YYe,dYYedr);
  dVfMSWbardr=-Conjugate(dVfMSWdr);

#pragma omp parallel for schedule(auto) private(Hf,Hfbar,UU,UUbar,kk,kkbar,dkk,dkkbar,dkkdr,dkkbardr,QQ,QQbar,AA,CC,dCCdr,BB,BBbar,Sfm,Sfmbar,JI) firstprivate(Ha,HB,dvdr,phase)
  for(i=0;i<=NE-1;i++){
    Hf  = HfV[matter][i]+VfMSW;
    kk  = k(Hf);
    dkk = deltak(Hf);
    CC  = CofactorMatrices(Hf,kk);
    AA  = MixingMatrixFactors(CC,C0[matter][i],A0[matter][i]);
    UU  = U(dkk,CC,AA);
    BB  = B(Y[matter][i][msw]);
    Sa[i][si] = B(Y[matter][i][si]);
    UWBW[i] = UU * W(Y[matter][i][msw]) * BB * W(Y[matter][i][si]);
    
    Hfbar = HfV[antimatter][i] + VfMSWbar;
    kkbar = kbar(Hfbar);
    dkkbar = deltakbar(Hfbar);
    CC = CofactorMatrices(Hfbar,kkbar);
    AA = MixingMatrixFactors(CC,C0[antimatter][i],A0[antimatter][i]);
    UUbar = Conjugate(U(dkkbar,CC,AA));
    BBbar = B(Y[antimatter][i][msw]);
    Sabar[i][si] = B(Y[antimatter][i][si]);
    UWBWbar[i] = UUbar * W(Y[antimatter][i][msw]) *BBbar * W(Y[antimatter][i][si]);
    
    // ****************
    // Matter section *
    // ****************
    phase[0] = M_2PI*(Y[matter][i][msw][4]-Y[matter][i][msw][5]);
    Ha[0][1]=0.;
    for(int j=0;j<=NF-2;j++)
      for(int k=j+1;k<=NF-1;k++)
	for(flavour f=e;f<=mu;f++)
	  Ha[j][k]+= conj(UU[f][j])*dVfMSWdr[f][f]*UU[f][k];
    
    Ha[0][1] *= I*cgs::constants::hbarc/dkk[0]*exp(I*phase[0]);
    Ha[1][0] = conj(Ha[0][1]);
    
    // HB = -I/cgs::constants::hbarc*Ha*BB;
    HB[0][0]=-I/cgs::constants::hbarc*( Ha[0][1]*BB[1][0] );
    HB[0][1]=-I/cgs::constants::hbarc*( Ha[0][1]*BB[1][1] );
    
    dvdr[0]=real(HB[0][1]);
    dvdr[1]=imag(HB[0][1]);
    dvdr[2]=real(HB[0][0]);
    dvdr[3]=imag(HB[0][0]);
    
    JI = JInverse(Y[matter][i][msw]);
    
    for(int j=0;j<=2;j++){
      K[matter][i][msw][j]=0.;
      for(int k=j;k<=3;k++) K[matter][i][msw][j] += JI[j][k]*dvdr[k];
      K[matter][i][msw][j]*=dr;
    }
    
    K[matter][i][msw][3] = 0.;
    dkkdr = dkdr(UU,dVfMSWdr);
    dCCdr = CofactorMatricesDerivatives(Hf,dVfMSWdr,dkk,dkkdr);
    QQ = Q(UU,dkk,CC,dCCdr);
    
    K[matter][i][msw][4] = (kk[0]+QQ[0])*dr/M_2PI/cgs::constants::hbarc;
    K[matter][i][msw][5] = (kk[1]+QQ[1])*dr/M_2PI/cgs::constants::hbarc;

    // ********************
    // Antimatter section *
    // ********************
    phase[0] = M_2PI*(Y[antimatter][i][msw][4]-Y[antimatter][i][msw][5]);
    Ha[0][1] = 0.;
    for(int j=0;j<=NF-2;j++)
      for(int k=j+1;k<=NF-1;k++)
	for(flavour f=e;f<=mu;f++)
	  Ha[j][k]+=conj(UUbar[f][j])*dVfMSWbardr[f][f]*UUbar[f][k];
    
    Ha[0][1] *= I*cgs::constants::hbarc/dkkbar[0]*exp(I*phase[0]);
    Ha[1][0] = conj(Ha[0][1]);
    
    //HB=-I/cgs::constants::hbarc*Ha*BBbar;
    HB[0][0]=-I/cgs::constants::hbarc*( Ha[0][1]*BBbar[1][0] );
    HB[0][1]=-I/cgs::constants::hbarc*( Ha[0][1]*BBbar[1][1] );
    
    dvdr[0]=real(HB[0][1]);
    dvdr[1]=imag(HB[0][1]);
    dvdr[2]=real(HB[0][0]);
    dvdr[3]=imag(HB[0][0]);

    JI = JInverse(Y[antimatter][i][msw]);

    for(int j=0;j<=2;j++){
      K[antimatter][i][msw][j] = 0.; 
      for(int k=j;k<=3;k++) K[antimatter][i][msw][j] += JI[j][k]*dvdr[k];
      K[antimatter][i][msw][j] *= dr;
    }

    K[antimatter][i][msw][3] = 0.;
    dkkbardr = dkdr(UUbar,dVfMSWbardr);
    dCCdr = CofactorMatricesDerivatives(Hfbar,dVfMSWbardr,dkkbar,dkkbardr);
    QQbar = Q(UUbar,dkkbar,CC,dCCdr);

    K[antimatter][i][msw][4] = (kkbar[0]+QQbar[0])*dr/M_2PI/cgs::constants::hbarc;
    K[antimatter][i][msw][5] = (kkbar[1]+QQbar[1])*dr/M_2PI/cgs::constants::hbarc;

    // *****************************************************************
    // contribution to the self-interaction potential from this energy *
    // *****************************************************************
    Sfm    = UWBW   [i]*Sa   [i][si];
    Sfmbar = UWBWbar[i]*Sabar[i][si];
    VfSIE[i] =     Sfm   *pmatrixm0[    matter][i]*Adjoint(Sfm   )
      - Conjugate( Sfmbar*pmatrixm0[antimatter][i]*Adjoint(Sfmbar) );

  }//end for loop over i

  // ************************************
  // compute self-interaction potential *
  // ************************************
  for(i=0;i<=NE-1;i++){
    VfSI[e ][e ]+=VfSIE[i][e ][e ];
    VfSI[e ][mu]+=VfSIE[i][e ][mu];
    VfSI[mu][e ]+=VfSIE[i][mu][e ];
    VfSI[mu][mu]+=VfSIE[i][mu][mu];
  }

  complex<double> Tr=VfSI[e][e]+VfSI[mu][mu];
  VfSI[e][e]+=Tr;
  VfSI[mu][mu]+=Tr;

  //  VfSI*=NSI*CSI(r);
  VfSIbar=-Conjugate(VfSI);

  // *********************
  // SI part of solution *
  // *********************

#pragma omp parallel for schedule(auto) private(JI) firstprivate(Ha,HB,dvdr)
  for(i=0;i<=NE-1;i++){
    //*********
    // Matter *
    //*********
    Ha = Adjoint(UWBW[i])*VfSI*UWBW[i];

    K[matter][i][si][4]=dr*real(Ha[0][0])/(M_2PI*cgs::constants::hbarc);
    K[matter][i][si][5]=dr*real(Ha[1][1])/(M_2PI*cgs::constants::hbarc);
    
    HB[0][0]=-I/cgs::constants::hbarc*( Ha[0][1]*Sa[i][si][1][0] );
    HB[0][1]=-I/cgs::constants::hbarc*( Ha[0][1]*Sa[i][si][1][1] );
    
    dvdr[0]=real(HB[0][1]);
    dvdr[1]=imag(HB[0][1]);
    dvdr[2]=real(HB[0][0]);
    dvdr[3]=imag(HB[0][0]);
    
    JI=JInverse(Y[matter][i][si]);
    
    for(int j=0;j<=2;j++){
      K[matter][i][si][j]=0.;
      for(int k=j;k<=3;k++) K[matter][i][si][j]+=JI[j][k]*dvdr[k];
      K[matter][i][si][j]*=dr;
    }
    
    K[matter][i][si][3]=0.;
    
    //*************
    // Antimatter *
    //*************
    Ha=Adjoint(UWBWbar[i])*VfSIbar*UWBWbar[i];

    K[antimatter][i][si][4]=dr*real(Ha[0][0])/(M_2PI*cgs::constants::hbarc);
    K[antimatter][i][si][5]=dr*real(Ha[1][1])/(M_2PI*cgs::constants::hbarc);

    HB[0][0]=-I/cgs::constants::hbarc*( Ha[0][1]*Sabar[i][si][1][0] );
    HB[0][1]=-I/cgs::constants::hbarc*( Ha[0][1]*Sabar[i][si][1][1] );

    dvdr[0]=real(HB[0][1]);
    dvdr[1]=imag(HB[0][1]);
    dvdr[2]=real(HB[0][0]);
    dvdr[3]=imag(HB[0][0]);

    JI = JInverse(Y[antimatter][i][si]);

    for(int j=0;j<=2;j++){
      K[antimatter][i][si][j]=0.;
      for(int k=j;k<=3;k++) K[antimatter][i][si][j]+=JI[j][k]*dvdr[k];
      K[antimatter][i][si][j]*=dr;
    }

    K[antimatter][i][si][3]=0.;
  }

}// end of K function


//=======//
// Ebins //
//=======//
vector<double> Ebins(int NE){
  vector<double> energybin(NE);
  for(int i=0;i<NE;i++)
    energybin[i]=1.e6*cgs::units::eV * exp(log(2*1000000) + i*(log(37.48*1000000) - log(2*1000000))/(NE-1))/1000000;
  return energybin;
}

//===========//
// Outputvsr //
//===========//
void Outputvsr(ofstream &fout,
	       ofstream &foutP,
	       ofstream &foutf,
	       ofstream *foutPvsr,
	       ofstream *foutFvsr,
	       double r,
	       vector<vector<vector<vector<double> > > > Y,
	       vector<vector<vector<MATRIX<complex<double>,NF,NF> > > > C0,
	       vector<vector<vector<vector<double> > > > A0,
	       vector<vector<MATRIX<complex<double>,NF,NF> > > Scumulative){

  vector<MATRIX<complex<double>,NF,NF> > VfMSW(NM), dVfMSWdr(NM);
  vector<MATRIX<complex<double>,NF,NF> > VfSI(NM);

  vector<MATRIX<complex<double>,NF,NF> > rhomatrix(NM);

  double rrho=-1; //UNCOMMENT exp(lnrho(log(r)));
  double drrhodr=-1; //UNCOMMENT rrho*lnrho.Derivative(log(r))/r;

  double YYe=-1; //UNCOMMENT Ye(r);
  double dYYedr=-1; //UNCOMMENT Ye.Derivative(r);

  VfMSW[matter][e][e]=Ve(rrho,YYe);
  VfMSW[matter][mu][mu]=Vmu(rrho,YYe);
  VfMSW[antimatter]=-VfMSW[matter];

  dVfMSWdr[matter][e][e]=dVedr(rrho,drrhodr,YYe,dYYedr);
  dVfMSWdr[matter][mu][mu]=dVmudr(rrho,drrhodr,YYe,dYYedr);
  dVfMSWdr[antimatter]=-dVfMSWdr[matter];

  vector<vector<MATRIX<complex<double>,NF,NF> > >
    Hf(NM,vector<MATRIX<complex<double>,NF,NF> >(NE));
  vector<vector<MATRIX<complex<double>,NF,NF> > >
    UU(NM,vector<MATRIX<complex<double>,NF,NF> >(NE));
  vector<vector<vector<MATRIX<complex<double>,NF,NF> > > > 
    WW(NM,vector<vector<MATRIX<complex<double>,NF,NF> > >(NE,vector<MATRIX<complex<double>,NF,NF> >(NS)));
  vector<vector<vector<MATRIX<complex<double>,NF,NF> > > > 
    BB(NM,vector<vector<MATRIX<complex<double>,NF,NF> > >(NE,vector<MATRIX<complex<double>,NF,NF> >(NS)));
  vector<vector<MATRIX<complex<double>,NF,NF> > > 
    Sm(NM,vector<MATRIX<complex<double>,NF,NF> >(NE)), 
    Smf(NM,vector<MATRIX<complex<double>,NF,NF> >(NE)),
    Sf(NM,vector<MATRIX<complex<double>,NF,NF> >(NE));

  vector<vector<vector<double> > > kk(NM,vector<vector<double> >(NE));
  vector<vector<vector<double> > > dkk(NM,vector<vector<double> >(NE));
  vector<double> ePotentialSum(NE),ebarPotentialSum(NE),heavyPotentialSum(NE);
  double totalANuFlux(0.);
  double totalNuFlux(0.);
  double totalHeavyFlux(0.);
  vector<double> Pe(NE),Pebar(NE),Pheavy(NE);
  vector<double> Pvalues(6);
  vector<double> s(6);
  vector<double> predP((NE+2)*(2));

  for(int i=0;i<=NE-1;i++){
    getP(r);
    ePotentialSum[i]=real(pmatrixf0[matter][i][e][e]);//UNCOMMENT eP[i](r);
    ebarPotentialSum[i]=real(pmatrixf0[antimatter][i][e][e]);//UNCOMMENT eBarP[i](r);
    heavyPotentialSum[i]=real(pmatrixf0[matter][i][mu][mu]);
    for(state m=matter;m<=antimatter;m++){
      pmatrixm0[m][i] = Scumulative[m][i]
	* Adjoint(U0[m][i])
	* pmatrixf0[m][i]
	* U0[m][i]
	* Adjoint( Scumulative[m][i] );
    }
  }


  for(int i=0;i<=NE-1;i++){
    //---- matter
    Hf[matter][i]  = HfV[matter][i] + VfMSW[matter];
    kk[matter][i]  = k(Hf[matter][i]);
    dkk[matter][i] = deltak(Hf[matter][i]);
    UU[matter][i]  = U(dkk[matter][i],C0[matter][i],A0[matter][i]);
    
    BB[matter][i][msw] = B(Y[matter][i][msw]);
    WW[matter][i][msw] = W(Y[matter][i][msw]);
    BB[matter][i][si] = B(Y[matter][i][si]);
    WW[matter][i][si] = W(Y[matter][i][si]);
    
    Sm[matter][i] = WW[matter][i][msw]
      * BB[matter][i][msw]
      * WW[matter][i][si]
      * BB[matter][i][si]
      * Scumulative[matter][i];
    Smf[matter][i]= Sm[matter][i] * Adjoint(U0[matter][i]);
    Sf[matter][i] = UU[matter][i] * Smf[matter][i];
    
    //---- antimatter
    Hf[antimatter][i]  = HfV[antimatter][i] + VfMSW[antimatter];
    kk[antimatter][i]  = kbar(Hf[antimatter][i]);
    dkk[antimatter][i] = deltakbar(Hf[antimatter][i]);
    UU[antimatter][i]  = Conjugate(U(dkk[antimatter][i],C0[antimatter][i],A0[antimatter][i]));
    
    BB[antimatter][i][msw] = B(Y[antimatter][i][msw]);
    WW[antimatter][i][msw] = W(Y[antimatter][i][msw]);
    BB[antimatter][i][si] = B(Y[antimatter][i][si]);
    WW[antimatter][i][si] = W(Y[antimatter][i][si]);
    
    Sm[antimatter][i] = WW[antimatter][i][msw]
      * BB[antimatter][i][msw]
      * WW[antimatter][i][si]
      * BB[antimatter][i][si]
      * Scumulative[antimatter][i];
    Smf[antimatter][i]= Sm[antimatter][i] * Adjoint(U0[antimatter][i]);
    Sf[antimatter][i] = UU[antimatter][i] * Smf[antimatter][i];
    
    // compute contribution to self interaction potential
    // scattering matrix matter(electron - x) - scattering matrix Antimatter(antielectron - anti-x)
    // what is VfSI[antimatter]?
    VfSI[matter] += Sf[    matter][i]*pmatrixf0[    matter][i]*Adjoint(Sf[    matter][i])
      - Conjugate(  Sf[antimatter][i]*pmatrixf0[antimatter][i]*Adjoint(Sf[antimatter][i]) );
  }
  
  complex<double> Tr=VfSI[matter][e][e]+VfSI[matter][mu][mu];
  VfSI[matter][e][e]+=Tr;
  VfSI[matter][mu][mu]+=Tr;

  //VfSI[matter]*=NSI*CSI(r);
  VfSI[antimatter] = -Conjugate(VfSI[matter]);

  for(int i=0;i<=NE-1;i++){
    //foutPvsr[i]<<"\n"<<E[i]/(giga*cgs::units::eV);
    foutPvsr[i]<<"\n"<<r;

    foutPvsr[i]<<"\t"<<norm(Sm[matter][i][0][0])<<"\t"<<norm(Sm[matter][i][0][1]);
    foutPvsr[i]<<"\t"<<norm(Sm[matter][i][1][0])<<"\t"<<norm(Sm[matter][i][1][1]);

    foutPvsr[i]<<"\t"<<norm(Sm[antimatter][i][0][0])<<"\t"<<norm(Sm[antimatter][i][0][1]);
    foutPvsr[i]<<"\t"<<norm(Sm[antimatter][i][1][0])<<"\t"<<norm(Sm[antimatter][i][1][1]);

    //foutPvsr[i]<<"\t"<<norm(BB[matter][i][msw][0][0])<<"\t"<<norm(BB[matter][i][msw][0][1]);
    //foutPvsr[i]<<"\t"<<norm(BB[matter][i][msw][1][0])<<"\t"<<norm(BB[matter][i][msw][1][1]);

    //foutPvsr[i]<<"\t"<<norm(BB[antimatter][i][msw][0][0])<<"\t"<<norm(BB[antimatter][i][msw][0][1]);
    //foutPvsr[i]<<"\t"<<norm(BB[antimatter][i][msw][1][0])<<"\t"<<norm(BB[antimatter][i][msw][1][1]);

    //foutPvsr[i]<<"\t"<<norm(BB[matter][i][si][0][0])<<"\t"<<norm(BB[matter][i][si][0][1]);
    //foutPvsr[i]<<"\t"<<norm(BB[matter][i][si][1][0])<<"\t"<<norm(BB[matter][i][si][1][1]);

    //foutPvsr[i]<<"\t"<<norm(BB[antimatter][i][si][0][0])<<"\t"<<norm(BB[antimatter][i][si][0][1]);
    //foutPvsr[i]<<"\t"<<norm(BB[antimatter][i][si][1][0])<<"\t"<<norm(BB[antimatter][i][si][1][1]);

    foutPvsr[i]<<"\t"<<norm(Sf[matter][i][e][e])<<"\t"<<norm(Sf[matter][i][e][mu]);
    foutPvsr[i]<<"\t"<<norm(Sf[matter][i][mu][e])<<"\t"<<norm(Sf[matter][i][mu][mu]);

    foutPvsr[i]<<"\t"<<norm(Sf[antimatter][i][e][e])<<"\t"<<norm(Sf[antimatter][i][e][mu]);
    foutPvsr[i]<<"\t"<<norm(Sf[antimatter][i][mu][e])<<"\t"<<norm(Sf[antimatter][i][mu][mu]);

    foutPvsr[i]<<"\t"<<kk[matter][i][0]<<"\t"<<kk[matter][i][1];
    //foutPvsr[i]<<"\t"<<dkk[matter][i][0];
    foutPvsr[i]<<"\t"<<kk[antimatter][i][0]<<"\t"<<kk[antimatter][i][1];
    //foutPvsr[i]<<"\t"<<dkk[antimatter][i][0];
    if(kV[i][1]>kV[i][0]){
      foutPvsr[i]<<"\t"<<kk[matter][i][1]-kk[matter][i][0];
      foutPvsr[i]<<"\t"<<kk[antimatter][i][1]-kk[antimatter][i][0];
    }
    if(kV[i][0]>kV[i][1]){
      foutPvsr[i]<<"\t"<<kk[matter][i][0]-kk[matter][i][1];
      foutPvsr[i]<<"\t"<<kk[antimatter][i][0]-kk[antimatter][i][1];
    }

    //foutPvsr[i]<<"\t"<<real(UU[matter][i][e][0])<<"\t"<<real(UU[matter][i][e][1]);
    //foutPvsr[i]<<"\t"<<real(UU[matter][i][mu][0])<<"\t"<<real(UU[matter][i][mu][1]);

    //foutPvsr[i]<<"\t"<<real(UU[antimatter][i][e][0])<<"\t"<<real(UU[antimatter][i][e][1]);
    //foutPvsr[i]<<"\t"<<real(UU[antimatter][i][mu][0])<<"\t"<<real(UU[antimatter][i][mu][1]);

    foutPvsr[i]<<"\t"<<YYe;
    foutPvsr[i]<<"\t"<<real(VfMSW[matter][e][e]);
    foutPvsr[i]<<"\t"<<real(VfSI[matter][e ][e ])<<"\t"<<imag(VfSI[matter][e ][e ]);
    foutPvsr[i]<<"\t"<<real(VfSI[matter][e ][mu])<<"\t"<<imag(VfSI[matter][e ][mu]);
    foutPvsr[i]<<"\t"<<real(VfSI[matter][mu][mu])<<"\t"<<imag(VfSI[matter][mu][mu]);

    //vector<double> dkkdr = dkdr(UU[matter][i],dVfMSWdr[matter]);
    //vector<double> QQ = Q(Hf[matter][i],dVfMSWdr[matter],UU[matter][i],kk[matter][i],dkk[matter][i],dkkdr);
    //MATRIX<complex<double>,NF,NF> GG=Gamma(dVfMSWdr[matter],UU[matter][i],kk[matter][i],QQ);
    //foutPvsr[i]<<"\t"<<abs(GG[0][1]);

    //dkkdr = dkdr(UU[antimatter][i],dVfMSWdr[antimatter]);
    //QQ = Q(Hf[antimatter][i],dVfMSWdr[antimatter],UU[antimatter][i],kk[antimatter][i],dkk[antimatter][i],dkkdr);
    //GG=Gamma(dVfMSWdr[antimatter],UU[antimatter][i],kk[antimatter][i],QQ);
    //foutPvsr[i]<<"\t"<<abs(GG[0][1]);

    foutPvsr[i].flush();

    // **************

    rhomatrix[matter]=Smf[matter][i]*pmatrixf0[matter][i]*Adjoint(Smf[matter][i]);
    rhomatrix[antimatter]=Smf[antimatter][i]*pmatrixf0[antimatter][i]*Adjoint(Smf[antimatter][i]);

    //foutFvsr[i]<<"\n"<<E[i]/(giga*cgs::units::eV);
    //foutFvsr[i]<<"\n"<<r;
    //foutFvsr[i]<<"\t"<<real( norm(UV[matter][e][0])*rhomatrix[matter][0][0] +norm(UV[matter][e][1])*rhomatrix[matter][1][1] ) * (0.2*mega*cgs::units::eV)*cgs::constants::c*pow(Rnu/(10.*kilo*cgs::units::pc),2.)
    //        <<"\t"<<real( norm(UV[matter][mu][0])*rhomatrix[matter][0][0] +norm(UV[matter][mu][1])*rhomatrix[matter][1][1] ) * (0.2*mega*cgs::units::eV)*cgs::constants::c*pow(Rnu/(10.*kilo*cgs::units::pc),2.);
    //foutFvsr[i]<<"\t"<<real( norm(UV[antimatter][e][0])*rhomatrix[antimatter][0][0] +norm(UV[antimatter][e][1])*rhomatrix[antimatter][1][1] ) * (0.2*mega*cgs::units::eV)*cgs::constants::c*pow(Rnu/(10.*kilo*cgs::units::pc),2.)
    //        <<"\t"<<real( norm(UV[antimatter][mu][0])*rhomatrix[antimatter][0][0] +norm(UV[antimatter][mu][1])*rhomatrix[antimatter][1][1] ) * (0.2*mega*cgs::units::eV)*cgs::constants::c*pow(Rnu/(10.*kilo*cgs::units::pc),2.);

    foutFvsr[i].flush();
  }

  /////////////////////////////////////////////////////////////////////

  for(int i=0;i<NE;i++){
    fout<<norm(Sf[    matter][i][e ][e ])<<"\t"<<norm(Sf[antimatter][i][e ][e ])<<"\t";
    fout<<real(Sf[    matter][i][e ][e ])<<"\t"<<imag(Sf[    matter][i][e ][e ])<<"\t";
    fout<<real(Sf[antimatter][i][e ][e ])<<"\t"<<imag(Sf[antimatter][i][e ][e ])<<"\t";
    fout<<real(Sf[    matter][i][e ][mu])<<"\t"<<imag(Sf[    matter][i][e ][mu])<<"\t";
    fout<<real(Sf[antimatter][i][e ][mu])<<"\t"<<imag(Sf[antimatter][i][e ][mu])<<"\t";
    fout<<real(Sf[    matter][i][mu][e ])<<"\t"<<imag(Sf[    matter][i][mu][e ])<<"\t";
    fout<<real(Sf[antimatter][i][mu][e ])<<"\t"<<imag(Sf[antimatter][i][mu][e ])<<"\t";
    fout<<real(Sf[    matter][i][mu][mu])<<"\t"<<imag(Sf[    matter][i][mu][mu])<<"\t";
    fout<<real(Sf[antimatter][i][mu][mu])<<"\t"<<imag(Sf[antimatter][i][mu][mu])<<"\t";
    fout<<norm(Sf[    matter][i][mu][mu])<<"\t"<<norm(Sf[antimatter][i][mu][mu])<<"\t";
    fout<<dm21/(4*E[i])*cgs::constants::c4<<"\t";

    fout.flush();

    Pe    [i] = norm(Sf[    matter][i][e ][e ]);
    Pebar [i] = norm(Sf[antimatter][i][e ][e ]);
    Pheavy[i] = norm(Sf[    matter][i][mu][mu]);
  }

  foutP<<r<<"\t"<<Ve(rrho,YYe)<<"\t";//1,2
  foutP<<real(VfSI[    matter][e ][e ])<<"\t"<<imag(VfSI[    matter][e ][e ])<<"\t";
  foutP<<real(VfSI[    matter][mu][mu])<<"\t"<<imag(VfSI[    matter][mu][mu])<<"\t";//3,4,5,6
  fout <<real(VfSI[    matter][e ][e ])<<"\t"<<imag(VfSI[    matter][e ][e ])<<"\t";
  fout <<real(VfSI[    matter][mu][mu])<<"\t"<<imag(VfSI[    matter][mu][mu]);
  fout <<real(VfSI[antimatter][e ][e ])<<"\t"<<imag(VfSI[antimatter][e ][e ])<<"\t";
  fout <<real(VfSI[antimatter][mu][mu])<<"\t"<<imag(VfSI[antimatter][mu][mu]);

  Pvalues = averageProbability(Pe,Pebar,Pheavy,ebarPotentialSum,ePotentialSum,heavyPotentialSum);
  totalNuFlux = Pvalues[3];
  totalANuFlux =Pvalues[4];
  totalHeavyFlux = Pvalues[5];
  foutP<<totalNuFlux<<"\t";//Nu,7
  foutP<<totalANuFlux<<"\t";//ANu,8
  foutP<<Pvalues[5]<<"\t";//Heavy,9
  foutP<<Pvalues[0]<<"\t"<<Pvalues[1]<<"\t"<<Pvalues[2]<<"\t";//Pe,Pebar,Pheavy;10,11,12

  predP=predictProbability(Pvalues[3],Pvalues[4],Ve(rrho,YYe),E,ebarPotentialSum,ePotentialSum,heavyPotentialSum);
  foutP<<predP[0]<<"\t"<<predP[1+NE]<<"\t";//13,14
  for(int i=0;i<NE;i++) foutP<<predP[1+i]<<"\t"<<predP[(NE+1)+i+1]<<"\t";//15,16,...2*(NE-1)+15,
  foutP<<predP[(NE+1)*2]<<"\t"<<predP[(NE+1)*2+1]<<"\t";//2*(NE-1)+17,2*(NE-1)+18
  fout<<endl;
  foutP<<endl;
  foutP.flush();
  fout.flush();

  foutf << r << "\t";
  for(int i=0; i<NE; i++)
    for(state m=matter; m<=antimatter; m++)
      for(flavour f1=e; f1<=mu; f1++)
	for(flavour f2=e; f2<=mu; f2++) {
	  foutf << real( fmatrixf[m][i][f1][f2] ) << "\t";
	  foutf << imag( fmatrixf[m][i][f1][f2] ) << "\t";
	}
  foutf << endl;
  foutf.flush();
}

#include "headers/update.h"
