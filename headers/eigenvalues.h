#ifndef EIGENVALUES_H
#define EIGENVALUES_H

#include "misc.h"

//===//
// D //
//===//
double D(MATRIX<complex<double>,NF,NF> Hf){
  return 4.*norm(Hf[e][mu]) +norm(Hf[e][e]-Hf[mu][mu]);
}

//=====//
// f,g //
//=====//
double f(double x){
  return 0.5*x*(1.-0.25*x*(1.-0.5*x*(1.-0.625*x)));
}
double g(double x){
  return 0.5*x*(1.-0.75*x*(1.-0.83333333*x*(1.-0.875*x)));
}

//====//
// k1 //
//====//
double k1(double T,double sqrtD){
  return 0.5*(T+a1*sqrtD);
}
double k1bar(double Tbar,double sqrtDbar){
  return 0.5*(Tbar+a1*sqrtDbar);
}
double asymptotick1(double Hee,double Hmm,double x,int s){
  return (1+a1*s)/2.*Hee +(1-a1*s)/2.*Hmm +a1*s/2.*(Hee-Hmm)*f(x);
}
double asymptotick1bar(double Hee,double Hmm,double x,int s){
  return (1+a1*s)/2.*Hee +(1-a1*s)/2.*Hmm +a1*s/2.*(Hee-Hmm)*f(x);
}

//====//
// k2 //
//====//
double k2(double T,double sqrtD){
  return 0.5*(T+a2*sqrtD);
}
double k2bar(double Tbar,double sqrtDbar){
  return 0.5*(Tbar+a2*sqrtDbar);
}
double asymptotick2(double Hee,double Hmm,double x,int s){
  return (1+a2*s)/2.*Hee +(1-a2*s)/2.*Hmm +a2*s/2.*(Hee-Hmm)*f(x);
}
double asymptotick2bar(double Hee,double Hmm,double x,int s){
  return (1+a2*s)/2.*Hee +(1-a2*s)/2.*Hmm +a2*s/2.*(Hee-Hmm)*f(x);
}

// *********************************************************************

//===//
// k //
//===//
array<double,NF> k(MATRIX<complex<double>,NF,NF> Hf){
  array<double,NF> k;
  double t=real(Trace(Hf)), sqrtd=sqrt(D(Hf));
  
  k[0]=k1(t,sqrtd);
  k[1]=k2(t,sqrtd);
  return k;
}

array<double,NF> k(double t,double sqrtd){
  array<double,NF> k;
  k[0]=k1(t,sqrtd);
  k[1]=k2(t,sqrtd);
  return k;
}

//======//
// kbar //
//======//
array<double,NF> kbar(MATRIX<complex<double>,NF,NF> Hf){
  array<double,NF> k;
  double tbar=real(Trace(Hf)), sqrtdbar=sqrt(D(Hf));
  
  k[0]=k1bar(tbar,sqrtdbar);
  k[1]=k2bar(tbar,sqrtdbar);
  
  return k;
}

array<double,NF> kbar(double tbar,double sqrtdbar){
  array<double,NF> k;
  k[0]=k1bar(tbar,sqrtdbar);
  k[1]=k2bar(tbar,sqrtdbar);
  return k;
}

//=============//
// asymptotick //
//=============//
array<double,NF> asymptotick(MATRIX<complex<double>,NF,NF> Hf){
  array<double,NF> k;
  double x=4.*norm(Hf[e][mu])/norm(Hf[e][e]-Hf[mu][mu]);
  int s=static_cast<int>(Sign(real(Hf[e][e]-Hf[mu][mu])));
  k[0]=asymptotick1(real(Hf[e][e]),real(Hf[mu][mu]),x,s);
  k[1]=asymptotick2(real(Hf[e][e]),real(Hf[mu][mu]),x,s);
  return k;
}

//================//
// asymptotickbar //
//================//
array<double,NF> asymptotickbar(MATRIX<complex<double>,NF,NF> Hf){
  array<double,NF> k;
  double x=4.*norm(Hf[e][mu])/norm(Hf[e][e]-Hf[mu][mu]); int s=static_cast<int>(Sign(real(Hf[e][e]-Hf[mu][mu])));
  k[0]=asymptotick1bar(real(Hf[e][e]),real(Hf[mu][mu]),x,s);
  k[1]=asymptotick2bar(real(Hf[e][e]),real(Hf[mu][mu]),x,s);
  return k;
}

//========//
// deltak //
//========//
array<double,NF-1> deltak(MATRIX<complex<double>,NF,NF> Hf){
  array<double,NF-1> dk;
  double d=D(Hf), sqrtd=sqrt(d);
  dk[0]=a1*sqrtd;
  return dk;
}

array<double,NF-1> deltak(double sqrtd){
  array<double,NF-1> dk;
  dk[0]=a1*sqrtd;
  return dk;
}

//===========//
// deltakbar //
//===========//
array<double,NF-1> deltakbar(MATRIX<complex<double>,NF,NF> Hfbar){
  array<double,NF-1> dk;
  double dbar=D(Hfbar), sqrtdbar=sqrt(dbar);
  dk[0]=a1*sqrtdbar;
  return dk;
}

array<double,NF-1> deltakbar(double sqrtdbar){
  array<double,NF-1> dk;
  dk[0]=a1*sqrtdbar;
  return dk;
}

//======//
// dkdr //
//======//
array<double,NF> dkdr(MATRIX<complex<double>,NF,NF> U,MATRIX<complex<double>,NF,NF> dHfdr){
  array<double,NF> dkdx;
  for(flavour f=e;f<=mu;f++)
    for(flavour g=e;g<=mu;g++){
      dkdx[0]+=real(conj(U[f][0])*dHfdr[f][g]*U[g][0]); 
      dkdx[1]+=real(conj(U[f][1])*dHfdr[f][g]*U[g][1]);
    }
  return dkdx;
}  

#endif
