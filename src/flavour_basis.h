#ifndef FLAVORBASIS_H
#define FLAVORBASIS_H

#include "misc.h"

array<array<MATRIX<complex<double>,NF,NF>,NE>,NM>
  Evaluate_VfVac(const array<array<double,NF>,NE>& kV,
	       const array<MATRIX<complex<double>,NF,NF>,NM>& UV){
  
  MATRIX<complex<double>,NF,NF> KV;
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> VfVac;
  
  for(int i=0;i<=NE-1;i++){
    for(int j=0;j<=NF-1;j++) KV[j][j]=kV[i][j];
    VfVac[matter][i]=UV[matter]*KV*Adjoint(UV[matter]);
    VfVac[antimatter][i]=Conjugate(VfVac[matter][i]);
  }
  return VfVac;
}

#endif
