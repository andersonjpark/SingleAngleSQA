#ifndef FLAVORBASIS_H
#define FLAVORBASIS_H

#include "misc.h"

array<array<MATRIX<complex<double>,NF,NF>,NE>,NM>
  Evaluate_HfV(const array<array<double,NF>,NE>& kV,
	       const array<MATRIX<complex<double>,NF,NF>,NM>& UV){
  
  MATRIX<complex<double>,NF,NF> KV;
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> HfV;
  
  for(int i=0;i<=NE-1;i++){
    for(int j=0;j<=NF-1;j++) KV[j][j]=kV[i][j];
    HfV[matter][i]=UV[matter]*KV*Adjoint(UV[matter]);
    HfV[antimatter][i]=Conjugate(HfV[matter][i]);
  }
  return HfV;
}

#endif
