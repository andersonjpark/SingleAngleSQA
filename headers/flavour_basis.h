#ifndef FLAVORBASIS_H
#define FLAVORBASIS_H

#include "misc.h"

void Evaluate_HfV(const array<array<double,NF>,NE>& kV){
  MATRIX<complex<double>,NF,NF> KV;
  
  for(int i=0;i<=NE-1;i++){
    for(int j=0;j<=NF-1;j++) KV[j][j]=kV[i][j];
    HfV[matter][i]=UV[matter]*KV*Adjoint(UV[matter]);
    HfV[antimatter][i]=Conjugate(HfV[matter][i]);
  }
}

#endif
