
//=========//
// UpdateC //
//=========//
array<array<array<MATRIX<complex<double>,NF,NF>,NF>,NE>,NM> UpdateC(const State& s, const DISCONTINUOUS& lnrho, const DISCONTINUOUS& Ye){
  vector<MATRIX<complex<double>,NF,NF> > VfMSW(NM);
  MATRIX<complex<double>,NF,NF> Hf,Hfbar;
  array<double,NF> kk,kkbar;

  array<array<array<MATRIX<complex<double>,NF,NF>,NF>,NE>,NM> CC;
  VfMSW[matter][e][e] = Ve(s.rho,s.Ye); 
  VfMSW[matter][mu][mu] = Vmu(s.rho,s.Ye); 

  VfMSW[antimatter] = -VfMSW[matter];

  int i;
#pragma omp parallel for schedule(auto) private(Hf,Hfbar,kk,kkbar)
  for(i=0;i<=NE-1;i++){
    Hf = s.HfV[matter][i] + VfMSW[matter]; 
    kk = k(Hf);
    CC[matter][i] = CofactorMatrices(Hf,kk);
    
    Hfbar = s.HfV[antimatter][i]+VfMSW[antimatter];
    kkbar = kbar(Hfbar);
    CC[antimatter][i] = CofactorMatrices(Hfbar,kkbar);
  }

  return CC;
}

//=========//
// UpdateA //
//=========//
array<array<array<array<double,NF>,NF>,NE>,NM>
UpdateA(array<array<array<MATRIX<complex<double>,NF,NF>,NF>,NE>,NM>& C,
	array<array<array<MATRIX<complex<double>,NF,NF>,NF>,NE>,NM>& C0,
	array<array<array<array<double,NF>,NF>,NE>,NM>& A0){

  array<array<array<array<double,NF>,NF>,NE>,NM> A;
  
  int i;
#pragma omp parallel for schedule(auto)
  for(i=0;i<=NE-1;i++){
    A[matter][i]=MixingMatrixFactors(C[matter][i],C0[matter][i],A0[matter][i]);
    A[antimatter][i]=MixingMatrixFactors(C[antimatter][i],C0[antimatter][i],A0[antimatter][i]);
  }
  
  return A;
}

