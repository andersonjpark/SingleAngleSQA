//==========//
// UpdateSM //
//==========//
vector<vector<MATRIX<complex<double>,NF,NF> > > 
UpdateSm(double rminus,
	 double rplus,
	 double Yeminus,
	 double Yeplus,
	 vector<vector<vector<vector<double> > > > Y,
	 vector<vector<vector<MATRIX<complex<double>,NF,NF> > > > C0,
	 vector<vector<vector<vector<double> > > > A0,vector<vector<MATRIX<complex<double>,NF,NF> > > Smprior){

  vector<MATRIX<complex<double>,NF,NF> > VfMSW(NM,MATRIX<complex<double>,NF,NF>(NF,NF));
  MATRIX<complex<double>,NF,NF> Hf,Hfbar;
  MATRIX<complex<double>,NF,NF> UU,UUbar; 

  vector<double> kk,kkbar,dkk,dkkbar;

  vector<MATRIX<complex<double>,NF,NF> > CC(NF);
  vector<vector<double> > AA(NF,vector<double>(NF));

  vector<vector<MATRIX<complex<double>,NF,NF> > > Sm(NM,vector<MATRIX<complex<double>,NF,NF> >(NE));

  double rrho, YYe;

  // multiply SMSW by the mixing matrix at rminus
  rrho = exp(lnrho(log(rminus)));
  YYe = Ye(rminus); 
    
  VfMSW[matter][e][e] = Ve(rrho,YYe); 
  VfMSW[matter][mu][mu] = Vmu(rrho,YYe); 

  VfMSW[antimatter] = -Conjugate(VfMSW[matter]);

  for(int i=0;i<=NE-1;i++){
    Hf=HfV[matter][i]+VfMSW[matter];
    kk=k(Hf);
    dkk=deltak(Hf);
    UU=U(dkk,C0[matter][i],A0[matter][i]);    
    Sm[matter][i] = UU
      * W(Y[matter][i][msw])
      * B(Y[matter][i][msw])
      * W(Y[matter][i][si])
      * B(Y[matter][i][si])
      * Smprior[matter][i];
    
    Hfbar = HfV[antimatter][i]+VfMSW[antimatter];
    kkbar = kbar(Hfbar);
    dkkbar = deltakbar(Hfbar);
    UUbar = Conjugate(U(dkkbar,C0[antimatter][i],A0[antimatter][i]));
    Sm[antimatter][i] = UUbar
      * W(Y[antimatter][i][msw])
      * B(Y[antimatter][i][msw])
      * W(Y[antimatter][i][si])
      * B(Y[antimatter][i][si])
      * Smprior[antimatter][i];
  }
  
  // multiply SMSW by the adjoint of the mixing matrix at rplus
  rrho = exp(lnrho(log(rplus)));
  YYe = Ye(rplus); 
  
  VfMSW[matter][e][e] = Ve(rrho,YYe); 
  VfMSW[matter][mu][mu] = Vmu(rrho,YYe); 

  VfMSW[antimatter] = -Conjugate(VfMSW[matter]);

  for(int i=0;i<=NE-1;i++){
    Hf = HfV[matter][i]+VfMSW[matter];
    kk = k(Hf);
    dkk = deltak(Hf);
    CC = CofactorMatrices(Hf,kk);
    AA = MixingMatrixFactors(CC,C0[matter][i],A0[matter][i]);
    UU = U(dkk,CC,AA);    
    Sm[matter][i] = Adjoint(UU)*MATRIX<complex<double>,NF,NF>(Sm[matter][i]); 
    
    Hfbar = HfV[antimatter][i] + VfMSW[antimatter];
    kkbar = kbar(Hfbar);
    dkkbar = deltakbar(Hfbar);
    CC = CofactorMatrices(Hfbar,kkbar);
    AA = MixingMatrixFactors(CC,C0[antimatter][i],A0[antimatter][i]);
    UUbar = Conjugate(U(dkkbar,CC,AA));
    Sm[antimatter][i] = Adjoint(UUbar)*MATRIX<complex<double>,NF,NF>(Sm[antimatter][i]);
  }
  
  return Sm;
}

//=========//
// UpdateC //
//=========//
vector<vector<vector<MATRIX<complex<double>,NF,NF> > > > UpdateC(double r,double Ye){
  vector<MATRIX<complex<double>,NF,NF> > VfMSW(NM);
  MATRIX<complex<double>,NF,NF> Hf,Hfbar;
  vector<double> kk,kkbar;

  vector<vector<vector<MATRIX<complex<double>,NF,NF> > > > 
    CC(NM,vector<vector<MATRIX<complex<double>,NF,NF> > >(NE,vector<MATRIX<complex<double>,NF,NF> >(NF)));

  double rrho = exp(lnrho(log(r)));
    
  VfMSW[matter][e][e] = Ve(rrho,Ye); 
  VfMSW[matter][mu][mu] = Vmu(rrho,Ye); 

  VfMSW[antimatter] = -VfMSW[matter];

  int i;
#pragma omp parallel for schedule(auto) private(Hf,Hfbar,kk,kkbar)
  for(i=0;i<=NE-1;i++){
    Hf = HfV[matter][i] + VfMSW[matter]; 
    kk = k(Hf);
    CC[matter][i] = CofactorMatrices(Hf,kk);
    
    Hfbar = HfV[antimatter][i]+VfMSW[antimatter];
    kkbar = kbar(Hfbar);
    CC[antimatter][i] = CofactorMatrices(Hfbar,kkbar);
  }

  return CC;
}

//=========//
// UpdateA //
//=========//
vector<vector<vector<vector<double> > > >
UpdateA(vector<vector<vector<MATRIX<complex<double>,NF,NF> > > > C,
	vector<vector<vector<MATRIX<complex<double>,NF,NF> > > > C0,
	vector<vector<vector<vector<double> > > > A0){

  vector<vector<vector<vector<double> > > > 
    A(NM,vector<vector<vector<double> > >(NE,vector<vector<double> >(NF,vector<double>(NF))));
  
  int i;
#pragma omp parallel for schedule(auto)
  for(i=0;i<=NE-1;i++){
    A[matter][i]=MixingMatrixFactors(C[matter][i],C0[matter][i],A0[matter][i]);
    A[antimatter][i]=MixingMatrixFactors(C[antimatter][i],C0[antimatter][i],A0[antimatter][i]);
  }
  
  return A;
}

