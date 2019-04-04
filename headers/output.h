#ifndef OUTPUT_H
#define OUTPUT_H

void Outputvsr(ofstream &fout,
	       ofstream &foutP,
	       ofstream &foutf,
	       ofstream &foutdangledr,
	       const State& s,
	       const array<DISCONTINUOUS,NE>& eP,
	       const array<DISCONTINUOUS,NE>& eBarP,
	       const array<DISCONTINUOUS,NE>& xP){

  array<double,NE> ePotentialSum,ebarPotentialSum,heavyPotentialSum;
  array<double,NE> Pe,Pebar,Pheavy;
  for(int i=0;i<=NE-1;i++){
    ePotentialSum[i]=eP[i](s.r);
    ebarPotentialSum[i]=eBarP[i](s.r);
    heavyPotentialSum[i]=xP[i](s.r);
    Pe    [i] = norm(s.Sf[matter][i][e ][e ]);
    Pebar [i] = norm(s.Sf[antimatter][i][e ][e ]);
    Pheavy[i] = norm(s.Sf[matter][i][mu][mu]);
  }


  fout << s.r << "\t";
  for(int i=0; i<NE; i++)
    for(state m=matter; m<=antimatter; m++){
      for(flavour f1=e; f1<=mu; f1++)
	for(flavour f2=e; f2<=mu; f2++) {
	  fout << real(s.Sf[m][i][f1][f2] ) << "\t";
	  fout << imag(s.Sf[m][i][f1][f2] ) << "\t";
	}
    }
  fout << endl;
  fout.flush();

  foutP<<s.r<<"\t"<<Ve(s.rho,s.Ye)<<"\t";//1,2
  foutP<<real(s.VfSI[    matter][e ][e ])<<"\t"<<real(s.VfSI[    matter][mu][mu])<<"\t";
  foutP<<real(s.VfSI[antimatter][e ][e ])<<"\t"<<real(s.VfSI[antimatter][mu][mu])<<"\t";//3,4,5,6
  array<double,6> Pvalues = averageProbability(Pe,Pebar,Pheavy,ebarPotentialSum,ePotentialSum,heavyPotentialSum);
  double totalNuFlux = Pvalues[3];
  double totalANuFlux =Pvalues[4];
  double totalHeavyFlux = Pvalues[5];
  foutP<<totalNuFlux<<"\t";//Nu,7
  foutP<<totalANuFlux<<"\t";//ANu,8
  foutP<<Pvalues[5]<<"\t";//Heavy,9
  foutP<<Pvalues[0]<<"\t"<<Pvalues[1]<<"\t"<<Pvalues[2]<<"\t";//Pe,Pebar,Pheavy;10,11,12

  array<double,(NE+2)*2> predP=predictProbability(Pvalues[3],Pvalues[4],Ve(s.rho,s.Ye),E,ebarPotentialSum,ePotentialSum,heavyPotentialSum);
  foutP<<predP[0]<<"\t"<<predP[1+NE]<<"\t";//13,14
  for(int i=0;i<NE;i++) foutP<<predP[1+i]<<"\t"<<predP[(NE+1)+i+1]<<"\t";//15,16,...2*(NE-1)+15,
  foutP<<predP[(NE+1)*2]<<"\t"<<predP[(NE+1)*2+1]<<"\t";//2*(NE-1)+17,2*(NE-1)+18
  foutP<<endl;
  foutP.flush();

  foutf << s.r << "\t";
  for(int i=0; i<NE; i++)
    for(state m=matter; m<=antimatter; m++)
      for(flavour f1=e; f1<=mu; f1++)
	for(flavour f2=e; f2<=mu; f2++) {
	  foutf << real( s.fmatrixf[m][i][f1][f2] ) << "\t";
	  foutf << imag( s.fmatrixf[m][i][f1][f2] ) << "\t";
	}
  

  foutf << endl;
  foutf.flush();

  foutdangledr << s.r << "\t";
  for(state m=matter; m<=antimatter; m++)
    for(int i=0; i<NE; i++)
      foutdangledr << s.dtheta_dr_osc[i][m] << "\t";
  for(state m=matter; m<=antimatter; m++)
    for(int i=0; i<NE; i++)
      foutdangledr << s.dphi_dr_osc[i][m] << "\t";
  for(state m=matter; m<=antimatter; m++)
    for(int i=0; i<NE; i++)
      foutdangledr << s.dtheta_dr_interact[i][m] << "\t";
  for(state m=matter; m<=antimatter; m++)
    for(int i=0; i<NE; i++)
      foutdangledr << s.dphi_dr_interact[i][m] << "\t";
  foutdangledr << endl;
  foutdangledr.flush();
}

#endif
