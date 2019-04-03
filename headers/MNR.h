#include<vector>
#include<iostream>
/*************************************************************************
*	Yonglin Zhu yzhu14@ncsu.edu
*	March,2016
*
*	Prediction of survival probability of electron with corrections including 
*	hierarchies.
*	neutrino and anti-neutrino during Matter-Neutrino Resonance(MNR)
*
*
*%\cite{Vaananen:2015hfa}
\bibitem{Vaananen:2015hfa} 
  D.~Vaananen and G.~C.~McLaughlin,
  %``Uncovering the Matter-Neutrino Resonance,''
  Phys.\ Rev.\ D {\bf 93}, no. 10, 105044 (2016)
  doi:10.1103/PhysRevD.93.105044
  [arXiv:1510.00751 [hep-ph]].
  %%CITATION = doi:10.1103/PhysRevD.93.105044;%%
  %3 citations counted in INSPIRE as of 29 Jul 2016
**************************************************************************/
//prediction without considering the hierarchy
array<double,(NE+2)*(2)>
predictProbability(double totalNu,
		   double totalANu,
		   double Ve){

  array<double,(NE+2)*(2)> predP;
  double alpha(totalANu/totalNu);
  double mu(totalNu);
  predP[0]=(1+(alpha*alpha-(Ve/mu)*(Ve/mu)-1)/(2*(Ve/mu)))/2;
  predP[1]=(1+(alpha*alpha+(Ve/mu)*(Ve/mu)-1)/(2*alpha*(Ve/mu)))/2;
  
  return predP;
}

//prediction considering the hierarchy
array<double,(NE+2)*(2)>
predictProbability(double totalNu,
		   double totalANu,
		   double Ve,
		   const vector<double>& E,
		   const array<double,NE>& ebarPotentialSum,
		   const array<double,NE>& ePotentialSum,
		   const array<double,NE>& heavyPotentialSum){
  
  array<double,(NE+2)*(2)> predP;
  double alpha(totalANu/totalNu);
  double mu(totalNu);
  array<double,NE> Delta;
  array<double,NE> epsilon;
  double sume(0.0);

  for (int i=0;i<=NE-1;i++){
    //predicted Probability at each energy level,to get <P()>.
    Delta[i] = dm21 / (2*E[i]) * cgs::constants::c4;      	
    epsilon[i] = 2. * Delta[i] * cos(2.*theta12V) / (Ve + Delta[i]*cos(2.*theta12V));
    sume += epsilon[i];
    predP[       i+1] = 0.5 * (1. + (alpha*alpha - (Ve/mu)*(Ve/mu) - (1.-epsilon[i])*(1.-epsilon[i])) / (2.*(Ve/mu)*(1.-epsilon[i])));
    predP[(NE+1)+i+1] = 0.5 * (1. + (alpha*alpha + (Ve/mu)*(Ve/mu) - (1.-epsilon[i])*(1.-epsilon[i])) / (2.*alpha*(Ve/mu)));
    //to get <Pee()>*totalNu, <Pebarebar()>*totalANu
    predP[0] += predP[i+1] * ePotentialSum[i];
    predP[NE+1] += predP[(NE+1)+i+1] * ebarPotentialSum[i];
  }

  double meane = sume / epsilon.size();
  //to get <Pee()>, <Pebarebar()>
  predP[0] /= totalNu;
  predP[1+NE] /= totalANu;
  //to get Pee(<>),Pebarebar(<>)
  predP[(NE+1)*2]   = (1.+(alpha*alpha-(Ve/mu)*(Ve/mu)-(1.-meane)*(1.-meane))/(2.*(Ve/mu)*(1.-meane)))/2.;
  predP[(NE+1)*2+1] = (1.+(alpha*alpha+(Ve/mu)*(Ve/mu)-(1.-meane)*(1.-meane))/(2.*alpha*(Ve/mu)))/2.;
  
  return predP;
}

