#include <vector>
/*************************************************************************
*       Yonglin Zhu yzhu14@ncsu.edu
*       March,2016
*
*	Energy bins for Brett Deaton's simple disk model
*	Average Probability weighted by self-interaction potential at each 
*	energy level.(in single angle approximation only)
**************************************************************************/
//Energy bins
/*
vector<double>
Ebins(int NE){
  vector<double> energybin(NE);
  for(int i=0;i<NE;i++)
    energybin[i]=50/32+i*15*50/16/(NE-1);//brett simple disk
  return energybin;
}
*/

//====================//
// averageProbability //
//====================//
array<double,6>
averageProbability(array<double,NE> Pe,
		   array<double,NE> Pebar,
		   array<double,NE> Pheavy,
		   array<double,NE> ebarPotentialSum,
		   array<double,NE> ePotentialSum,
		   array<double,NE> heavyPotentialSum){

  array<double,6> Pvalues;
  double totalANu(0.);
  double totalNu(0.);
  double totalHeavy(0.);
  int NE(Pe.size());

  for(int i=0;i<NE;i++){
    Pvalues[0] += Pe[i] * ePotentialSum[i];
    Pvalues[1] += Pebar[i] * ebarPotentialSum[i];
    Pvalues[2] += Pheavy[i] * heavyPotentialSum[i];
    totalANu   += ebarPotentialSum[i];
    totalNu    += ePotentialSum[i];
    totalHeavy += heavyPotentialSum[i];
  }
  
  Pvalues[0] /= totalNu;
  Pvalues[1] /= totalANu;
  Pvalues[2] /= totalHeavy;
  Pvalues[3]  = totalNu;
  Pvalues[4]  = totalANu;
  Pvalues[5]  = totalHeavy;

  return Pvalues;
}

//patch for to_string
namespace patch{
  template < typename T > std::string to_string(const T& n){
    std::ostringstream stm;
    stm<<n;
    return stm.str();
  }
}
