#if !defined(_FLUX)
#define _FLUX

//=========================//
// Initialize_Luminosities //
//=========================//
void Initialize_Luminosities(string spectrapath, double t){
  DISCONTINUOUS Lnue,Lanue,Lnux,Lanux;
  
  Lnue.Open((spectrapath  + string("L.nue.dat")).c_str());
  Lanue.Open((spectrapath + string("L.anue.dat")).c_str());
  Lnux.Open((spectrapath  + string("L.nux.dat")).c_str());
  Lanux.Open((spectrapath + string("L.anux.dat")).c_str());
  
  L[matter][e]       = Lnue(t);  // ergs/s      
  L[antimatter][e]   = Lanue(t); 
  L[matter][mu]      = Lnux(t);        
  L[antimatter][mu]  = Lanux(t); 
  //L[matter][tau]     = L[matter][mu]; 
  //L[antimatter][tau] = L[antimatter][mu];
  
  //L[matter][e]=0.;         L[antimatter][e]=0.; // ergs/s
  //L[matter][mu]=0.;        L[antimatter][mu]=0.; 
  //L[matter][tau]=0.;       L[antimatter][tau]=0.;
  
  cout << endl ;
  cout << "Le\t"      << L[matter    ][e ];
  cout << "\tLae\t"   << L[antimatter][e ];
  cout << "\tLnux\t"  << L[matter    ][mu];
  cout << "\tLanux\t" << L[antimatter][mu];
  cout << " erg/s." << endl;     
}

//=========================//
// Initialize_MeanEnergies //
//=========================//
void Initialize_MeanEnergies(string spectrapath, double t){
  DISCONTINUOUS meanEnue,meanEanue,meanEnux,meanEanux;

  meanEnue.Open((spectrapath  + string("meanE.nue.dat")).c_str());
  meanEanue.Open((spectrapath + string("meanE.anue.dat")).c_str());
  meanEnux.Open((spectrapath  + string("meanE.nux.dat")).c_str());
  meanEanux.Open((spectrapath + string("meanE.anux.dat")).c_str());
  
  meanE[matter][e]       = meanEnue(t);         
  meanE[antimatter][e]   = meanEanue(t); // MeV
  meanE[matter][mu]      = meanEnux(t);       
  meanE[antimatter][mu]  = meanEanux(t);  
  //meanE[matter][tau]     = meanE[matter][mu]; 
  //meanE[antimatter][tau] = meanE[antimatter][mu];
  
  cout << endl;
  cout << "Ee\t"      << meanE[matter    ][e ];
  cout << "\tEae\t"   << meanE[antimatter][e ];
  cout << "\tEnux\t"  << meanE[matter    ][mu];
  cout << "\tEanux\t" << meanE[antimatter][mu];
  cout << " MeV." << endl; 
  
  for(state m=matter;m<=antimatter;m++)
    for(flavour f=e;f<=mu;f++)
      meanE[m][f] *= 1e6*cgs::units::eV;
}

//============================//
// Initialize_PinchParameters //
//============================//
void Initialize_PinchParameters(string spectrapath, double t){
  DISCONTINUOUS rmsEnue,rmsEanue,rmsEnux,rmsEanux;
  
  rmsEnue.Open((spectrapath  + string("rmsE.nue.dat")).c_str());
  rmsEanue.Open((spectrapath + string("rmsE.anue.dat")).c_str());
  rmsEnux.Open((spectrapath  + string("rmsE.nux.dat")).c_str());
  rmsEanux.Open((spectrapath + string("rmsE.anux.dat")).c_str());
  
  double CE;
  CE = rmsEnue(t)*1e6*cgs::units::eV/meanE[matter][e] ;
  pinch[matter][e] = (2.-CE*CE)/(CE*CE-1.);
  
  CE = rmsEanue(t)*1e6*cgs::units::eV/meanE[antimatter][e] ;
  pinch[antimatter][e] = (2.-CE*CE)/(CE*CE-1.);
  
  CE = rmsEnux(t)*1e6*cgs::units::eV/meanE[matter][mu];
  pinch[matter][mu]  = (2.-CE*CE)/(CE*CE-1.); 
  //pinch[matter][tau] = pinch[matter][mu];
  
  CE = rmsEanux(t)*1e6*cgs::units::eV/meanE[antimatter][mu];
  pinch[antimatter][mu]  = (2.-CE*CE)/(CE*CE-1.); 
  //pinch[antimatter][tau] = pinch[antimatter][mu];
  
  cout << endl;
  cout << "rmsE\t"       << rmsEnue(t);
  cout << "\trmsEae\t"   << rmsEanue(t);
  cout << "\trmsEnux\t"  << rmsEnux(t);
  cout << "\trmsEanux\t" << rmsEanux(t);
  cout << " eV." << endl;
  
  cout << "pinch e\t"    << pinch[matter    ][e];
  cout << "\tpinch ae\t" << pinch[antimatter][e];
  cout << "\tpinch x\t"  << pinch[matter    ][mu];
  cout << "\tpinch ax\t" << pinch[antimatter][mu];
  cout << endl;
}


//=======//
// Cflux //
//=======//
double Cflux(double r){
  return 1. - sqrt(1.+Rnu/r)*sqrt(1.-Rnu/r);
}

#endif





