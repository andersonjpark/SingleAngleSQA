#if !defined(_FLUX)
#define _FLUX

void Initialize_Luminosities(string spectrapath, double t);
void Initialize_MeanEnergies(string spectrapath, double t);
void Initialize_PinchParameters(string spectrapath, double t);

double Fnue0(double E);
double Fnumu0(double E);
double Fnutau0(double E);

double Fanue0(double E);
double Fanumu0(double E);
double Fanutau0(double E);

double Cflux(double r);

// *************************************************************************

void Initialize_Luminosities(string spectrapath, double t)
     { DISCONTINUOUS Lnue,Lanue,Lnux,Lanux;

       Lnue.Open((spectrapath  + string("L.nue.dat")).c_str());
       Lanue.Open((spectrapath + string("L.anue.dat")).c_str());
       Lnux.Open((spectrapath  + string("L.nux.dat")).c_str());
       Lanux.Open((spectrapath + string("L.anux.dat")).c_str());
  
       L[matter][e]       = Lnue(t);  // ergs/s      
       L[antimatter][e]   = Lanue(t); 
       L[matter][mu]      = Lnux(t);        
       L[antimatter][mu]  = Lanux(t); 
//       L[matter][tau]     = L[matter][mu]; 
  //     L[antimatter][tau] = L[antimatter][mu];

       //L[matter][e]=0.;         L[antimatter][e]=0.; // ergs/s
       //L[matter][mu]=0.;        L[antimatter][mu]=0.; 
       //L[matter][tau]=0.;       L[antimatter][tau]=0.;
      
       cout << endl << "Le\t" << L[matter][e] << "\tLae\t" << L[antimatter][e] << "\tLnux\t"  << L[matter][mu] << "\tLanux\t" << L[antimatter][mu]  << " erg/s." << endl;     
      }

// ****************************************************

void Initialize_MeanEnergies(string spectrapath, double t)
     { DISCONTINUOUS meanEnue,meanEanue,meanEnux,meanEanux;

       meanEnue.Open((spectrapath  + string("meanE.nue.dat")).c_str());
       meanEanue.Open((spectrapath + string("meanE.anue.dat")).c_str());
       meanEnux.Open((spectrapath  + string("meanE.nux.dat")).c_str());
       meanEanux.Open((spectrapath + string("meanE.anux.dat")).c_str());
	    
       meanE[matter][e]       = meanEnue(t);         
       meanE[antimatter][e]   = meanEanue(t); // MeV
       meanE[matter][mu]      = meanEnux(t);       
       meanE[antimatter][mu]  = meanEanux(t);  
    //   meanE[matter][tau]     = meanE[matter][mu]; 
    //   meanE[antimatter][tau] = meanE[antimatter][mu];

       cout << endl << "Ee\t" << meanE[matter][e] << "\tEae\t" << meanE[antimatter][e] << "\tEnux\t"  << meanE[matter][mu] << "\tEanux\t" << meanE[antimatter][mu] << " MeV." << endl; 

       for(state m=matter;m<=antimatter;m++){ for(flavour f=e;f<=mu;f++){ meanE[m][f] *= mega*cgs::units::eV;} }
      }

// *************************************************************************

void Initialize_PinchParameters(string spectrapath, double t)
     { DISCONTINUOUS rmsEnue,rmsEanue,rmsEnux,rmsEanux;
  
       rmsEnue.Open((spectrapath  + string("rmsE.nue.dat")).c_str());
       rmsEanue.Open((spectrapath + string("rmsE.anue.dat")).c_str());
       rmsEnux.Open((spectrapath  + string("rmsE.nux.dat")).c_str());
       rmsEanux.Open((spectrapath + string("rmsE.anux.dat")).c_str());

       double CE;
       CE = rmsEnue(t)*mega*cgs::units::eV/meanE[matter][e] ;
       pinch[matter][e] = (2.-CE*CE)/(CE*CE-1.);
  
       CE = rmsEanue(t)*mega*cgs::units::eV/meanE[antimatter][e] ;
       pinch[antimatter][e] = (2.-CE*CE)/(CE*CE-1.);
       
       CE = rmsEnux(t)*mega*cgs::units::eV/meanE[matter][mu];
       pinch[matter][mu]  = (2.-CE*CE)/(CE*CE-1.); 
//       pinch[matter][tau] = pinch[matter][mu];
  
       CE = rmsEanux(t)*mega*cgs::units::eV/meanE[antimatter][mu];
       pinch[antimatter][mu]  = (2.-CE*CE)/(CE*CE-1.); 
  //     pinch[antimatter][tau] = pinch[antimatter][mu];
  
       cout << endl << "rmsE\t" << rmsEnue(t) << "\trmsEae\t" << rmsEanue(t) << "\trmsEnux\t" << rmsEnux(t) << "\trmsEanux\t" << rmsEanux(t) << " eV." << endl;

       cout << "pinch e\t" << pinch[matter][e] << "\tpinch ae\t" << pinch[antimatter][e] << "\tpinch x\t" << pinch[matter][mu] << "\tpinch ax\t" << pinch[antimatter][mu] << endl;
      }

// ***************************************************************************

// differential spectrum in units of /erg/s
double Fnue0(double E)
       { return pow(pinch[matter][e]+1.,pinch[matter][e]+1.)/meanE[matter][e]/Gamma(pinch[matter][e]+1.)*pow(E/meanE[matter][e],pinch[matter][e])*exp(-(pinch[matter][e]+1.)*E/meanE[matter][e]);}

double Fanue0(double E)
       { return pow(pinch[antimatter][e]+1.,pinch[antimatter][e]+1.)/meanE[antimatter][e]/Gamma(pinch[antimatter][e]+1.)*pow(E/meanE[antimatter][e],pinch[antimatter][e])*exp(-(pinch[antimatter][e]+1.)*E/meanE[antimatter][e]);}

double Fnumu0(double E)
       { return pow(pinch[matter][mu]+1.,pinch[matter][mu]+1.)/meanE[matter][mu]/Gamma(pinch[matter][mu]+1.)*pow(E/meanE[matter][mu],pinch[matter][mu])*exp(-(pinch[matter][mu]+1.)*E/meanE[matter][mu]);}

double Fanumu0(double E)
       { return pow(pinch[antimatter][mu]+1.,pinch[antimatter][mu]+1.)/meanE[antimatter][mu]/Gamma(pinch[antimatter][mu]+1.)*pow(E/meanE[antimatter][mu],pinch[antimatter][mu])*exp(-(pinch[antimatter][mu]+1.)*E/meanE[antimatter][mu]);}
/*
double Fnutau0(double E)
       { return pow(pinch[matter][tau]+1.,pinch[matter][tau]+1.)/meanE[matter][tau]/Gamma(pinch[matter][tau]+1.)*pow(E/meanE[matter][tau],pinch[matter][tau])*exp(-(pinch[matter][tau]+1.)*E/meanE[matter][tau]);}

double Fanutau0(double E)
       { return pow(pinch[antimatter][tau]+1.,pinch[antimatter][tau]+1.)/meanE[antimatter][tau]/Gamma(pinch[antimatter][tau]+1.)*pow(E/meanE[antimatter][tau],pinch[antimatter][tau])*exp(-(pinch[antimatter][tau]+1.)*E/meanE[antimatter][tau]);}
*/
// ***************************************************************************

double Cflux(double r){ return 1. - sqrt(1.+Rnu/r)*sqrt(1.-Rnu/r);}

#endif





