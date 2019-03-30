#ifndef PARAMETERS_H
#define PARAMETERS_H

const int NM=2;
enum state { matter,antimatter};
state operator++(state &n,int){ state tmp=n; n=(state)( (int)n+1 ); return tmp;};

const int NF=2;
enum flavour { e, mu };
flavour operator++(flavour &n,int){ flavour tmp=n; n=(flavour)( (int)n+1 ); return tmp;};

// number of parametrs needed to describe neutrino S matrix
const int NY=6; 

// number of nuclei followed 
const int NN=1; 

// number of energy bins, min and max energy
const int NE = 16;
double Emin, Emax;
vector<double> E; 

// the problem is broken up into seperate 'solutions'
const int NS=2; 
enum solution { msw, si};
solution operator++(solution &n,int){ solution tmp=n; n=(solution)( (int)n+1 ); return tmp;};

/* // vacuum Hamiltonian and mixing matrices */
/* vector<vector<MATRIX<complex<double>,NF,NF> > > HfV(NM); */
/* vector<MATRIX<complex<double>,NF,NF> > UV(NM); */

/* // vacuum values of the off-diagonal elements of the cofactor matrices */
/* vector<vector<MATRIX<complex<double>,NF,NF> > > CV; */

/* // mixing matrix element prefactors */
/* vector<vector<vector<double> > > AV; */

/* // initial mixing matrices, needs to be global */
/* vector<vector<MATRIX<complex<double>,NF,NF> > > U0(NM);  */

// neutron star radius
double Rnu;

// time snapshot
double t;

// vectors of mean energies, luminosities, temperatures, pinch paramaters
vector<vector<double> > meanE(NM,vector<double>(NF));
vector<vector<double> > L(NM,vector<double>(NF));
vector<vector<double> > pinch(NM,vector<double>(NF));

// pointers to initial spectra functions in found flux.h
vector<vector<double(*)(double)> > F0(NM,vector<double(*)(double)>(NF));
//

const double M_2PI = 2.*M_PI;
const complex<double> I = 1i;
// units, etc
namespace cgs{
  namespace units{
    const double cm = 1.; //cm
    const double eV = 1.60218e-12; //erg
  }
  namespace constants{
    const double hbar = 1.05457266e-27; // erg s
    const double c = 2.99792458e10; //cm/s
    const double c2 = c*c;
    const double c4 = c2*c2;
    const double hbarc = hbar*c;
    const double GF = 1.1663787e-5/*GeV^-2*//(1e9*1e9*units::eV*units::eV) * hbarc*hbarc*hbarc; //erg cm^3
    const double Mp = 1.6726219e-24; // g
  }
}

// mass of mass state1, delta m^2 differences
const double eV = 1.60218e-12; //erg
const double c = 2.99792458e10; // cm/s
const double c2 = c*c;
const double c4 = c2*c2;
const double m1 = 0.; // g, mass of state 1
const double dm21 = 2.43e-3/*eV^2*/ *eV*eV/c4; // g^2, delta m2^2-m1^2
const double a1 = dm21>0 ? -1. : +1.; // sign of eigenvalue 1 of vacuum Hamiltonian
const double a2 = dm21>0 ? +1. : -1.; // sign of eigenvalue 2 of vacuum Hamiltonian

// mixing angles
const double theta12V = 9. * M_PI/180.; // rad
const double c12V = cos(theta12V);
const double s12V = sin(theta12V);
const double alphaV[NF] = {0,0};
const double betaV[NF-1] = {0};
const double sin2thetaW = 0.23122;



#endif
