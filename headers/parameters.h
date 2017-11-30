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

// mass of mass state1, delta m^2 differences
double m1,dm21;

// number of energy bins, min and max energy
int NE;
double Emin, Emax;
vector<double> E; 

// the problem is broken up into seperate 'solutions'
const int NS=2; 
enum solution { msw, si};
solution operator++(solution &n,int){ solution tmp=n; n=(solution)( (int)n+1 ); return tmp;};

// vacuum eigenvalues
vector<vector<double> > kV;
vector<int> ordering(NF);
int a1,a2;

// vacuum Hamiltonian and mixing matrices
vector<vector<MATRIX<complex<double>,NF,NF> > > HfV(NM);
vector<MATRIX<complex<double>,NF,NF> > UV(NM);

// vacuum values of the off-diagonal elements of the cofactor matrices
vector<vector<MATRIX<complex<double>,NF,NF> > > CV;

// mixing matrix element prefactors
vector<vector<vector<double> > > AV;

// initial mixing matrices, needs to be global
vector<vector<MATRIX<complex<double>,NF,NF> > > U0(NM); 

double theta12V;
vector<double> alphaV(NF), betaV(NF-1);
double c12V,s12V;

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



