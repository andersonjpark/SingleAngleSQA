const int Nstates=2;
enum state { matter,antimatter};
state operator++(state &n,int){ state tmp=n; n=(state)( (int)n+1 ); return tmp;};

const int Nflavours=2;
enum flavour { e,mu};
flavour operator++(flavour &n,int){ flavour tmp=n; n=(flavour)( (int)n+1 ); return tmp;};

// vacuum eigenvalues
double m1,dm21;

//yzhu Nov23 int NE;
double Emin,Emax;
vector<double> E;
vector<double> matterFlux;
vector<double> antiMatterFlux;
vector<double> heavyFlux;

vector<vector<double> > kV;

vector<vector<MATRIX<complex<double> > > > HfV(Nstates);
MATRIX<complex<double> > Vf0(Nflavours,Nflavours);
vector<MATRIX<complex<double> > > U0;
vector<MATRIX<complex<double> > > Ubar0;

// parameters: square mass splitting, vacuum angle, momentum
double thetaV, alpha1V, alpha2V, betaV; 
// cos and sin of the vacuum angles
double cV,sV;


double Rnu;
double t,tscale;
vector<vector<double> > meanE(Nstates,vector<double>(Nflavours));
vector<vector<double> > T(Nstates,vector<double>(Nflavours));
vector<vector<double> > L(Nstates,vector<double>(Nflavours));
vector<vector<double> > alpha(Nstates,vector<double>(Nflavours));
vector<vector<double> > eta(Nstates,vector<double>(Nflavours));
vector<vector<double> > F2(Nstates,vector<double>(Nflavours)), F3(Nstates,vector<double>(Nflavours));

vector<vector<double(*)(double)> > F0(Nstates,vector<double(*)(double)>(Nflavours));

enum eigenvalue_ordering { normal,inverted } hierarchy;

int a1,a2;



