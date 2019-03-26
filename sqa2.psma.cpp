
#include <complex>
using std::complex;
using std::polar;
using std::abs;
using std::arg;
using std::real;
using std::imag;
#include <cstdarg>
using std::va_list;
#include <cstdlib>

#include<iostream>
using std::cout;
#include<ostream>
using std::ostream;
using std::endl;
using std::flush;
#include<fstream>
using std::ifstream;
using std::ofstream;
#include<sstream>
using std::stringstream;

#include<algorithm>
using std::min;
using std::max;
using std::swap;
using std::lower_bound;
#include<string>
using std::string;
#include <utility>
using std::pair;

#include<functional>
#include<limits>
using std::numeric_limits;

#include<vector>
using std::vector;

#include "mstl.h"
using interpolation::DISCONTINUOUS;

// global variables
DISCONTINUOUS rho, lnrho, Ye, temperature; // rho is the mass density
bool do_interact;

// headers
#include "headers/parameters.h"
#include "headers/potentials.h"
#include "headers/single_angle.h"
#include "headers/flavour_basis.h"
#include "headers/eigenvalues.h"
#include "headers/mixing_angles.h"
#include "headers/adiabatic_basis.h"
#include "headers/flux.h"
#include "headers/jacobians.h"
#include "headers/multiEnergy.h"
#include "headers/MNR.h"

//vector<vector<MATRIX<complex<double>,NF,NF> > > rhomatrixf0(NM), rhomatrixm0(NM);
vector<vector<MATRIX<complex<double>,NF,NF> > > pmatrixf0(NM), pmatrixm0(NM);
vector<vector<MATRIX<complex<double>,NF,NF> > > fmatrixf(NM), fmatrixm(NM);
vector<DISCONTINUOUS> eP,eBarP,xP;
vector<DISCONTINUOUS> eD,eBarD,xD;
vector<vector<double> > dphi_dr_interact, dtheta_dr_interact;
vector<vector<double> > dphi_dr_osc,      dtheta_dr_osc;

MATRIX<complex<double>,NF,NF> B(vector<double> y);
void K(double r,
       double dr,
       vector<vector<vector<vector<double> > > > &Y,
       vector<vector<vector<MATRIX<complex<double>,NF,NF> > > > &C0,
       vector<vector<vector<vector<double> > > > &A0,
       vector<vector<vector<vector<double> > > > &K);
void Outputvsr(ofstream &fout,
	       ofstream &foutP,
	       ofstream &foutf,
	       ofstream &foutdangledr,	       
	       ofstream *foutPvsr,
	       ofstream *foutFvsr,
	       double r,
	       vector<vector<vector<vector<double> > > > Y,
	       vector<vector<vector<MATRIX<complex<double>,NF,NF> > > > C0,
	       vector<vector<vector<vector<double> > > > A0,
	       vector<vector<MATRIX<complex<double>,NF,NF> > > Scumulative);

#include "headers/update.h"
#include "headers/project/albino.h"
//#include "headers/project/test_case_B.h"

void getP(const double r,
	  const vector<vector<MATRIX<complex<double>,NF,NF> > > U0, 
	  const vector<vector<MATRIX<complex<double>,NF,NF> > > Scumulative, 
	  vector<vector<MATRIX<complex<double>,NF,NF> > >& pmatrixf0,
	  vector<vector<MATRIX<complex<double>,NF,NF> > >& pmatrixm0){

  for(int i=0;i<=NE-1;i++){
    for(state m=matter; m<=antimatter; m++){
      // get unoscillated potential
      getPunosc(r, m, i, pmatrixf0[m][i]);

      // oscillate the potential and put into the mass basis
      pmatrixm0[m][i] = Scumulative[m][i]
	* Adjoint(U0[m][i])
	* pmatrixf0[m][i]
	* U0[m][i]
	* Adjoint(Scumulative[m][i]);
      pmatrixf0[m][i] =  U0[m][i] * pmatrixm0[m][i] * Adjoint(U0[m][i]);
    }
  }
}

void interact(vector<vector<MATRIX<complex<double>,NF,NF> > >& fmatrixf,
	      vector<vector<MATRIX<complex<double>,NF,NF> > >& Scumulative,
	      vector<vector<MATRIX<complex<double>,NF,NF> > >& U0,
	      double rho, double T, double Ye, double r, double dr){
  // save old fmatrix
  vector<vector<MATRIX<complex<double>,NF,NF> > > old_fmatrixf = fmatrixf;

  // let neutrinos interact
  if(do_interact) my_interact(fmatrixf, Scumulative, rho, T, Ye, r, dr);
  for(int i=0; i<NE; i++){
    for(state m=matter; m<=antimatter; m++){
      fmatrixm[m][i] = Adjoint(U0[m][i]) * fmatrixf[m][i] * U0[m][i];
    }
  }
  
  // loop through getting rotation matrices
  double hold[4];
  double hnew[4];
  for(int i=0; i<NE; i++){
    for(state m=matter; m<=antimatter; m++){
      pauli_decompose(old_fmatrixf[m][i], hold);
      pauli_decompose(    fmatrixf[m][i], hnew);

      // get the theta and phi contribution
      double oldmag2   = hold[0]*hold[0] + hold[1]*hold[1] + hold[2]*hold[2];
      double newmag2   = hnew[0]*hnew[0] + hnew[1]*hnew[1] + hnew[2]*hnew[2];
      if(oldmag2==0 or newmag2==0){
	continue;
      }
      dtheta_dr_interact[i][m] = (acos(hnew[2]/sqrt(newmag2)) - acos(hold[2]/sqrt(oldmag2))) / dr;
      dphi_dr_interact[i][m] = (atan2(hnew[1],hnew[0]) - atan2(hold[1],hold[0])) / dr;

      // get the axis of rotation
      double lrot[3];
      lrot[0] =   hold[1]*hnew[2] - hold[2]*hnew[1] ;
      lrot[1] = -(hold[0]*hnew[2] - hold[2]*hnew[0]);
      lrot[2] =   hold[0]*hnew[1] - hold[1]*hnew[0] ;
      double lmag2 = lrot[0]*lrot[0] + lrot[1]*lrot[1] + lrot[2]*lrot[2];
      double lmag = sqrt(lmag2);
      double sinalpha = sqrt(lmag2 / (oldmag2*newmag2));
      if(lmag > 0)
	for(unsigned i=0; i<3; i++) lrot[i] /= sqrt(lmag);
      else continue;

      // get the rotation operator in the flavor basis
      double alpha = asin(sinalpha);
      complex<double> Rcoeff[4];
      for(int i=0; i<3; i++) Rcoeff[i] = -I * sin(alpha/2.) * lrot[i];
      Rcoeff[3] = cos(alpha/2.);
      
      MATRIX<complex<double>,NF,NF> R;
      pauli_reconstruct(Rcoeff, R);

      // apply to Scumulative
      Scumulative[m][i] = MATRIX<complex<double>,NF,NF>(Adjoint(U0[m][i]) * R * U0[m][i] * Scumulative[m][i]);
    }
  }
}

#include "headers/nulib_interface.h"

//======//
// MAIN //
//======//
int main(int argc, char *argv[]){
  try{ 
    int in=1;
    string inputfilename,potential_directory;
    ofstream fout,foutC,foutP,foutS, foutf, foutdangledr;
    string outputfilename,rhofilename, Yefilename, vfilename, spectrapath, nulibfilename, temperaturefilename;
    string outputfilenamestem;
    string note;
    
    inputfilename=string(argv[in]);
    ifstream fin(inputfilename.c_str());
    
    // load the nulib table
    fin>>nulibfilename;
    fin>>potential_directory;
    cout << nulibfilename << endl;
    nulib_init(nulibfilename, 0);

    fin>>rhofilename;
    fin>>Yefilename;
    fin>>temperaturefilename;
    //fin>>spectrapath;
    fin>>outputfilename;
    outputfilenamestem = outputfilename+"/";
    
    double rmin, rmax;
    fin>>rmin>>rmax; // cm
    
    fin>>NE; //>>Emin>>Emax; // MeV
    
    //fin>>Rnu; // cm
    //fin>>t;   // s
    
    m1=0.;
    fin>>dm21;
    fin>>theta12V;
    
    alphaV[0]=0.;
    alphaV[1]=0.;
    
    betaV[0]=0.;
    
    double accuracy;
    fin>>accuracy;
    fin>>note;
    int iNT,fNT;
    string id;
    // fin>>iNT;
    // fin>>fNT;
    // fin>>id;//number of tracer

    double artificial_scaling;
    fin >> artificial_scaling;
    rmin *= artificial_scaling;
    rmax *= artificial_scaling;
    fin >> do_interact;
    
    cout<<"\n\n*********************************************************\n";
    cout<<"\nrho\t"<<rhofilename;
    cout<<"\nYe\t"<<Yefilename;
    cout<<"\nT\t"<<temperaturefilename;
    cout<<"\noutput\t"<<outputfilename;
    cout<<"\nrmin\t"<<rmin<<"\trmax\t"<<rmax;
    cout<<"\nartificial_scaling\t"<<artificial_scaling;
    cout<<"\ndo_interact\t"<<do_interact;
    //cout<<"\nRnu\t"<<Rnu<<"\nt\t"<<t;
    
    //cout << "\n\nNE\t" << NE << "\tEmin\t" << Emin << "\tEmax\t" << Emax;
    
    cout<<"\n\nm1\t"<<m1<<"\tdm21^2\t"<<dm21;
    cout<<"\ntheta12V\t"<<theta12V;
    cout<<"\nalpha1V\t"<<alphaV[0]<<"\talpha2V\t"<<alphaV[1];
    cout<<"\nbeta1V\t"<<betaV[0];
    
    cout<<"\naccuracy\t"<<accuracy<<"\n";
    cout.flush();
    
    // load rho and Ye data
    rho.Open(rhofilename,'#');
    Ye.Open(Yefilename,'#');
    temperature.Open(temperaturefilename,'#');
    vector<double> newX = rho.X();
    for(unsigned i=0; i<newX.size(); i++)
      newX[i] *= artificial_scaling;
    rho.SetX(newX);
    Ye.SetX(newX);
    temperature.SetX(newX);
    rmin=max(rmin,max(rho.XMin(),Ye.XMin()) );
    rmin = max(rmin, temperature.XMin() );
    rmax=min(rmax,min(rho.XMax(),Ye.XMax()) );
    rmax = min(rmax, temperature.XMax() );

    lnrho=rho;
    lnrho.TransformX(log);
    lnrho.TransformY(log);
    
    // load and compute spectral data    
    eP.resize(NE);
    eBarP.resize(NE);
    xP.resize(NE);
    eD.resize(NE);
    eBarD.resize(NE);
    xD.resize(NE);
    dtheta_dr_interact = vector< vector<double> > (NE, vector<double>(2,0) );
    dphi_dr_interact = vector< vector<double> > (NE, vector<double>(2,0) );
    dtheta_dr_osc      = vector< vector<double> > (NE, vector<double>(2,0) );
    dphi_dr_osc      = vector< vector<double> > (NE, vector<double>(2,0) );
    
    // load and compute spectral data
    for(int i=0;i<=NE-1;i++){
      eP   [i].Open(potential_directory+"/potential_s1_g"+patch::to_string(i+1)+note+".txt",'#');
      eBarP[i].Open(potential_directory+"/potential_s2_g"+patch::to_string(i+1)+note+".txt",'#');
      xP   [i].Open(potential_directory+"/potential_s3_g"+patch::to_string(i+1)+note+".txt",'#');
      eD   [i].Open(potential_directory+"/density_s1_g"+patch::to_string(i+1)+note+".txt",'#');
      eBarD[i].Open(potential_directory+"/density_s2_g"+patch::to_string(i+1)+note+".txt",'#');
      xD   [i].Open(potential_directory+"/density_s3_g"+patch::to_string(i+1)+note+".txt",'#');
      newX = eP[i].X();
      for(unsigned i=0; i<newX.size(); i++)
	newX[i] *= artificial_scaling;
      eP[i].SetX(newX);
      eBarP[i].SetX(newX);
      xP[i].SetX(newX);
      eD[i].SetX(newX);
      eBarD[i].SetX(newX);
      xD[i].SetX(newX);
    }

    // output filestreams: the arrays of ofstreams cannot use the vector container - bug in g++
    foutS.open((outputfilename+"/S.dat").c_str());
    foutS.precision(12);
    foutS.flush();
    foutP.open((outputfilename+"/2p.dat").c_str());
    foutP.precision(12);
    foutP.flush();
    fout.open((outputfilename+"/out.dat").c_str());
    fout.precision(12);
    foutf.open((outputfilename+"/f.dat").c_str());
    foutf.precision(12);
    foutf << "# 1:r ";
    fout << "# 1:r ";
    for(int i=0; i<NE; i++)
      for(state m=matter; m<=antimatter; m++)
	for(flavour f1=e; f1<=mu; f1++)
	  for(flavour f2=e; f2<=mu; f2++) {
	    int istart = 2*( f2 + f1*2 + m*2*2 + i*2*2*2) + 2;
	    foutf << istart   << ":ie"<<i<<"m"<<m<<"f"<<f1<<f2<<"R\t";
	    foutf << istart+1 << ":ie"<<i<<"m"<<m<<"f"<<f1<<f2<<"I\t";
	    fout  << istart   << ":ie"<<i<<"m"<<m<<"f"<<f1<<f2<<"R\t";
	    fout  << istart+1 << ":ie"<<i<<"m"<<m<<"f"<<f1<<f2<<"I\t";
    }
    foutf << endl;
    fout << endl;
    fout.flush();
    foutf.flush();
    
    foutdangledr.open((outputfilename+"/dangledr.dat").c_str());
    foutdangledr.precision(12);
    foutdangledr << "# 1:r ";
    for(state m=matter; m<=antimatter; m++)
      for(int i=0; i<NE; i++)
	foutdangledr << 2+0*NE*2 + i + m*NE << ":OscThetaie"<<i<<"m"<<m<<"\t";
    for(state m=matter; m<=antimatter; m++)
      for(int i=0; i<NE; i++)
	foutdangledr << 2+1*NE*2 + i + m*NE << ":OscPhiie"<<i<<"m"<<m<<"\t";
    for(state m=matter; m<=antimatter; m++)
      for(int i=0; i<NE; i++)
	foutdangledr << 2+2*NE*2 + i + m*NE << ":InteractThetaie"<<i<<"m"<<m<<"\t";
    for(state m=matter; m<=antimatter; m++)
      for(int i=0; i<NE; i++)
	foutdangledr << 2+3*NE*2 + i + m*NE << ":InteractPhiie"<<i<<"m"<<m<<"\t";
    foutdangledr << endl;
    foutdangledr.flush();
    
    ofstream foutPvsr[NE],foutFvsr[NE];
    ofstream fPvsE, fFvsE;
    
    for(int i=0;i<=NE-1;i++){
      foutPvsr[i].open(outputfilename+"/g"+std::to_string(i)+"_Pvsr.dat");
      foutPvsr[i].precision(12);
      
      foutFvsr[i].open(outputfilename+"/g"+std::to_string(i)+"_Fvsr.dat");
      foutFvsr[i].precision(12);
    }
    
    // unit conversion to cgs
    //Emin *= 1.*mega*cgs::units::eV;
    //Emax *= 1.*mega*cgs::units::eV;
    m1   *= 1.*cgs::units::eV/cgs::constants::c2;
    dm21 *= 1.*cgs::units::eV*cgs::units::eV/cgs::constants::c4;
    theta12V *= M_PI/180.;
    c12V = cos(theta12V);
    s12V = sin(theta12V);
    
    // *************************************************
    // set up global variables defined in parameters.h *
    // *************************************************

    // vectors of energies and vacuum eigenvalues
    kV = vector<vector<double> >(NE,vector<double>(NF));
    set_Ebins(E);
    set_kV(kV);

    // determine eigenvalue ordering
    ordering[0]=0; ordering[1]=1;
    vector<double> tempkV(kV[0]);
    Sort(tempkV,ordering,ascending);
    
    if(kV[0][1]>kV[0][0]){
      a1=-1;
      a2=+1;
      cout<<"\n\nNormal hierarchy";
    }
    else{ 
      if(kV[0][1]<kV[0][0]){
	a1=+1; 
	a2=-1; 
	cout<<"\n\nInverted hierarchy";
      }
      else{ 
	cout<<endl<<endl<<"Neither normal or Inverted"<<endl; 
	abort();
      }
    }
    
    // vaccum mixing matrices and Hamiltonians
    Evaluate_UV();
    
    HfV[matter] = vector<MATRIX<complex<double>,NF,NF> >(NE);
    HfV[antimatter] = vector<MATRIX<complex<double>,NF,NF> >(NE);
    Evaluate_HfV();
    
    // cofactor matrices in vacuum
    CV=vector<vector<MATRIX<complex<double>,NF,NF> > >(NE,vector<MATRIX<complex<double>,NF,NF> >(NF));
    Evaluate_CV();
    
    // mixing matrix element prefactors in vacuum
    AV=vector<vector<vector<double> > >(NE,vector<vector<double> >(NF,vector<double>(NF)));
    Evaluate_AV();
    
    // **************************************
    // quantities evaluated at inital point *
    // **************************************
    
    // MSW potential matrix
    double rrho = get_rho(rmin);
    double YYe  = get_Ye(rmin);
    
    MATRIX<complex<double>,NF,NF> VfMSW0, Hf0;
    vector<double> k0, deltak0;
    
    VfMSW0[e][e]=Ve(rrho,YYe);
    VfMSW0[mu][mu]=Vmu(rrho,YYe);
    
    // cofactor matrices at initial point - will be recycled as cofactor matrices at beginning of every step
    vector<vector<vector<MATRIX<complex<double>,NF,NF> > > > 
      C0(NM,vector<vector<MATRIX<complex<double>,NF,NF> > >(NE,vector<MATRIX<complex<double>,NF,NF> >(NF)));

    // mixing matrix element prefactors at initial point - will be recycled like C0
    vector<vector<vector<vector<double> > > > 
      A0(NM,vector<vector<vector<double> > >(NE,vector<vector<double> >(NF,vector<double>(NF))));
    
    // accumulated S matrices from prior integration domains
    vector<vector<MATRIX<complex<double>,NF,NF> > > 
      Scumulative(NM,vector<MATRIX<complex<double>,NF,NF> >(NE,UnitMatrix<complex<double> >(NF)));

    // mixing angles to MSW basis at initial point
    U0[matter] = vector<MATRIX<complex<double>,NF,NF> >(NE);
    U0[antimatter] = vector<MATRIX<complex<double>,NF,NF> >(NE);
    
    for(int i=0;i<=NE-1;i++){
      Hf0=HfV[matter][i]+VfMSW0;
      k0=k(Hf0);
      deltak0=deltak(Hf0);
      C0[matter][i]=CofactorMatrices(Hf0,k0);
      
      for(int j=0;j<=NF-1;j++){
	if(real(C0[matter][i][j][mu][e]*CV[i][j][mu][e]) < 0.)
	  A0[matter][i][j][e]=-AV[i][j][e];
	else A0[matter][i][j][e]=AV[i][j][e];
	A0[matter][i][j][mu]=AV[i][j][mu];
      }
      U0[matter][i]=U(deltak0,C0[matter][i],A0[matter][i]);
      
      Hf0=HfV[antimatter][i]-VfMSW0;
      k0=kbar(Hf0);
      deltak0=deltakbar(Hf0);
      C0[antimatter][i]=CofactorMatrices(Hf0,k0);
      for(int j=0;j<=NF-1;j++){
	if(real(C0[antimatter][i][j][mu][e]*CV[i][j][mu][e]) < 0.)
	  A0[antimatter][i][j][e]=-AV[i][j][e];
	else A0[antimatter][i][j][e]=AV[i][j][e];
	A0[antimatter][i][j][mu]=AV[i][j][mu];
      }
      U0[antimatter][i]=Conjugate(U(deltak0,C0[antimatter][i],A0[antimatter][i]));
    }
    
    vector<vector<vector<double> > > P0 (NM,vector<vector<double> >(NF,vector <double>(NE)));

    // density matrices at initial point, rhomatrixm0 - but not rhomatrixf0
    // will be updated whenever discontinuities are crossed and/or S is reset
    pmatrixf0[matter]=pmatrixf0[antimatter]=vector<MATRIX<complex<double>,NF,NF> >(NE);
    pmatrixm0[matter]=pmatrixm0[antimatter]=vector<MATRIX<complex<double>,NF,NF> >(NE);
    fmatrixf[matter]=fmatrixf[antimatter]=vector<MATRIX<complex<double>,NF,NF> >(NE);
    fmatrixm[matter]=fmatrixm[antimatter]=vector<MATRIX<complex<double>,NF,NF> >(NE);

    // yzhu14 density/potential matrices art rmin
    double rho0 = rho(rmin);
    double T0 = temperature(rmin);
    double ye0 = Ye(rmin);
    initialize(fmatrixf,rmin,rho0,T0,ye0);
    getP(rmin,U0,Scumulative,pmatrixf0,pmatrixm0);

    // ***************************************
    // quantities needed for the calculation *
    // ***************************************
    double r,r0,dr,drmin,dr_this_step;
    int ND;
    vector<double> rs;
    
    double maxerror,increase=3.;
    bool repeat, finish, resetflag, output;
    int counterout,step;
    
    // comment out if not following as a function of r
    fin>>step;
    
    // self-interaction integration factor
    //NSI=M_SQRT2*cgs::constants::GF*(Emax-Emin)/(NE-1.);
    
    // ***************************************
    // variables followed as a function of r *
    // ***************************************
    
    vector<vector<vector<vector<double> > > > 
      Y(NM,vector<vector<vector<double> > >(NE,vector<vector<double> >(NS,vector<double>(NY))));
    vector<vector<vector<vector<double> > > > 
      dY(NM,vector<vector<vector<double> > >(NE,vector<vector<double> >(NS,vector<double>(NY))));
    vector<vector<vector<vector<double> > > > 
      Y0(NM,vector<vector<vector<double> > >(NE,vector<vector<double> >(NS,vector<double>(NY))));
    vector<vector<vector<vector<double> > > > 
      Yerror(NM,vector<vector<vector<double> > >(NE,vector<vector<double> >(NS,vector<double>(NY))));
    
    // cofactor matrices
    vector<vector<vector<MATRIX<complex<double>,NF,NF> > > > C=C0;;
    // mixing matrix prefactors
    vector<vector<vector<vector<double> > > > A=A0;
        
    // ************************
    // Runge-Kutta quantities *
    // ************************
    int NRK,NRKOrder;
    const double *AA=NULL,**BB=NULL,*CC=NULL,*DD=NULL;
    RungeKuttaCashKarpParameters(NRK,NRKOrder,AA,BB,CC,DD);
    
    vector<vector<vector<vector<vector<double> > > > > 
      Ks(NRK,vector<vector<vector<vector<double> > > >
	 (NM,vector<vector<vector<double> > >(NE,vector<vector<double> >(NS,vector<double>(NY)))));
    
    // temporaries
    MATRIX<complex<double>,NF,NF> SSMSW,SSSI,SThisStep;
    
    // *********************
    // integration domians *
    // *********************
    
    rho.FindDomains();
    
    rs.push_back(rmin);
    for(int d=1;d<=rho.NDomains();d++){
      r=rho.Discontinuity(d);
      if(r>rmin && r<rmax) rs.push_back(r);
    }
    rs.push_back(rmax);
    
    sort(rs.begin(),rs.end());
    ND=rs.size()-1;
    
    // **********************
    // start of calculation *
    // **********************
    
    for(int d=0;d<=ND-1;d++){
      if(d==0) rmin=rs.front();
      else rmin=rs[d]+1.*cgs::units::cm;
      
      if(d==ND-1) rmax=rs.back();
      else rmax=rs[d+1]-1.*cgs::units::cm;
      
      cout<<"\n"<<d<<"\t"<<rmin<<"\t"<<rmax << endl;
      cout << endl;
      cout << "r(km)  dr(cm) rho(g/ccm) T(MeV) Ye" << endl;
      cout.flush();
      
      // *****************************************
      // initialize at beginning of every domain *
      // *****************************************
      r=rmin;
      dr=1e-3*cgs::units::cm;
      drmin=4.*r*numeric_limits<double>::epsilon();
      
      for(state m=matter;m<=antimatter;m++)
	for(int i=0;i<=NE-1;i++){
	  Y[m][i][msw][0] = Y[m][i][si][0] = M_PI/2.;
	  Y[m][i][msw][1] = Y[m][i][si][1] = M_PI/2.;
	  Y[m][i][msw][2] = Y[m][i][si][2] = 0.;
	  Y[m][i][msw][3] = Y[m][i][si][3] = 1.; // The determinant of the S matrix
	  Y[m][i][msw][4] = Y[m][i][si][4] = 0.;
	  Y[m][i][msw][5] = Y[m][i][si][5] = 0.;
	}
      
      // *************************************************
      // comment out if not following as a function of r *
      // *************************************************
	
      finish=output=false;
      counterout=1;
      Outputvsr(fout,foutP,foutf,foutdangledr,foutPvsr,foutFvsr,r,Y,C,A,Scumulative);
	
      for(state m=matter; m<=antimatter; m++)
	for(int i=0; i<NE; i++)
	  for(flavour f1=e; f1<=mu; f1++)
	    for(flavour f2=e; f2<=mu; f2++)
	      assert(fmatrixf[m][i][f1][f2] == fmatrixf[m][i][f1][f2]);
      
      // ***********************
      // start the loop over r *
      // ***********************
      do{ 
	double intkm = int(r/1e5)*1e5;
	if(r - intkm <= dr){
	  cout << r/1e5 << " " << dr << " " << rho(r) << " " << temperature(r) << " " << Ye(r) << endl;
	  cout.flush();
	}

	if(r+dr>rmax){
	  dr=rmax-r;
	  finish=true;
	  output=true;
	}
	  
	r0=r;
	Y0=Y;
	C0=C;
	A0=A;
	  
	for(state m=matter; m<=antimatter; m++)
	  for(int i=0; i<NE; i++)
	    for(flavour f1=e; f1<=mu; f1++)
	      for(flavour f2=e; f2<=mu; f2++)
		assert(fmatrixf[m][i][f1][f2] == fmatrixf[m][i][f1][f2]);
	
	// beginning of RK section
	do{ 
	  repeat=false;
	  for(int k=0;k<=NRK-1;k++){
	    r=r0+AA[k]*dr;
	    Y=Y0;

	    for(state m = matter; m <= antimatter; m++)
	      for(int i=0;i<=NE-1;i++)
		for(solution x=msw;x<=si;x++)
		  for(int j=0;j<=NY-1;j++)
		    for(int l=0;l<=k-1;l++)
		      Y[m][i][x][j] += BB[k][l] * Ks[l][m][i][x][j];

	    getP(r,U0,Scumulative,pmatrixf0,pmatrixm0);
	    K(r,dr,Y,C,A,Ks[k]);
	  }
	  
	  // increment all quantities and update C and A arrays
	  r=r0+dr;
	  for(state m=matter;m<=antimatter;m++){
	    for(int i=0;i<=NE-1;i++){
	      for(solution x=msw;x<=si;x++){
		for(int j=0;j<=NY-1;j++){
		  Y[m][i][x][j] = Y0[m][i][x][j];
		  Yerror[m][i][x][j] = 0.;
		  for(int k=0;k<=NRK-1;k++){
		    assert(CC[k] == CC[k]);
		    assert(Ks[k][m][i][x][j] == Ks[k][m][i][x][j]);
		    Y[m][i][x][j] += CC[k] * Ks[k][m][i][x][j];
		    Yerror[m][i][x][j] += (CC[k]-DD[k]) * Ks[k][m][i][x][j];
		    assert(Y[m][i][x][j] == Y[m][i][x][j]);
		  }
		}
	      }
	    }
	  }
	  
	  C=UpdateC(r,get_Ye(r));
	  A=UpdateA(C,C0,A0);
	    
	  // find largest error
	  maxerror=0.;
	  for(state m=matter;m<=antimatter;m++)
	    for(int i=0;i<=NE-1;i++)
	      for(solution x=msw;x<=si;x++)
		for(int j=0;j<=NY-1;j++)
		  maxerror = max( maxerror, fabs(Yerror[m][i][x][j]) );
	    
	  // decide whether to accept step, if not adjust step size
	  dr_this_step = dr;
	  if(maxerror>accuracy){
	    dr *= 0.9 * pow(accuracy/maxerror, 1./(NRKOrder-1.));
	    if(dr > drmin) repeat=true;
	  }

	  // reset integration variables to those at beginning of step
	  if(repeat==true){
	    r=r0;
	    Y=Y0;
	    C=C0;
	    A=A0;
	  }
	  
	}while(repeat==true); // end of RK section

	// interact with the matter
	interact(fmatrixf, Scumulative, U0, rho(r), temperature(r), Ye(r), r, dr_this_step);
	for(state m=matter; m<=antimatter; m++)
	  for(int i=0; i<NE; i++)
	    for(flavour f1=e; f1<=mu; f1++)
	      for(flavour f2=e; f2<=mu; f2++)
		assert(fmatrixf[m][i][f1][f2] == fmatrixf[m][i][f1][f2]);

	// accumulate S and reset variables
	vector<vector<MATRIX<complex<double>,NF,NF> > > old_fmatrixf = fmatrixf;
	for(state m=matter;m<=antimatter;m++){
	  for(int i=0;i<=NE-1;i++){
	    SSMSW = W(Y[m][i][msw])*B(Y[m][i][msw]);
	    SSSI  = W(Y[m][i][si ])*B(Y[m][i][si ]);
	    SThisStep = SSMSW*SSSI;
	    Scumulative[m][i]=MATRIX<complex<double>,NF,NF>(SThisStep*Scumulative[m][i] );

	    // convert fmatrix from flavor basis to mass basis
	    // oscillate fmatrix in mass basis
	    // convert back to flavor basis.
	    // don't need to modify pmatrix since it's re-read at each timestep
	    for(flavour f1=e; f1<=mu; f1++)
	      for(flavour f2=e; f2<=mu; f2++){
		assert(fmatrixm[m][i][f1][f2] == fmatrixm[m][i][f1][f2]);
		assert(fmatrixf[m][i][f1][f2] == fmatrixf[m][i][f1][f2]);
	      }
	    fmatrixm[m][i] = SThisStep
	      * Adjoint( U0[m][i] ) 
	      * fmatrixf[m][i] 
	      * U0[m][i]
	      * Adjoint( SThisStep );
	    fmatrixf[m][i] = U0[m][i]
	      * fmatrixm[m][i] 
	      * Adjoint(  U0[m][i] );
	    for(flavour f1=e; f1<=mu; f1++)
	      for(flavour f2=e; f2<=mu; f2++){
		assert(fmatrixm[m][i][f1][f2] == fmatrixm[m][i][f1][f2]);
		assert(fmatrixf[m][i][f1][f2] == fmatrixf[m][i][f1][f2]);
	      }
	    
	    // reset the evolution matrix to identity
	    Y[m][i][msw][0]=Y[m][i][msw][1]=M_PI/2.;
	    Y[m][i][msw][2]=0.;
	    Y[m][i][msw][3]=1.;
	    Y[m][i][msw][4]=Y[m][i][msw][5]=0.;
	    
	    Y[m][i][si][0]=Y[m][i][si][1]=M_PI/2.;
	    Y[m][i][si][2]=0.;
	    Y[m][i][si][3]=1.;
	    Y[m][i][si][4]=Y[m][i][si][5]=0.;

	    // get rate of change of fmatrix from oscillation
	    double hold[4], hnew[4];
	    pauli_decompose(old_fmatrixf[m][i], hold);
	    pauli_decompose(    fmatrixf[m][i], hnew);
	    double oldmag   = sqrt(hold[0]*hold[0] + hold[1]*hold[1] + hold[2]*hold[2]);
	    double newmag   = sqrt(hnew[0]*hnew[0] + hnew[1]*hnew[1] + hnew[2]*hnew[2]);
	    double costheta = (hold[0]*hnew[0] + hold[1]*hnew[1] + hold[2]*hnew[2]) / (newmag*oldmag);
	    assert(costheta-1. < 1e-10);
	    costheta = min(1.,costheta);
	    dtheta_dr_osc[i][m] = (acos(hnew[2]/newmag) - acos(hold[2]/oldmag)) / dr;
	    dphi_dr_osc[i][m] = (atan2(hnew[1],hnew[0]) - atan2(hold[1],hold[0])) / dr;
	  }
	}

	// comment out if not following as a function of r
	if(counterout==step){
	  output=true;
	  counterout=1;
	}
	else counterout++;
	
	if(output==true || finish==true){
	  Outputvsr(fout,foutP,foutf,foutdangledr,foutPvsr,foutFvsr,r,Y,C,A,Scumulative);
	  output=false;
	}

	// adjust step size based on RK error
	// could be moved up to RK section but better left here 
	// in case adjustments are necessary based on new S matrices
	dr = min(dr*pow(accuracy/maxerror,1./max(1,NRKOrder)),increase*dr);
	drmin = 4.*r*numeric_limits<double>::epsilon();
	dr = max(dr,drmin);

      } while(finish==false);

      // if this is not the last domain then carry the S matrix across the discontinuities
      if(d<=ND-2){
	double rminus=rs[d+1]-1.*cgs::units::cm;
	double rplus=rs[d+1]+1.*cgs::units::cm;
	Scumulative = UpdateSm(rminus,rplus,Ye(rminus),Ye(rplus),Y,C,A,Scumulative);
	C0=C;
	C=UpdateC(rplus,Ye(rplus));
	A=UpdateA(C,C0,A0);
      }
      else{ // output at the end of the code
	Outputvsr(fout,foutP,foutf,foutdangledr,foutPvsr,foutFvsr,rs[d+1],Y,C,A,Scumulative);
      }

    }// end of r loop

    for(int i=0;i<=NE-1;i++){
      foutPvsr[i].close();
      foutFvsr[i].close();
    }
    fPvsE.close();
    fFvsE.close();

  }
  catch(OUT_OF_RANGE<unsigned int> OOR){ OOR.Message();}
  catch(OUT_OF_RANGE<int> OOR){ OOR.Message();}
  catch(OUT_OF_RANGE<double> OOR){ OOR.Message();}
  catch(EMPTY E){ E.Message();}
  catch(FUNCTION_ERROR FE){ FE.Message();}
  catch(BASIC_ERROR BE){ BE.Message();}
  catch(...){ UNKNOWN_ERROR("main");}

  cout<<"\nFinished\n\a"; cout.flush();

  return 0;
}


//===//
// B //
//===//
MATRIX<complex<double>,NF,NF> B(vector<double> y){
  MATRIX<complex<double>,NF,NF> s;
  double cPsi1=cos(y[0]),sPsi1=sin(y[0]), cPsi2=cos(y[1]),sPsi2=sin(y[1]), cPsi3=cos(y[2]),sPsi3=sin(y[2]);
  
  s[0][1] = cPsi1 + I*sPsi1*cPsi2;
  sPsi1 *= sPsi2;
  s[0][0] = sPsi1 * (cPsi3 + I*sPsi3);

  s[1][0] = -y[3]*conj(s[0][1]);
  s[1][1] =  y[3]*conj(s[0][0]);

  return s;
}

//===//
// K //
//===//
void K(double r,
       double dr,
       vector<vector<vector<vector<double> > > > &Y,
       vector<vector<vector<MATRIX<complex<double>,NF,NF> > > > &C0,
       vector<vector<vector<vector<double> > > > &A0,
       vector<vector<vector<vector<double> > > > &K){

  MATRIX<complex<double>,NF,NF> VfSI,VfSIbar;  // self-interaction potential
  vector<MATRIX<complex<double>,NF,NF> > VfSIE(NE); // contribution to self-interaction potential from each energy
  MATRIX<complex<double>,NF,NF> VfMSW,VfMSWbar;
  MATRIX<complex<double>,NF,NF> dVfMSWdr,dVfMSWbardr;
  MATRIX<complex<double>,NF,NF> Hf,Hfbar, UU,UUbar;
  vector<double> kk,kkbar, dkk,dkkbar, dkkdr,dkkbardr, QQ,QQbar;
  vector<MATRIX<complex<double>,NF,NF> > CC,dCCdr;
  vector<vector<double> > AA;
  MATRIX<complex<double>,NF,NF> BB,BBbar;
  MATRIX<complex<double>,NF,NF> Sfm,Sfmbar;
  vector<vector<MATRIX<complex<double>,NF,NF> > > 
    Sa(NE,vector<MATRIX<complex<double>,NF,NF> >(NS)), 
    Sabar(NE,vector<MATRIX<complex<double>,NF,NF> >(NS));
  vector<MATRIX<complex<double>,NF,NF> > UWBW(NE);
  vector<MATRIX<complex<double>,NF,NF> > UWBWbar(NE);
  double rrho,drrhodr, YYe,dYYedr;
  MATRIX<double,3,4> JI;
  int i;
  MATRIX<complex<double>,NF,NF> Ha;
  MATRIX<complex<double>,NF,NF> HB;
  vector<double> phase(1);
  vector<double> dvdr(4);
  // *************
  rrho=get_rho(r);
  drrhodr=get_drhodr(rrho,r);
  YYe=get_Ye(r);
  dYYedr=get_dYedr(r);
  VfMSW[e][e]=Ve(rrho,YYe);
  VfMSW[mu][mu]=Vmu(rrho,YYe);
  VfMSWbar=-Conjugate(VfMSW);
  dVfMSWdr[e][e]=dVedr(rrho,drrhodr,YYe,dYYedr);
  dVfMSWdr[mu][mu]=dVmudr(rrho,drrhodr,YYe,dYYedr);
  dVfMSWbardr=-Conjugate(dVfMSWdr);

#pragma omp parallel for schedule(auto) private(Hf,Hfbar,UU,UUbar,kk,kkbar,dkk,dkkbar,dkkdr,dkkbardr,QQ,QQbar,AA,CC,dCCdr,BB,BBbar,Sfm,Sfmbar,JI) firstprivate(Ha,HB,dvdr,phase)
  for(i=0;i<=NE-1;i++){
    Hf  = HfV[matter][i]+VfMSW;
    kk  = k(Hf);
    dkk = deltak(Hf);
    CC  = CofactorMatrices(Hf,kk);
    AA  = MixingMatrixFactors(CC,C0[matter][i],A0[matter][i]);
    UU  = U(dkk,CC,AA);
    BB  = B(Y[matter][i][msw]);
    Sa[i][si] = B(Y[matter][i][si]);
    UWBW[i] = UU * W(Y[matter][i][msw]) * BB * W(Y[matter][i][si]);
    
    Hfbar = HfV[antimatter][i] + VfMSWbar;
    kkbar = kbar(Hfbar);
    dkkbar = deltakbar(Hfbar);
    CC = CofactorMatrices(Hfbar,kkbar);
    AA = MixingMatrixFactors(CC,C0[antimatter][i],A0[antimatter][i]);
    UUbar = Conjugate(U(dkkbar,CC,AA));
    BBbar = B(Y[antimatter][i][msw]);
    Sabar[i][si] = B(Y[antimatter][i][si]);
    UWBWbar[i] = UUbar * W(Y[antimatter][i][msw]) *BBbar * W(Y[antimatter][i][si]);
    
    // ****************
    // Matter section *
    // ****************
    phase[0] = M_2PI*(Y[matter][i][msw][4]-Y[matter][i][msw][5]);
    Ha[0][1]=0.;
    for(int j=0;j<=NF-2;j++)
      for(int k=j+1;k<=NF-1;k++)
	for(flavour f=e;f<=mu;f++)
	  Ha[j][k]+= conj(UU[f][j])*dVfMSWdr[f][f]*UU[f][k];
    
    Ha[0][1] *= I*cgs::constants::hbarc/dkk[0]*exp(I*phase[0]);
    Ha[1][0] = conj(Ha[0][1]);
    
    // HB = -I/cgs::constants::hbarc*Ha*BB;
    HB[0][0]=-I/cgs::constants::hbarc*( Ha[0][1]*BB[1][0] );
    HB[0][1]=-I/cgs::constants::hbarc*( Ha[0][1]*BB[1][1] );
    
    dvdr[0]=real(HB[0][1]);
    dvdr[1]=imag(HB[0][1]);
    dvdr[2]=real(HB[0][0]);
    dvdr[3]=imag(HB[0][0]);
    
    JI = JInverse(Y[matter][i][msw]);
    
    for(int j=0;j<=2;j++){
      K[matter][i][msw][j]=0.;
      for(int k=j;k<=3;k++) K[matter][i][msw][j] += JI[j][k]*dvdr[k];
      K[matter][i][msw][j]*=dr;
    }
    
    K[matter][i][msw][3] = 0.;
    dkkdr = dkdr(UU,dVfMSWdr);
    dCCdr = CofactorMatricesDerivatives(Hf,dVfMSWdr,dkk,dkkdr);
    QQ = Q(UU,dkk,CC,dCCdr);
    
    K[matter][i][msw][4] = (kk[0]+QQ[0])*dr/M_2PI/cgs::constants::hbarc;
    K[matter][i][msw][5] = (kk[1]+QQ[1])*dr/M_2PI/cgs::constants::hbarc;

    // ********************
    // Antimatter section *
    // ********************
    phase[0] = M_2PI*(Y[antimatter][i][msw][4]-Y[antimatter][i][msw][5]);
    Ha[0][1] = 0.;
    for(int j=0;j<=NF-2;j++)
      for(int k=j+1;k<=NF-1;k++)
	for(flavour f=e;f<=mu;f++)
	  Ha[j][k]+=conj(UUbar[f][j])*dVfMSWbardr[f][f]*UUbar[f][k];
    
    Ha[0][1] *= I*cgs::constants::hbarc/dkkbar[0]*exp(I*phase[0]);
    Ha[1][0] = conj(Ha[0][1]);
    
    //HB=-I/cgs::constants::hbarc*Ha*BBbar;
    HB[0][0]=-I/cgs::constants::hbarc*( Ha[0][1]*BBbar[1][0] );
    HB[0][1]=-I/cgs::constants::hbarc*( Ha[0][1]*BBbar[1][1] );
    
    dvdr[0]=real(HB[0][1]);
    dvdr[1]=imag(HB[0][1]);
    dvdr[2]=real(HB[0][0]);
    dvdr[3]=imag(HB[0][0]);

    JI = JInverse(Y[antimatter][i][msw]);

    for(int j=0;j<=2;j++){
      K[antimatter][i][msw][j] = 0.; 
      for(int k=j;k<=3;k++) K[antimatter][i][msw][j] += JI[j][k]*dvdr[k];
      K[antimatter][i][msw][j] *= dr;
    }

    K[antimatter][i][msw][3] = 0.;
    dkkbardr = dkdr(UUbar,dVfMSWbardr);
    dCCdr = CofactorMatricesDerivatives(Hfbar,dVfMSWbardr,dkkbar,dkkbardr);
    QQbar = Q(UUbar,dkkbar,CC,dCCdr);

    K[antimatter][i][msw][4] = (kkbar[0]+QQbar[0])*dr/M_2PI/cgs::constants::hbarc;
    K[antimatter][i][msw][5] = (kkbar[1]+QQbar[1])*dr/M_2PI/cgs::constants::hbarc;

    // *****************************************************************
    // contribution to the self-interaction potential from this energy *
    // *****************************************************************
    Sfm    = UWBW   [i]*Sa   [i][si];
    Sfmbar = UWBWbar[i]*Sabar[i][si];
    VfSIE[i] =     Sfm   *pmatrixm0[    matter][i]*Adjoint(Sfm   )
      - Conjugate( Sfmbar*pmatrixm0[antimatter][i]*Adjoint(Sfmbar) );

  }//end for loop over i

  // ************************************
  // compute self-interaction potential *
  // ************************************
  for(i=0;i<=NE-1;i++){
    VfSI[e ][e ]+=VfSIE[i][e ][e ];
    VfSI[e ][mu]+=VfSIE[i][e ][mu];
    VfSI[mu][e ]+=VfSIE[i][mu][e ];
    VfSI[mu][mu]+=VfSIE[i][mu][mu];
  }

  complex<double> Tr=VfSI[e][e]+VfSI[mu][mu];
  VfSI[e][e]+=Tr;
  VfSI[mu][mu]+=Tr;

  //  VfSI*=NSI*CSI(r);
  VfSIbar=-Conjugate(VfSI);

  // *********************
  // SI part of solution *
  // *********************

#pragma omp parallel for schedule(auto) private(JI) firstprivate(Ha,HB,dvdr)
  for(i=0;i<=NE-1;i++){
    //*********
    // Matter *
    //*********
    Ha = Adjoint(UWBW[i])*VfSI*UWBW[i];

    K[matter][i][si][4]=dr*real(Ha[0][0])/(M_2PI*cgs::constants::hbarc);
    K[matter][i][si][5]=dr*real(Ha[1][1])/(M_2PI*cgs::constants::hbarc);
    
    HB[0][0]=-I/cgs::constants::hbarc*( Ha[0][1]*Sa[i][si][1][0] );
    HB[0][1]=-I/cgs::constants::hbarc*( Ha[0][1]*Sa[i][si][1][1] );
    
    dvdr[0]=real(HB[0][1]);
    dvdr[1]=imag(HB[0][1]);
    dvdr[2]=real(HB[0][0]);
    dvdr[3]=imag(HB[0][0]);
    
    JI=JInverse(Y[matter][i][si]);
    
    for(int j=0;j<=2;j++){
      K[matter][i][si][j]=0.;
      for(int k=j;k<=3;k++) K[matter][i][si][j]+=JI[j][k]*dvdr[k];
      K[matter][i][si][j]*=dr;
    }
    
    K[matter][i][si][3]=0.;
    
    //*************
    // Antimatter *
    //*************
    Ha=Adjoint(UWBWbar[i])*VfSIbar*UWBWbar[i];

    K[antimatter][i][si][4]=dr*real(Ha[0][0])/(M_2PI*cgs::constants::hbarc);
    K[antimatter][i][si][5]=dr*real(Ha[1][1])/(M_2PI*cgs::constants::hbarc);

    HB[0][0]=-I/cgs::constants::hbarc*( Ha[0][1]*Sabar[i][si][1][0] );
    HB[0][1]=-I/cgs::constants::hbarc*( Ha[0][1]*Sabar[i][si][1][1] );

    dvdr[0]=real(HB[0][1]);
    dvdr[1]=imag(HB[0][1]);
    dvdr[2]=real(HB[0][0]);
    dvdr[3]=imag(HB[0][0]);

    JI = JInverse(Y[antimatter][i][si]);

    for(int j=0;j<=2;j++){
      K[antimatter][i][si][j]=0.;
      for(int k=j;k<=3;k++) K[antimatter][i][si][j]+=JI[j][k]*dvdr[k];
      K[antimatter][i][si][j]*=dr;
    }

    K[antimatter][i][si][3]=0.;
  }

}// end of K function


//===========//
// Outputvsr //
//===========//
void Outputvsr(ofstream &fout,
	       ofstream &foutP,
	       ofstream &foutf,
	       ofstream &foutdangledr,
	       ofstream *foutPvsr,
	       ofstream *foutFvsr,
	       double r,
	       vector<vector<vector<vector<double> > > > Y,
	       vector<vector<vector<MATRIX<complex<double>,NF,NF> > > > C0,
	       vector<vector<vector<vector<double> > > > A0,
	       vector<vector<MATRIX<complex<double>,NF,NF> > > Scumulative){

  vector<MATRIX<complex<double>,NF,NF> > VfMSW(NM), dVfMSWdr(NM);
  vector<MATRIX<complex<double>,NF,NF> > VfSI(NM);

  vector<MATRIX<complex<double>,NF,NF> > rhomatrix(NM);

  double rrho=get_rho(r);
  double drrhodr=get_drhodr(rrho,r);

  double YYe=get_Ye(r);
  double dYYedr=get_dYedr(r);

  VfMSW[matter][e][e]=Ve(rrho,YYe);
  VfMSW[matter][mu][mu]=Vmu(rrho,YYe);
  VfMSW[antimatter]=-VfMSW[matter];

  dVfMSWdr[matter][e][e]=dVedr(rrho,drrhodr,YYe,dYYedr);
  dVfMSWdr[matter][mu][mu]=dVmudr(rrho,drrhodr,YYe,dYYedr);
  dVfMSWdr[antimatter]=-dVfMSWdr[matter];

  vector<vector<MATRIX<complex<double>,NF,NF> > >
    Hf(NM,vector<MATRIX<complex<double>,NF,NF> >(NE));
  vector<vector<MATRIX<complex<double>,NF,NF> > >
    UU(NM,vector<MATRIX<complex<double>,NF,NF> >(NE));
  vector<vector<vector<MATRIX<complex<double>,NF,NF> > > > 
    WW(NM,vector<vector<MATRIX<complex<double>,NF,NF> > >(NE,vector<MATRIX<complex<double>,NF,NF> >(NS)));
  vector<vector<vector<MATRIX<complex<double>,NF,NF> > > > 
    BB(NM,vector<vector<MATRIX<complex<double>,NF,NF> > >(NE,vector<MATRIX<complex<double>,NF,NF> >(NS)));
  vector<vector<MATRIX<complex<double>,NF,NF> > > 
    Sm(NM,vector<MATRIX<complex<double>,NF,NF> >(NE)), 
    Smf(NM,vector<MATRIX<complex<double>,NF,NF> >(NE)),
    Sf(NM,vector<MATRIX<complex<double>,NF,NF> >(NE));

  vector<vector<vector<double> > > kk(NM,vector<vector<double> >(NE));
  vector<vector<vector<double> > > dkk(NM,vector<vector<double> >(NE));
  vector<double> ePotentialSum(NE),ebarPotentialSum(NE),heavyPotentialSum(NE);
  double totalANuFlux(0.);
  double totalNuFlux(0.);
  double totalHeavyFlux(0.);
  vector<double> Pe(NE),Pebar(NE),Pheavy(NE);
  vector<double> Pvalues(6);
  vector<double> s(6);
  vector<double> predP((NE+2)*(2));

  MATRIX<complex<double>,NF,NF> p_unosc;
  for(int i=0;i<=NE-1;i++){
    getPunosc(r, matter, i, p_unosc);
    ePotentialSum[i]=real(p_unosc[e][e]);
    heavyPotentialSum[i]=real(p_unosc[mu][mu]);

    getPunosc(r, antimatter, i, p_unosc);
    ebarPotentialSum[i]=real(p_unosc[e][e]);
  }


  for(int i=0;i<=NE-1;i++){
    //---- matter
    Hf[matter][i]  = HfV[matter][i] + VfMSW[matter];
    kk[matter][i]  = k(Hf[matter][i]);
    dkk[matter][i] = deltak(Hf[matter][i]);
    UU[matter][i]  = U(dkk[matter][i],C0[matter][i],A0[matter][i]);
    
    BB[matter][i][msw] = B(Y[matter][i][msw]);
    WW[matter][i][msw] = W(Y[matter][i][msw]);
    BB[matter][i][si] = B(Y[matter][i][si]);
    WW[matter][i][si] = W(Y[matter][i][si]);
    
    Sm[matter][i] = WW[matter][i][msw]
      * BB[matter][i][msw]
      * WW[matter][i][si]
      * BB[matter][i][si]
      * Scumulative[matter][i];
    Smf[matter][i]= Sm[matter][i] * Adjoint(U0[matter][i]);
    Sf[matter][i] = UU[matter][i] * Smf[matter][i];
    
    //---- antimatter
    Hf[antimatter][i]  = HfV[antimatter][i] + VfMSW[antimatter];
    kk[antimatter][i]  = kbar(Hf[antimatter][i]);
    dkk[antimatter][i] = deltakbar(Hf[antimatter][i]);
    UU[antimatter][i]  = Conjugate(U(dkk[antimatter][i],C0[antimatter][i],A0[antimatter][i]));
    
    BB[antimatter][i][msw] = B(Y[antimatter][i][msw]);
    WW[antimatter][i][msw] = W(Y[antimatter][i][msw]);
    BB[antimatter][i][si] = B(Y[antimatter][i][si]);
    WW[antimatter][i][si] = W(Y[antimatter][i][si]);
    
    Sm[antimatter][i] = WW[antimatter][i][msw]
      * BB[antimatter][i][msw]
      * WW[antimatter][i][si]
      * BB[antimatter][i][si]
      * Scumulative[antimatter][i];
    Smf[antimatter][i]= Sm[antimatter][i] * Adjoint(U0[antimatter][i]);
    Sf[antimatter][i] = UU[antimatter][i] * Smf[antimatter][i];
    
    // compute contribution to self interaction potential
    // scattering matrix matter(electron - x) - scattering matrix Antimatter(antielectron - anti-x)
    // what is VfSI[antimatter]?
    // MATRIX<complex<double>,NF,NF> p_unosc_matter, p_unosc_antimatter;
    // getPunosc(r, matter,     i,     p_unosc_matter);
    // getPunosc(r, antimatter, i, p_unosc_antimatter);
    // VfSI[matter] += Sf[    matter][i]*p_unosc_matter    *Adjoint(Sf[    matter][i])
    //   - Conjugate(  Sf[antimatter][i]*p_unosc_antimatter*Adjoint(Sf[antimatter][i]) );
    VfSI[matter] += pmatrixf0[matter][i] - Conjugate(pmatrixf0[antimatter][i]);
  }
  
  complex<double> Tr=VfSI[matter][e][e]+VfSI[matter][mu][mu];
  VfSI[matter][e][e]+=Tr;
  VfSI[matter][mu][mu]+=Tr;

  //VfSI[matter]*=NSI*CSI(r);
  VfSI[antimatter] = -Conjugate(VfSI[matter]);


  fout << r << "\t";
  for(int i=0; i<NE; i++)
    for(state m=matter; m<=antimatter; m++){
      for(flavour f1=e; f1<=mu; f1++)
	for(flavour f2=e; f2<=mu; f2++) {
	  fout << real(Scumulative[m][i][f1][f2] ) << "\t";
	  fout << imag(Scumulative[m][i][f1][f2] ) << "\t";
	}
    }
  fout << endl;
  fout.flush();
  for(int i=0;i<NE;i++){
    fout.flush();
    Pe    [i] = norm(Sf[    matter][i][e ][e ]);
    Pebar [i] = norm(Sf[antimatter][i][e ][e ]);
    Pheavy[i] = norm(Sf[    matter][i][mu][mu]);
  }

  foutP<<r<<"\t"<<Ve(rrho,YYe)<<"\t";//1,2
  foutP<<real(VfSI[    matter][e ][e ])<<"\t"<<real(VfSI[    matter][mu][mu])<<"\t";
  foutP<<real(VfSI[antimatter][e ][e ])<<"\t"<<real(VfSI[antimatter][mu][mu])<<"\t";//3,4,5,6
  Pvalues = averageProbability(Pe,Pebar,Pheavy,ebarPotentialSum,ePotentialSum,heavyPotentialSum);
  totalNuFlux = Pvalues[3];
  totalANuFlux =Pvalues[4];
  totalHeavyFlux = Pvalues[5];
  foutP<<totalNuFlux<<"\t";//Nu,7
  foutP<<totalANuFlux<<"\t";//ANu,8
  foutP<<Pvalues[5]<<"\t";//Heavy,9
  foutP<<Pvalues[0]<<"\t"<<Pvalues[1]<<"\t"<<Pvalues[2]<<"\t";//Pe,Pebar,Pheavy;10,11,12

  predP=predictProbability(Pvalues[3],Pvalues[4],Ve(rrho,YYe),E,ebarPotentialSum,ePotentialSum,heavyPotentialSum);
  foutP<<predP[0]<<"\t"<<predP[1+NE]<<"\t";//13,14
  for(int i=0;i<NE;i++) foutP<<predP[1+i]<<"\t"<<predP[(NE+1)+i+1]<<"\t";//15,16,...2*(NE-1)+15,
  foutP<<predP[(NE+1)*2]<<"\t"<<predP[(NE+1)*2+1]<<"\t";//2*(NE-1)+17,2*(NE-1)+18
  foutP<<endl;
  foutP.flush();

  foutf << r << "\t";
  for(int i=0; i<NE; i++)
    for(state m=matter; m<=antimatter; m++)
      for(flavour f1=e; f1<=mu; f1++)
	for(flavour f2=e; f2<=mu; f2++) {
	  foutf << real( fmatrixf[m][i][f1][f2] ) << "\t";
	  foutf << imag( fmatrixf[m][i][f1][f2] ) << "\t";
	}
  

  foutf << endl;
  foutf.flush();

  foutdangledr << r << "\t";
  for(state m=matter; m<=antimatter; m++)
    for(int i=0; i<NE; i++)
      foutdangledr << dtheta_dr_osc[i][m] << "\t";
  for(state m=matter; m<=antimatter; m++)
    for(int i=0; i<NE; i++)
      foutdangledr << dphi_dr_osc[i][m] << "\t";
  for(state m=matter; m<=antimatter; m++)
    for(int i=0; i<NE; i++)
      foutdangledr << dtheta_dr_interact[i][m] << "\t";
  for(state m=matter; m<=antimatter; m++)
    for(int i=0; i<NE; i++)
      foutdangledr << dphi_dr_interact[i][m] << "\t";
  foutdangledr << endl;
  foutdangledr.flush();
}
