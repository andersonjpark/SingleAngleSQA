/*
//  Copyright (c) 2018, James Kneller and Sherwood Richers
//
//  This file is part of IsotropicSQA.
//
//  IsotropicSQA is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  IsotropicSQA is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with IsotropicSQA.  If not, see <http://www.gnu.org/licenses/>.
//
*/

#include <complex>
using std::complex;
using std::abs;
#include<iostream>
using std::cout;
#include<ostream>
using std::endl;
using std::flush;
#include<fstream>
using std::ifstream;
using std::ofstream;
#include<algorithm>
using std::min;
using std::max;
#include<string>
using std::string;
#include<limits>
using std::numeric_limits;
#include<vector>
using std::vector;
#include<array>
using std::array;
#include<hdf5.h>

#include "headers/DISCONTINUOUS.h"


// headers
#include "headers/MATRIX.h"
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
#include "headers/misc.h"
#include "headers/State.h"


MATRIX<complex<double>,NF,NF> B(array<double,NY> y);
array<array<array<array<double,NY>,NS>,NE>,NM> K(double dr,
       array<array<array<array<double,NY>,NS>,NE>,NM>& Y,
       array<array<array<MATRIX<complex<double>,NF,NF>,NF>,NE>,NM>& C0,
       array<array<array<array<double,NF>,NF>,NE>,NM> &A0,
       State& s);
void Outputvsr(ofstream &fout,
	       ofstream &foutP,
	       ofstream &foutf,
	       ofstream &foutdangledr,
	       const array<array<array<array<double,NY>,NS>,NE>,NM>& Y,
	       const array<array<array<MATRIX<complex<double>,NF,NF>,NF>,NE>,NM>& C0,
	       const array<array<array<array<double,NF>,NF>,NE>,NM>& A0,
	       const State& s,
	       const array<DISCONTINUOUS,NE>& eP,
	       const array<DISCONTINUOUS,NE>& eBarP,
	       const array<DISCONTINUOUS,NE>& xP);

#include "headers/update.h"
#include "headers/project/albino.h"
//#include "headers/project/test_case_B.h"

void interact(double dr, State& s,
	      const array<DISCONTINUOUS,NE>& eD,
	      const array<DISCONTINUOUS,NE>& eBarD,
	      const array<DISCONTINUOUS,NE>& xD){
  
  // save old fmatrix
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> old_fmatrixf = s.fmatrixf;

  // let neutrinos interact
  my_interact(s.fmatrixf, dr, s, eD,eBarD,xD);
  for(int i=0; i<NE; i++){
    for(state m=matter; m<=antimatter; m++){
      s.fmatrixm[m][i] = Adjoint(s.U0[m][i]) * s.fmatrixf[m][i] * s.U0[m][i];
    }
  }
  
  // loop through getting rotation matrices
  double hold[4];
  double hnew[4];
  for(int i=0; i<NE; i++){
    for(state m=matter; m<=antimatter; m++){
      pauli_decompose(old_fmatrixf[m][i], hold);
      pauli_decompose(  s.fmatrixf[m][i], hnew);

      // get the theta and phi contribution
      double oldmag2   = hold[0]*hold[0] + hold[1]*hold[1] + hold[2]*hold[2];
      double newmag2   = hnew[0]*hnew[0] + hnew[1]*hnew[1] + hnew[2]*hnew[2];
      if(oldmag2==0 or newmag2==0){
	continue;
      }
      s.dtheta_dr_interact[i][m] = (acos(hnew[2]/sqrt(newmag2)) - acos(hold[2]/sqrt(oldmag2))) / dr;
      s.dphi_dr_interact[i][m] = (atan2(hnew[1],hnew[0]) - atan2(hold[1],hold[0])) / dr;

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
      s.Scumulative[m][i] = MATRIX<complex<double>,NF,NF>(Adjoint(s.U0[m][i]) * R * s.U0[m][i] * s.Scumulative[m][i]);
    }
  }
}

#include "headers/nulib_interface.h"

//======//
// MAIN //
//======//
int main(int argc, char *argv[]){
    string inputfilename;
    ofstream fout,foutC,foutP,foutS, foutf, foutdangledr;
    
    inputfilename=string(argv[1]);
    ifstream fin(inputfilename.c_str());
    
    // load the nulib table
    const string nulibfilename = get_parameter<string>(fin,"nulibfilename");
    const string potential_directory = get_parameter<string>(fin,"potential_directory");
    const string rhofilename = get_parameter<string>(fin, "rhofilename");
    const string Yefilename = get_parameter<string>(fin, "Yefilename");
    const string temperaturefilename = get_parameter<string>(fin, "temparaturefilename");
    const string outputfilename = get_parameter<string>(fin, "outputfilename");
    const double rmin = get_parameter<double>(fin, "rmin"); // cm
    const double rmax = get_parameter<double>(fin, "rmax"); // cm
    const double accuracy = get_parameter<double>(fin, "accuracy");
    const bool do_interact = get_parameter<bool>(fin, "do_interact");
    const int out_every = get_parameter<int>(fin, "out_every");
    fin.close();
    
    const string outputfilenamestem = outputfilename+"/";

    nulib_init(nulibfilename, 0);

    // interpolation variables
    DISCONTINUOUS lnrho, Ye, temperature; // rho is the mass density
    array<DISCONTINUOUS,NE> eP,eBarP,xP;
    array<DISCONTINUOUS,NE> eD,eBarD,xD;

    
    // load rho and Ye data
    lnrho.Open(rhofilename,'#');
    Ye.Open(Yefilename,'#');
    temperature.Open(temperaturefilename,'#');
    assert(rmin >= lnrho.XMin());
    assert(rmin <= lnrho.XMax());

    lnrho = lnrho.copy_logy();
    
    // load and compute spectral data
    for(int i=0;i<=NE-1;i++){
      eP   [i].Open(potential_directory+"/potential_s1_g"+patch::to_string(i+1)+"_.txt",'#');
      eBarP[i].Open(potential_directory+"/potential_s2_g"+patch::to_string(i+1)+"_.txt",'#');
      xP   [i].Open(potential_directory+"/potential_s3_g"+patch::to_string(i+1)+"_.txt",'#');
      eD   [i].Open(potential_directory+"/density_s1_g"+patch::to_string(i+1)+"_.txt",'#');
      eBarD[i].Open(potential_directory+"/density_s2_g"+patch::to_string(i+1)+"_.txt",'#');
      xD   [i].Open(potential_directory+"/density_s3_g"+patch::to_string(i+1)+"_.txt",'#');
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
    
    ofstream fPvsE, fFvsE;
    
    // *************************************************
    // set up global variables defined in parameters.h *
    // *************************************************

    set_Ebins(E);
    State s(E);
    // vectors of energies and vacuum eigenvalues
    //s.kV = set_kV(E);
    
    // vaccum mixing matrices and Hamiltonians
    //Evaluate_UV();
    
    //HfV[matter] = vector<MATRIX<complex<double>,NF,NF> >(NE);
    //HfV[antimatter] = vector<MATRIX<complex<double>,NF,NF> >(NE);
    //Evaluate_HfV(s.kV);
    
    // cofactor matrices in vacuum
    //CV = Evaluate_CV(s.kV);
    
    // mixing matrix element prefactors in vacuum
    //AV=vector<vector<vector<double> > >(NE,vector<vector<double> >(NF,vector<double>(NF)));
    //Evaluate_AV(s.kV);
    
    // **************************************
    // quantities evaluated at inital point *
    // **************************************
    
    // MSW potential matrix
    s.r=rmin;
    s.update_background(lnrho,temperature,Ye,eD,eBarD,xD,eP,eBarP,xP);
    
    
    // cofactor matrices at initial point - will be recycled as cofactor matrices at beginning of every step
    array<array<array<MATRIX<complex<double>,NF,NF>,NF>,NE>,NM> C0, C;

    // mixing matrix element prefactors at initial point - will be recycled like C0
    array<array<array<array<double,NF>,NF>,NE>,NM> A0, A;
    
    // mixing angles to MSW basis at initial point
    for(state m=matter; m<=antimatter; m++){
      for(int i=0;i<=NE-1;i++){
	MATRIX<complex<double>,NF,NF> Hf0=s.HfV[m][i]+s.VfMSW[m];
	array<double,NF> k0=k(Hf0);
	array<double,NF-1> deltak0=deltak(Hf0);
	C0[m][i]=CofactorMatrices(Hf0,k0);
	
	for(int j=0;j<=NF-1;j++){
	  if(real(C0[m][i][j][mu][e]*s.CV[i][j][mu][e]) < 0.)
	    A0[m][i][j][e]=-s.AV[i][j][e];
	  else A0[m][i][j][e]=s.AV[i][j][e];
	  A0[m][i][j][mu]=s.AV[i][j][mu];
	}
	s.U0[m][i]=U(deltak0,C0[m][i],A0[m][i]);
      }
    }

    // yzhu14 density/potential matrices art rmin
    initialize(s,rmin,eD,eBarD,xD);

    // ***************************************
    // quantities needed for the calculation *
    // ***************************************
    double r0,dr,drmin,dr_this_step;
    
    double maxerror,increase=3.;
    bool repeat, finish, resetflag, output;
    int counterout;
    
    // comment out if not following as a function of r
    
    // ***************************************
    // variables followed as a function of r *
    // ***************************************
    
    array<array<array<array<double,NY>,NS>,NE>,NM> Y, dY, Y0, Yerror;
    
    // cofactor matrices
    C=C0;
    // mixing matrix prefactors
    A=A0;
        
    // ************************
    // Runge-Kutta quantities *
    // ************************
    
    array<array<array<array<array<double,NY>,NS>,NE>,NM>,NRK> Ks;
    
    // temporaries
    MATRIX<complex<double>,NF,NF> SSMSW,SSSI,SThisStep;
    
    
      // *****************************************
      // initialize at beginning of every domain *
      // *****************************************
      dr=1e-3*cgs::units::cm;
      drmin=4.*s.r*numeric_limits<double>::epsilon();
      
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
      s.update_background(lnrho,temperature,Ye,eD,eBarD,xD,eP,eBarP,xP);
      Outputvsr(fout,foutP,foutf,foutdangledr,Y,C,A,s,eP,eBarP,xP);
	
      for(state m=matter; m<=antimatter; m++)
	for(int i=0; i<NE; i++)
	  for(flavour f1=e; f1<=mu; f1++)
	    for(flavour f2=e; f2<=mu; f2++)
	      assert(s.fmatrixf[m][i][f1][f2] == s.fmatrixf[m][i][f1][f2]);
      
      // ***********************
      // start the loop over r *
      // ***********************
      do{ 
	double intkm = int(s.r/1e5)*1e5;
	if(s.r - intkm <= dr){
	  cout << s.r/1e5 << " " << dr << " " << s.rho << " " << s.T << " " << s.Ye << endl;
	  cout.flush();
	}

	if(s.r+dr>rmax){
	  dr=rmax-s.r;
	  finish=true;
	  output=true;
	}
	  
	r0=s.r;
	Y0=Y;
	C0=C;
	A0=A;
	  
	for(state m=matter; m<=antimatter; m++)
	  for(int i=0; i<NE; i++)
	    for(flavour f1=e; f1<=mu; f1++)
	      for(flavour f2=e; f2<=mu; f2++)
		assert(s.fmatrixf[m][i][f1][f2] == s.fmatrixf[m][i][f1][f2]);
	
	// beginning of RK section
	do{ 
	  repeat=false;
	  for(int k=0;k<=NRK-1;k++){
	    s.r=r0+AA[k]*dr;
	    Y=Y0;

	    for(state m = matter; m <= antimatter; m++)
	      for(int i=0;i<=NE-1;i++)
		for(solution x=msw;x<=si;x++)
		  for(int j=0;j<=NY-1;j++)
		    for(int l=0;l<=k-1;l++)
		      Y[m][i][x][j] += BB[k][l] * Ks[l][m][i][x][j];

	    s.update_background(lnrho,temperature,Ye,eD,eBarD,xD,eP,eBarP,xP);
	    Ks[k] = K(dr,Y,C,A,s);
	  }
	  
	  // increment all quantities and update C and A arrays
	  s.r=r0+dr;
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
	  
	  C=UpdateC(s,lnrho,Ye);
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
	    s.r=r0;
	    Y=Y0;
	    C=C0;
	    A=A0;
	  }
	  
	}while(repeat==true); // end of RK section

	// interact with the matter
	if(do_interact)
	  interact(dr_this_step, s,eD,eBarD,xD);
	for(state m=matter; m<=antimatter; m++)
	  for(int i=0; i<NE; i++)
	    for(flavour f1=e; f1<=mu; f1++)
	      for(flavour f2=e; f2<=mu; f2++)
		assert(s.fmatrixf[m][i][f1][f2] == s.fmatrixf[m][i][f1][f2]);

	// accumulate S and reset variables
	array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> old_fmatrixf = s.fmatrixf;
	for(state m=matter;m<=antimatter;m++){
	  for(int i=0;i<=NE-1;i++){
	    SSMSW = W(Y[m][i][msw])*B(Y[m][i][msw]);
	    SSSI  = W(Y[m][i][si ])*B(Y[m][i][si ]);
	    SThisStep = SSMSW*SSSI;
	    s.Scumulative[m][i]=MATRIX<complex<double>,NF,NF>(SThisStep*s.Scumulative[m][i] );

	    // convert fmatrix from flavor basis to mass basis
	    // oscillate fmatrix in mass basis
	    // convert back to flavor basis.
	    // don't need to modify pmatrix since it's re-read at each timestep
	    for(flavour f1=e; f1<=mu; f1++)
	      for(flavour f2=e; f2<=mu; f2++){
		assert(s.fmatrixm[m][i][f1][f2] == s.fmatrixm[m][i][f1][f2]);
		assert(s.fmatrixf[m][i][f1][f2] == s.fmatrixf[m][i][f1][f2]);
	      }
	    s.fmatrixm[m][i] = SThisStep
	      * Adjoint( s.U0[m][i] ) 
	      * s.fmatrixf[m][i] 
	      * s.U0[m][i]
	      * Adjoint( SThisStep );
	    s.fmatrixf[m][i] = s.U0[m][i]
	      * s.fmatrixm[m][i] 
	      * Adjoint(  s.U0[m][i] );
	    for(flavour f1=e; f1<=mu; f1++)
	      for(flavour f2=e; f2<=mu; f2++){
		assert(s.fmatrixm[m][i][f1][f2] == s.fmatrixm[m][i][f1][f2]);
		assert(s.fmatrixf[m][i][f1][f2] == s.fmatrixf[m][i][f1][f2]);
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
	    pauli_decompose(  s.fmatrixf[m][i], hnew);
	    double oldmag   = sqrt(hold[0]*hold[0] + hold[1]*hold[1] + hold[2]*hold[2]);
	    double newmag   = sqrt(hnew[0]*hnew[0] + hnew[1]*hnew[1] + hnew[2]*hnew[2]);
	    double costheta = (hold[0]*hnew[0] + hold[1]*hnew[1] + hold[2]*hnew[2]) / (newmag*oldmag);
	    assert(costheta-1. < 1e-10);
	    costheta = min(1.,costheta);
	    s.dtheta_dr_osc[i][m] = (acos(hnew[2]/newmag) - acos(hold[2]/oldmag)) / dr;
	    s.dphi_dr_osc[i][m] = (atan2(hnew[1],hnew[0]) - atan2(hold[1],hold[0])) / dr;
	  }
	}

	// comment out if not following as a function of r
	if(counterout==out_every){
	  output=true;
	  counterout=1;
	}
	else counterout++;
	
	if(output==true || finish==true){
	  s.update_background(lnrho,temperature,Ye,eD,eBarD,xD,eP,eBarP,xP);
	  Outputvsr(fout,foutP,foutf,foutdangledr,Y,C,A,s,eP,eBarP,xP);
	  output=false;
	}

	// adjust step size based on RK error
	// could be moved up to RK section but better left here 
	// in case adjustments are necessary based on new S matrices
	dr = min(dr*pow(accuracy/maxerror,1./max(1,NRKOrder)),increase*dr);
	drmin = 4.*s.r*numeric_limits<double>::epsilon();
	dr = max(dr,drmin);

      } while(finish==false);

      s.update_background(lnrho,temperature,Ye,eD,eBarD,xD,eP,eBarP,xP);
      Outputvsr(fout,foutP,foutf,foutdangledr,Y,C,A,s,eP,eBarP,xP);
    fPvsE.close();
    fFvsE.close();

  cout<<"\nFinished\n\a"; cout.flush();

  return 0;
}


//===//
// B //
//===//
MATRIX<complex<double>,NF,NF> B(array<double,NY> y){
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
array<array<array<array<double,NY>,NS>,NE>,NM> K(double dr,
       array<array<array<array<double,NY>,NS>,NE>,NM>& Y,
       array<array<array<MATRIX<complex<double>,NF,NF>,NF>,NE>,NM>& C0,
       array<array<array<array<double,NF>,NF>,NE>,NM> &A0,
       State& s){

  array<array<array<array<double,NY>,NS>,NE>,NM> K;
  array<MATRIX<complex<double>,NF,NF>,NM> VfSI;  // self-interaction potential
  array<MATRIX<complex<double>,NF,NF>,NE> VfSIE; // contribution to self-interaction potential from each energy
  array<array<array<MATRIX<complex<double>,NF,NF>,NS>,NE>,NM> Sa;
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> UWBW;
  
#pragma omp parallel for
  for(int i=0;i<=NE-1;i++){
    array<array<double,NF>,NM> kk;
    array<array<double,4>,NM> dvdr;

    array<array<MATRIX<complex<double>,NF,NF>,NF>,NM> CC;
    array<MATRIX<complex<double>,NF,NF>,NM> Ha,HB;
    array<array<MATRIX<complex<double>,NF,NF>,NF>,NM> dCCdr;
    
    for(int m=matter; m<=antimatter; m++){
      MATRIX<complex<double>,NF,NF> Hf = s.HfV[m][i]+s.VfMSW[m];
      kk[m] = k(Hf);
      array<double,NF-1> dkk = deltak(Hf);
      CC[m]  = CofactorMatrices(Hf,kk[m]);
      array<array<double,NF>,NF> AA = MixingMatrixFactors(CC[m],C0[m][i],A0[m][i]);
      MATRIX<complex<double>,NF,NF> UU = U(dkk,CC[m],AA);
      MATRIX<complex<double>,NF,NF> BB = B(Y[m][i][msw]);
      Sa[m][i][si] = B(Y[m][i][si]);
      UWBW[m][i] = UU * W(Y[m][i][msw]) * BB * W(Y[m][i][si]);

      array<double,NF-1> phase;
      phase[0] = M_2PI*(Y[m][i][msw][4]-Y[m][i][msw][5]);
      Ha[m][0][1]=0.;
      for(int j=0;j<=NF-2;j++)
	for(int k=j+1;k<=NF-1;k++)
	  for(flavour f=e;f<=mu;f++)
	    Ha[m][j][k]+= conj(UU[f][j])*s.dVfMSWdr[m][f][f]*UU[f][k];
    
      Ha[m][0][1] *= I*cgs::constants::hbarc/dkk[0]*exp(I*phase[0]);
      Ha[m][1][0] = conj(Ha[m][0][1]);
    
      // HB = -I/cgs::constants::hbarc*Ha*BB;
      HB[m][0][0]=-I/cgs::constants::hbarc*( Ha[m][0][1]*BB[1][0] );
      HB[m][0][1]=-I/cgs::constants::hbarc*( Ha[m][0][1]*BB[1][1] );

      dvdr[m][0]=real(HB[m][0][1]);
      dvdr[m][1]=imag(HB[m][0][1]);
      dvdr[m][2]=real(HB[m][0][0]);
      dvdr[m][3]=imag(HB[m][0][0]);

      MATRIX<double,3,4> JI = JInverse(Y[m][i][msw]);
      
      array<double,NF> dkkdr = dkdr(UU,s.dVfMSWdr[m]);
      dCCdr[m] = CofactorMatricesDerivatives(Hf,s.dVfMSWdr[m],dkkdr);
      array<double,NF> QQ =  Q(UU,dkk,CC[m],dCCdr[m]);

      for(int j=0;j<=2;j++){
	K[m][i][msw][j]=0.;
	for(int k=j;k<=3;k++)
	  K[m][i][msw][j] += JI[j][k]*dvdr[m][k];
      }
      K[m][i][msw][3] = 0.;
      K[m][i][msw][4] = (kk[m][0]+QQ[0])/M_2PI/cgs::constants::hbarc;
      K[m][i][msw][5] = (kk[m][1]+QQ[1])/M_2PI/cgs::constants::hbarc;
      for(int j=0;j<NY;j++)
	K[m][i][msw][j]*=dr;
      
      // contribution to the self-interaction potential from this energy
      MATRIX<complex<double>,NF,NF> Sfm    = UWBW[m][i]*Sa[m][i][si];
      MATRIX<complex<double>,NF,NF> VfSIE = Sfm * s.pmatrixm0[m][i] * Adjoint(Sfm);
      if(m==antimatter) VfSIE = -Conjugate(VfSIE);
      #pragma omp critical
      VfSI[matter] += VfSIE;
    }
  }//end for loop over i

  #pragma omp single
  VfSI[antimatter]=-Conjugate(VfSI[matter]);

  // *********************
  // SI part of solution *
  // *********************

#pragma omp parallel for
  for(int i=0;i<=NE-1;i++){
    //*********
    // Matter *
    //*********
    MATRIX<complex<double>,NF,NF> Ha,HB;
    Ha = Adjoint(UWBW[matter][i])*VfSI[matter]*UWBW[matter][i];

    K[matter][i][si][4]=dr*real(Ha[0][0])/(M_2PI*cgs::constants::hbarc);
    K[matter][i][si][5]=dr*real(Ha[1][1])/(M_2PI*cgs::constants::hbarc);
    
    HB[0][0]=-I/cgs::constants::hbarc*( Ha[0][1]*Sa[matter][i][si][1][0] );
    HB[0][1]=-I/cgs::constants::hbarc*( Ha[0][1]*Sa[matter][i][si][1][1] );
    
    array<double,4> dvdr;
    dvdr[0]=real(HB[0][1]);
    dvdr[1]=imag(HB[0][1]);
    dvdr[2]=real(HB[0][0]);
    dvdr[3]=imag(HB[0][0]);
    
    MATRIX<double,3,4> JI = JInverse(Y[matter][i][si]);
    
    for(int j=0;j<=2;j++){
      K[matter][i][si][j]=0.;
      for(int k=j;k<=3;k++) K[matter][i][si][j]+=JI[j][k]*dvdr[k];
      K[matter][i][si][j]*=dr;
    }
    
    K[matter][i][si][3]=0.;
    
    //*************
    // Antimatter *
    //*************
    Ha=Adjoint(UWBW[antimatter][i])*VfSI[antimatter]*UWBW[antimatter][i];

    K[antimatter][i][si][4]=dr*real(Ha[0][0])/(M_2PI*cgs::constants::hbarc);
    K[antimatter][i][si][5]=dr*real(Ha[1][1])/(M_2PI*cgs::constants::hbarc);

    HB[0][0]=-I/cgs::constants::hbarc*( Ha[0][1]*Sa[antimatter][i][si][1][0] );
    HB[0][1]=-I/cgs::constants::hbarc*( Ha[0][1]*Sa[antimatter][i][si][1][1] );

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

  return K;
}// end of K function


//===========//
// Outputvsr //
//===========//
void Outputvsr(ofstream &fout,
	       ofstream &foutP,
	       ofstream &foutf,
	       ofstream &foutdangledr,
	       const array<array<array<array<double,NY>,NS>,NE>,NM>& Y,
	       const array<array<array<MATRIX<complex<double>,NF,NF>,NF>,NE>,NM>& C0,
	       const array<array<array<array<double,NF>,NF>,NE>,NM>& A0,
	       const State& s,
	       const array<DISCONTINUOUS,NE>& eP,
	       const array<DISCONTINUOUS,NE>& eBarP,
	       const array<DISCONTINUOUS,NE>& xP){

  array<MATRIX<complex<double>,NF,NF>,NM> VfSI;

  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> Hf, UU;
  array<array<array<MATRIX<complex<double>,NF,NF>,NM>,NE>,NS> WW,BB;
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> Sm, Smf, Sf;

  array<array<array<double,NF>,NE>,NM> kk;
  array<array<array<double,NF-1>,NE>,NM> dkk;
  array<double,NE> ePotentialSum,ebarPotentialSum,heavyPotentialSum;
  double totalANuFlux(0.);
  double totalNuFlux(0.);
  double totalHeavyFlux(0.);
  array<double,NE> Pe,Pebar,Pheavy;
  array<double,6> Pvalues;
  array<double,(NE+2)*2> predP;

  MATRIX<complex<double>,NF,NF> p_unosc;
  for(int i=0;i<=NE-1;i++){
    ePotentialSum[i]=eP[i](s.r);
    ebarPotentialSum[i]=eBarP[i](s.r);
    heavyPotentialSum[i]=xP[i](s.r);
  }


  for(int i=0;i<=NE-1;i++){
    //---- matter
    Hf[matter][i]  = s.HfV[matter][i] + s.VfMSW[matter];
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
      * s.Scumulative[matter][i];
    Smf[matter][i]= Sm[matter][i] * Adjoint(s.U0[matter][i]);
    Sf[matter][i] = UU[matter][i] * Smf[matter][i];
    
    //---- antimatter
    Hf[antimatter][i]  = s.HfV[antimatter][i] + s.VfMSW[antimatter];
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
      * s.Scumulative[antimatter][i];
    Smf[antimatter][i]= Sm[antimatter][i] * Adjoint(s.U0[antimatter][i]);
    Sf[antimatter][i] = UU[antimatter][i] * Smf[antimatter][i];
    
    // compute contribution to self interaction potential
    // scattering matrix matter(electron - x) - scattering matrix Antimatter(antielectron - anti-x)
    // what is VfSI[antimatter]?
    // MATRIX<complex<double>,NF,NF> p_unosc_matter, p_unosc_antimatter;
    // getPunosc(r, matter,     i,     p_unosc_matter);
    // getPunosc(r, antimatter, i, p_unosc_antimatter);
    // VfSI[matter] += Sf[    matter][i]*p_unosc_matter    *Adjoint(Sf[    matter][i])
    //   - Conjugate(  Sf[antimatter][i]*p_unosc_antimatter*Adjoint(Sf[antimatter][i]) );
    VfSI[matter] += s.pmatrixf0[matter][i] - Conjugate(s.pmatrixf0[antimatter][i]);
  }
  
  complex<double> Tr=VfSI[matter][e][e]+VfSI[matter][mu][mu];
  VfSI[matter][e][e]+=Tr;
  VfSI[matter][mu][mu]+=Tr;

  VfSI[antimatter] = -Conjugate(VfSI[matter]);


  fout << s.r << "\t";
  for(int i=0; i<NE; i++)
    for(state m=matter; m<=antimatter; m++){
      for(flavour f1=e; f1<=mu; f1++)
	for(flavour f2=e; f2<=mu; f2++) {
	  fout << real(s.Scumulative[m][i][f1][f2] ) << "\t";
	  fout << imag(s.Scumulative[m][i][f1][f2] ) << "\t";
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

  foutP<<s.r<<"\t"<<Ve(s.rho,s.Ye)<<"\t";//1,2
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

  predP=predictProbability(Pvalues[3],Pvalues[4],Ve(s.rho,s.Ye),E,ebarPotentialSum,ePotentialSum,heavyPotentialSum);
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
