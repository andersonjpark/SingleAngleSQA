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
#include "headers/flavour_basis.h"
#include "headers/eigenvalues.h"
#include "headers/mixing_angles.h"
#include "headers/adiabatic_basis.h"
#include "headers/jacobians.h"
#include "headers/multiEnergy.h"
#include "headers/MNR.h"
#include "headers/misc.h"
#include "headers/State.h"
#include "headers/time_derivatives.h"
#include "headers/output.h"
#include "headers/project/albino.h"
#include "headers/interact.h"
#include "headers/nulib_interface.h"

//======//
// MAIN //
//======//
int main(int argc, char *argv[]){
  string inputfilename;
  ofstream fout,foutP,foutf, foutdangledr;
    
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
    
  init_output(outputfilename, fout, foutP, foutf, foutdangledr);

  //nulib_init(nulibfilename, 0);

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

  // *************************************************
  // set up global variables defined in parameters.h *
  // *************************************************
  // vectors of energies and vacuum eigenvalues
  set_Ebins(E);
  const array<array<double,NF>,NE> kV = set_kV(E);
  const array<MATRIX<complex<double>,NF,NF>,NM> UV = Evaluate_UV();
  const array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> HfV = Evaluate_HfV(kV,UV);
  const array<array<MATRIX<complex<double>,NF,NF>,NF>,NE> CV  = Evaluate_CV(kV, HfV);
  const array<array<array<double,NF>,NF>,NE> AV = Evaluate_AV(kV,HfV,UV);
  
  State s(E);
    
  // **************************************
  // quantities evaluated at inital point *
  // **************************************
    
  // MSW potential matrix
  s.r=rmin;
  s.update_potential(lnrho,temperature,Ye,eP,eBarP,xP,HfV);
    
  // mixing angles to MSW basis at initial point
  for(state m=matter; m<=antimatter; m++){
    for(int i=0;i<=NE-1;i++){
      s.C[m][i]=CofactorMatrices(s.Hf[m][i],s.kk[m][i]);
	
      for(int j=0;j<=NF-1;j++){
	if(real(s.C[m][i][j][mu][e]*CV[i][j][mu][e]) < 0.)
	  s.A[m][i][j][e]=-AV[i][j][e];
	else s.A[m][i][j][e]=AV[i][j][e];
	s.A[m][i][j][mu]=AV[i][j][mu];
      }
      s.U0[m][i]=U(s.dkk[m][i],s.C[m][i],s.A[m][i]);
    }
  }

  // yzhu14 density/potential matrices art rmin
  initialize(s,rmin,eD,eBarD,xD);

  // ***************************************
  // quantities needed for the calculation *
  // ***************************************
  double dr,drmin,dr_this_step;
    
  double maxerror,increase=3.;
  bool repeat, finish, resetflag, output;
  int counterout;
    
  // comment out if not following as a function of r
    
  // ***************************************
  // variables followed as a function of r *
  // ***************************************
    
    
  // ************************
  // Runge-Kutta quantities *
  // ************************
    
  array<array<array<array<array<double,NY>,NS>,NE>,NM>,NRK> Ks;
    
    
  // *****************************************
  // initialize at beginning of every domain *
  // *****************************************
  dr=1e-3*cgs::units::cm;
  drmin=4.*s.r*numeric_limits<double>::epsilon();
      
  // *************************************************
  // comment out if not following as a function of r *
  // *************************************************
	
  finish=output=false;
  counterout=1;
  s.update_potential(lnrho,temperature,Ye,eP,eBarP,xP,HfV);
  Outputvsr(fout,foutP,foutf,foutdangledr,s,eP,eBarP,xP);
	
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
	  
    State sReset = s;
	  
    for(state m=matter; m<=antimatter; m++)
      for(int i=0; i<NE; i++)
	for(flavour f1=e; f1<=mu; f1++)
	  for(flavour f2=e; f2<=mu; f2++)
	    assert(s.fmatrixf[m][i][f1][f2] == s.fmatrixf[m][i][f1][f2]);
	
    // beginning of RK section
    do{ 
      repeat=false;
      for(int k=0;k<=NRK-1;k++){
	s = sReset;
	s.r +=AA[k]*dr;

	for(state m = matter; m <= antimatter; m++)
	  for(int i=0;i<=NE-1;i++)
	    for(solution x=msw;x<=si;x++)
	      for(int j=0;j<=NY-1;j++)
		for(int l=0;l<=k-1;l++)
		  s.Y[m][i][x][j] += BB[k][l] * Ks[l][m][i][x][j];

	s.update_potential(lnrho,temperature,Ye,eP,eBarP,xP,HfV);
	Ks[k] = K(dr,s);
      }
	  
      // increment all quantities and update C and A arrays
      s.r= sReset.r+dr;
      s.Y = sReset.Y;
      maxerror=0.;
      for(state m=matter;m<=antimatter;m++){
	for(int i=0;i<=NE-1;i++){
	  for(solution x=msw;x<=si;x++){
	    for(int j=0;j<=NY-1;j++){
	      double Yerror = 0.;
	      for(int k=0;k<=NRK-1;k++){
		assert(CC[k] == CC[k]);
		assert(Ks[k][m][i][x][j] == Ks[k][m][i][x][j]);
		s.Y[m][i][x][j] += CC[k] * Ks[k][m][i][x][j];
		Yerror += (CC[k]-DD[k]) * Ks[k][m][i][x][j];
		assert(s.Y[m][i][x][j] == s.Y[m][i][x][j]);
	      }
	      maxerror = max(maxerror, fabs(Yerror));
	    }
	  }
	}
      }
      s.update_potential(lnrho,temperature,Ye,eP,eBarP,xP,HfV);
	    
      // decide whether to accept step, if not adjust step size
      dr_this_step = dr;
      if(maxerror>accuracy){
	dr *= 0.9 * pow(accuracy/maxerror, 1./(NRKOrder-1.));
	if(dr > drmin) repeat=true;
      }

      // reset integration variables to those at beginning of step
      if(repeat==true) s = sReset;
	  
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
	s.Scumulative[m][i] = s.SThisStep[m][i] * s.Scumulative[m][i];

	// convert fmatrix from flavor basis to (reset-point) mass basis
	// evolve fmatrix from reset-point to current-point mass basis
	// convert fmatrix from (current-point) mass basis to flavor basis
	MATRIX<complex<double>,NF,NF> SfThisStep =
	  s.UU[m][i]
	  * s.SThisStep[m][i]
	  * Adjoint(sReset.UU[m][i]);
	s.fmatrixf[m][i] = SfThisStep * sReset.fmatrixf[m][i] * Adjoint(SfThisStep);
	    
	// reset the evolution matrix to identity
	s.Y[m][i] = YIdentity;

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
      s.update_potential(lnrho,temperature,Ye,eP,eBarP,xP,HfV);
      Outputvsr(fout,foutP,foutf,foutdangledr,s,eP,eBarP,xP);
      output=false;
    }

    // adjust step size based on RK error
    // could be moved up to RK section but better left here 
    // in case adjustments are necessary based on new S matrices
    dr = min(dr*pow(accuracy/maxerror,1./max(1,NRKOrder)),increase*dr);
    drmin = 4.*s.r*numeric_limits<double>::epsilon();
    dr = max(dr,drmin);

  } while(finish==false);

  s.update_potential(lnrho,temperature,Ye,eP,eBarP,xP,HfV);
  Outputvsr(fout,foutP,foutf,foutdangledr,s,eP,eBarP,xP);

  cout<<"\nFinished\n\a"; cout.flush();

  return 0;
}


