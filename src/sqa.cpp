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



// headers
#include "DISCONTINUOUS.h"
#include "MATRIX.h"
#include "parameters.h"
#include "potentials.h"
#include "flavour_basis.h"
#include "eigenvalues.h"
#include "mixing_angles.h"
#include "adiabatic_basis.h"
#include "jacobians.h"
#include "misc.h"
#include "State.h"
#include "time_derivatives.h"
#include "IO.h"
#include "albino.h"
#include "nulib_interface.h"
#include "evolve.h"

//======//
// MAIN //
//======//
int main(int argc, char *argv[]){
  string inputfilename;
    
  inputfilename=string(argv[1]);
  ifstream fin(inputfilename.c_str());
    
  // load the nulib table
  const string nulibfilename = get_parameter<string>(fin,"nulibfilename");
  const string input_directory = get_parameter<string>(fin,"input_directory");
  const double rmin = get_parameter<double>(fin, "rmin"); // cm
  const double rmax = get_parameter<double>(fin, "rmax"); // cm
  const double dr0 = get_parameter<double>(fin, "dr0"); // cm
  const double dr_block_max = get_parameter<double>(fin, "dr_block_max"); // cm
  const double accuracy = get_parameter<double>(fin, "accuracy");
  const bool do_oscillate = get_parameter<bool>(fin, "do_oscillate");
  const bool do_interact = get_parameter<bool>(fin, "do_interact");
  const double target_impact = get_parameter<double>(fin, "target_impact");
  const double increase = get_parameter<double>(fin, "increase");
  fin.close();
    
  //nulib_init(nulibfilename, 0);

  DISCONTINUOUS lnrho, Ye, temperature;
  array<array<array<DISCONTINUOUS,NF>,NE>,NM> P_unosc, D_unosc;
  load_input_data(input_directory, lnrho, Ye, temperature, P_unosc, D_unosc);
  assert(rmin >= lnrho.XMin());
  assert(rmin <= lnrho.XMax());

  // *************************************************
  // set up global variables defined in parameters.h *
  // *************************************************
  // vectors of energies and vacuum eigenvalues
  const array<double,NE> E = set_Ebins();
  const array<array<double,NF>,NE> kV = set_kV(E);
  const array<MATRIX<complex<double>,NF,NF>,NM> UV = Evaluate_UV();
  const array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> HfV = Evaluate_HfV(kV,UV);
  const array<array<MATRIX<complex<double>,NF,NF>,NF>,NE> CV  = Evaluate_CV(kV, HfV);
  const array<array<array<double,NF>,NF>,NE> AV = Evaluate_AV(kV,HfV,UV);
  
  // **************************************
  // quantities evaluated at inital point *
  // **************************************

  State s(E);
  s.r=rmin;
  s.update_potential(lnrho,temperature,Ye,P_unosc,HfV,s);
  for(state m=matter; m<=antimatter; m++){
    for(int i=0;i<=NE-1;i++){
      for(int j=0;j<=NF-1;j++){
	if(real(s.CC[m][i][j][mu][e]*CV[i][j][mu][e]) < 0.)
	  s.AA[m][i][j][e]=-AV[i][j][e];
	else s.AA[m][i][j][e]=AV[i][j][e];
	s.AA[m][i][j][mu]=AV[i][j][mu];
      }
      s.UU[m][i]=U(s.dkk[m][i],s.CC[m][i],s.AA[m][i]);
    }
  }
  initialize(s,rmin,E,D_unosc);
  const State s0 = s;
    
  // *****************************************
  // initialize at beginning of every domain *
  // *****************************************
  double dr_block = dr0;
  double dr_osc   = dr0;
  double dr_int   = dr0;
      
  // *************************************************
  // comment out if not following as a function of r *
  // *************************************************
	
  s.update_potential(lnrho,temperature,Ye,P_unosc,HfV,s0);
  FilePointers fp = setup_HDF5_file(s.E);
  write_data_HDF5(fp, s, dr_osc, dr_int, dr_block);
	
  // ***********************
  // start the loop over r *
  // ***********************
  bool finish = false;
  do{
    double r_end = s.r + dr_block * min(5., exponential_random());
    if(r_end>rmax){
      r_end = rmax;
      finish=true;
    }

    State sBlockStart = s;
    
    // oscillate
    if(do_oscillate){
      s.assert_noNaN(accuracy);
      s.r = sBlockStart.r;
      evolve_oscillations(s, s0, sBlockStart, r_end, dr_osc, lnrho, temperature, Ye, P_unosc, accuracy, increase, HfV);
    }

    // interact with the matter
    double impact = 0;
    if(do_interact){
      s.assert_noNaN(accuracy);
      s.r = sBlockStart.r;
      evolve_interactions(s, s0, sBlockStart, r_end, dr_int, lnrho, temperature, Ye, D_unosc, accuracy, increase, impact);
    }

    // output data
    write_data_HDF5(fp, s,dr_osc, dr_int, dr_block);
    
    // timestepping
    if(impact > target_impact)
      cout << "WARNING: impact="<<impact<< endl;
    double corrected_impact = impact / (s.r-sBlockStart.r) * dr_block;
    if(corrected_impact<.1*target_impact)
      dr_block *= min(increase, .1*target_impact/corrected_impact);
    if(corrected_impact>.1*target_impact)
      dr_block *= .1*target_impact/corrected_impact;
    dr_block = min(dr_block, dr_block_max);

  } while(finish==false);

  
  cout<<"\nFinished\n\a"; cout.flush();
  return 0;
}


