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
#include "nulib_interface.h"
#include "evolve.h"
#include "profile.h"

//======//
// MAIN //
//======//
int main(int argc, char *argv[]){
  string inputfilename;

  assert(argc==2);
  inputfilename=string(argv[1]);
  ifstream fin(inputfilename.c_str());

  // load the nulib table
  const string input_directory = get_parameter<string>(fin,"input_directory");
  const string eosfilename = get_parameter<string>(fin, "eosfilename");
  const double rhostart = get_parameter<double>(fin, "rhostart"); // cm
  const double dr0 = get_parameter<double>(fin, "dr0"); // cm
  const double dr_block_max = get_parameter<double>(fin, "dr_block_max"); // cm
  const double accuracy = get_parameter<double>(fin, "accuracy");
  const bool do_oscillate = get_parameter<bool>(fin, "do_oscillate");
  const bool do_interact = get_parameter<bool>(fin, "do_interact");
  const bool do_interact_rotation = get_parameter<bool>(fin, "do_interact_rotation");
  const bool do_SR = get_parameter<bool>(fin, "do_SR");
  const bool do_GR = get_parameter<bool>(fin, "do_GR");
  const bool do_two_loop_contribution = get_parameter<bool>(fin, "do_two_loop_contribution");
  const double target_impact = get_parameter<double>(fin, "target_impact");
  const double increase = get_parameter<double>(fin, "increase");
  const double initial_mixing = get_parameter<double>(fin, "initial_mixing");
  assert(do_oscillate or do_interact);
  fin.close();

  Profile profile(input_directory, rhostart, do_SR, do_GR);
  read_eos_table_((char*)eosfilename.c_str());
  cout << "m_ref = " << __nulib_MOD_m_ref << " (939=m_n for LS, 931=m_amu for Hempel)" << endl;

  // **************************************
  // quantities evaluated at inital point *
  // **************************************

  State s(profile, profile.rho.x[0], initial_mixing);
  State s0 = s; // ONLY used for oscillation stuff. s0.fmatrixf is meaningless

  // *****************************************
  // initialize at beginning of every domain *
  // *****************************************
  double dr_block = dr0;
  double dr_osc   = dr0;
  double dr_int   = dr0;

  // set up files
  string outputfilename = "output.h5";
  ifstream tmp_ifstream(outputfilename);
  FilePointers fp;
  if(tmp_ifstream){
    cout << "Recovering from " << outputfilename << endl;
    recover(outputfilename, s, dr_osc, dr_int, dr_block, fp);
  }
  else{
    cout << "Creating " << outputfilename << endl;
    setup_HDF5_file(outputfilename, profile.Ecom, profile.Etopcom, fp);
    write_data_HDF5(fp, s, dr_osc, dr_int, dr_block, true);
  }

  // ***********************
  // start the loop over r *
  // ***********************
  bool finish = false;
  size_t iter = 0;
  do{
    double r_end = s.r + dr_block * min(5., exponential_random());
    if(r_end>profile.rho.XMax()){
      r_end = profile.rho.XMax();
      finish=true;
    }
    if(s.r<0 and r_end>=0){
      r_end = 0;
      dr_block = dr0;
      dr_osc   = dr0;
      dr_int   = dr0;
    }

    State sBlockStart = s;

    // oscillate
    if(do_oscillate){
      s.assert_noNaN(accuracy);
      s.r = sBlockStart.r;
      evolve_oscillations(s, s0, sBlockStart, r_end, dr_osc, profile, accuracy, increase);
    }

    // interact with the matter
    double impact = 0;
    if(do_interact){
      s.assert_noNaN(accuracy);
      s.r = sBlockStart.r;
      evolve_interactions(s, sBlockStart, r_end, dr_int, profile, accuracy, increase, do_interact_rotation, impact);
      s.assert_noNaN(accuracy);
    }

    // output data
    bool do_header = (++iter)%100==0;
    write_data_HDF5(fp, s,dr_osc, dr_int, dr_block, do_header);

    // make sure fmatrixf is Hermitian and Scumulative is unitary
    #pragma omp parallel for collapse(2)
    for(int i=0; i<NE; i++){
      for(int m=matter; m<=antimatter; m++){
	Hermitize(s.fmatrixf[m][i], accuracy);
	unitarize(s.Scumulative[m][i], accuracy);
      }
    }

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
