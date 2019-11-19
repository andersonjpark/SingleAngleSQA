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

#ifndef MISC_H
#define MISC_H

#include <sstream>
using std::stringstream;

template<typename T>
T get_parameter(ifstream& fin, const char* name){
  stringstream line;
  string linestring;
  T to_set;
  std::getline(fin, linestring);
  line = stringstream(linestring);
  line >> to_set;
  cout << "PARAMETER " << name << " = " << to_set << endl;
  return to_set;
}

inline double Sign(const double input){
  return input>0 ? 1 : -1;
}

// Cash-Karp RK parameters
const int NRK=6;
const int NRKOrder=5;
static const double AA[]={ 0., 1./5., 3./10., 3./5., 1., 7./8. };
static const double b0[]={};
static const double b1[]={ 1./5. };
static const double b2[]={ 3./40.,9./40. };
static const double b3[]={ 3./10.,-9./10.,6./5. };
static const double b4[]={ -11./54.,5./2.,-70./27.,35./27. };
static const double b5[]={ 1631./55296.,175./512.,575./13824.,44275./110592.,253./4096. };
static const double* BB[]={ b0,b1,b2,b3,b4,b5 };
static const double CC[]={ 37./378.,0.,250./621.,125./594.,0.,512./1771. };
static const double DD[]={ 2825./27648.,0.,18575./48384.,13525./55296.,277./14336.,1./4. };

// random number generation
#include <random>
double uniform(){
  return (float)rand() / (float)RAND_MAX;
}
double exponential_random(){
  return -log(uniform());
}

// energy grid stuff
double Ebottom(int i, const array<double,NE>& Etop){
  assert(i>=0);
  assert(i<NE);
  double result = (i>0 ? Etop[i-1] : 0);
  assert(result < Etop[i]);
  return result;
}

double Vphase(double Elow, double Ehi){ // cm^-3
  assert(Ehi>Elow);
  assert(Elow>=0);
  double dE3 = pow(Ehi,3) - pow(Elow,3);
  double result = 4.*M_PI * dE3/3. / pow(2.*M_PI*cgs::constants::hbarc,3);
  return result;
}

double Vphase(int i, const array<double,NE>& Etop){ // cm^-3
  assert(i>=0);
  assert(i<NE);
  return Vphase(Ebottom(i,Etop), Etop[i]);
}

double Vphase_overlap(double Elow1, double Ehi1, double Elow2, double Ehi2){ // cm^-3
  assert(Elow1<Ehi1);
  assert(Elow2<Ehi2);
  assert(Elow1>=0);
  assert(Elow2>=0);
  double Elow = max(Elow1,Elow2);
  double Ehi  = min(Ehi1, Ehi2);
  if(Ehi<=Elow) return 0;
  else return Vphase(Elow, Ehi);
}

#endif
