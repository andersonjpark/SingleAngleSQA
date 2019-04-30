#ifndef _NULIB_INTERFACE_H
#define _NULIB_INTERFACE_H

#include "H5Cpp.h"
#include <string>
#include <cstdlib>

// physical constants
const double clight = 2.99792458e10; // cm/s
const double hplanck = 1.0545716e-27; // erg.s
const double MeV_to_ergs = 1.60217646e-6;
const int NEUTRINO_SCHEME = 2;
const long int NS_NULIB = 4;

inline bool hdf5_dataset_exists(const char* filename, const char* datasetname){
  bool exists = true;

  // Temporarily turn off error printing                                                                                   
  H5E_auto2_t func;
  void* client_data;
  H5::Exception::getAutoPrint(func,&client_data);
  H5::Exception::dontPrint();

  // See if dataset exists                                                                                                 
  H5::H5File file(filename, H5F_ACC_RDONLY);
  H5::DataSet dataset;
  try{
    dataset = file.openDataSet(datasetname);
  }
  catch(H5::FileIException& exception){
    exists = false;
  }

  // Turn error printing back on                                                                                           
  H5::Exception::setAutoPrint(func,client_data);
  file.close();

  return exists;
}

// module variables set in fortran NuLib code
extern int     __nulib_MOD_total_eos_variables;
extern double* __nulib_MOD_energies;
extern "C"{
  void set_eos_variables_(double* eos_variables);
  void read_eos_table_(char* filename);
  void total_absorption_opacities_(const int* neutrino_species,
				   const double* neutrino_energy,
				   double* absorption_opacity,
				   const double* eos_variables);
}


class EAS{
 public:

  vector<double> eos_variables;
  double munue_kT, eta;

  array<array<double,NE>,NS_NULIB> absorption_opacity, emissivities, scattering_opacity, delta;
  
  EAS(){
    int number_species = 6; // has to be 6 no matter how many are included in nux
    int number_groups = NE;
    eos_variables.resize(__nulib_MOD_total_eos_variables);
  }

  void update(const double rho /* g/ccm */, const double T /*MeV*/, const double Ye){
    // set the EOS variables
    for(int i=0; i<__nulib_MOD_total_eos_variables; i++) eos_variables[i] = 0;
    eos_variables[0] = rho;
    eos_variables[1] = T;
    eos_variables[2] = Ye;
    set_eos_variables_(&eos_variables[0]);
    double mue = eos_variables[10];
    double muhat = eos_variables[13];
    munue_kT = (mue-muhat) / T;
    eta = mue / T;
  }

  /* int index(const int is,const int ig,const int iv) const{ */
  /*   return is + ig*ns + iv*ns*ng; */
  /* } */

  /* int kernel_index(const int is,const int igin, const int igout) const{ */
  /*   return is + igin*ns + igout*ns*ng; */
  /* } */

  /* inline int nu4_kernel_index(const int ik, const int i1, const int i3) const{ */
  /*   return i3 + ng*i1 + ng*ng*ik; */
  /* } */
  /* inline int nu4_bin2(const int ik, const int i1, const int i3) const{ */
  /*   return i1+i3-ik; */
  /* } */
  
  double abs(int s, double E /*erg*/) const{ // 1/cm
    double EMeV = E / (1e6*eV);
    int s_nulib = s+1;

    double absopac;
    total_absorption_opacities_(&s_nulib, &EMeV, &absopac, &eos_variables[0]);

    return absopac;
  }
  double scat(int s, int ig) const{ // 1/cm
    return 1;
  }
  /* double delta(const int is,const int ig) const{ // 1/cm */
  /*   if(do_delta){ */
  /*     int ind = index(is,ig,3); */
  /*     assert(ind < eas.size()); */
  /*     return eas[ind]; */
  /*   } */
  /*   else return 0; */
  /* } */
  double fermidirac(int s, double E_kT) const{
    double mu_kT;
    if     (s==0) mu_kT =  munue_kT;
    else if(s==1) mu_kT = -munue_kT;
    else          mu_kT = 0;
    double result = 1./(1. + exp(E_kT-mu_kT));
    return result;
  }
  /* double Phi0scat(const int is,const int igin, const int igout) const{ // cm^3/s/sr */
  /*   double result = 0; */
  /*   if(igin == igout) */
  /*     result += scat(is,igin) */
  /* 	/(4.*M_PI*nu[igin]*nu[igin]*dnu[igin]/cgs::constants::c4); */
  /*   if(do_iscat) */
  /*     result += escat_kernel0[kernel_index(is,igin,igout)]; */
  /*   return result; */
  /* } */
  /* double Phi1scat(const int is,const int igin, const int igout) const{ // cm^3/s/sr */
  /*   double result = 0; */
  /*   if(igin == igout) */
  /*     result += scat(is,igin)*delta(is,igin)/3. */
  /* 	/(4.*M_PI*nu[igin]*nu[igin]*dnu[igin]/cgs::constants::c4); */
  /*   if(do_iscat) */
  /*     result += escat_kernel1[kernel_index(is,igin,igout)]; */
  /*   else return 0; */
  /* } */
  /* double Phi0pair(const int is,const int igin, const int igout) const{ // cm^3/s/sr */
  /*   double result = 0; */
  /*   if(do_pair) */
  /*     result += pair_kernel0[kernel_index(is,igin,igout)]; */
  /*   return result; */
  /* } */
  /* double Phi1pair(const int is,const int igin, const int igout) const{ // cm^3/s/sr */
  /*   double result = 0; */
  /*   if(do_pair) */
  /*     result += pair_kernel1[kernel_index(is,igin,igout)]; */
  /*   else return 0; */
  /* } */

};
EAS eas;

#endif
