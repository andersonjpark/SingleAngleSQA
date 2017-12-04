#ifndef _NULIB_INTERFACE_H
#define _NULIB_INTERFACE_H

#include "H5Cpp.h"
#include <string>
#include <cstdlib>

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
extern int     __nulibtable_MOD_nulibtable_number_species;
extern int     __nulibtable_MOD_nulibtable_number_easvariables;
extern int     __nulibtable_MOD_nulibtable_number_groups;
extern int     __nulibtable_MOD_nulibtable_nrho;
extern int     __nulibtable_MOD_nulibtable_ntemp;
extern int     __nulibtable_MOD_nulibtable_nye;
extern int     __nulibtable_MOD_nulibtable_nitemp;
extern int     __nulibtable_MOD_nulibtable_nieta;
extern double* __nulibtable_MOD_nulibtable_energies;
extern double* __nulibtable_MOD_nulibtable_ewidths;
extern double* __nulibtable_MOD_nulibtable_ebottom;
extern double* __nulibtable_MOD_nulibtable_etop;
extern double* __nulibtable_MOD_nulibtable_logrho;
extern double* __nulibtable_MOD_nulibtable_logtemp;
extern double* __nulibtable_MOD_nulibtable_ye;
extern double* __nulibtable_MOD_nulibtable_logitemp;
extern double* __nulibtable_MOD_nulibtable_logieta;
extern double  __nulibtable_MOD_nulibtable_logtemp_min;
extern double  __nulibtable_MOD_nulibtable_logtemp_max;
extern double  __nulibtable_MOD_nulibtable_logrho_min;
extern double  __nulibtable_MOD_nulibtable_logrho_max;
extern double  __nulibtable_MOD_nulibtable_ye_min;
extern double  __nulibtable_MOD_nulibtable_ye_max;
extern double  __nulibtable_MOD_nulibtable_logitemp_min;
extern double  __nulibtable_MOD_nulibtable_logitemp_max;
extern double  __nulibtable_MOD_nulibtable_logieta_min;
extern double  __nulibtable_MOD_nulibtable_logieta_max;
extern int     __nulib_MOD_total_eos_variables;

// These are fortran functions and module variables in nulib.a                                                                   
extern "C"{
  void nulibtable_range_species_range_energy_(
		  double*, //rho
		  double*, //temp
		  double*, //ye
		  double*, //eas_species_energy (3D array)
		  int*,    //number of species (3,5,6)
		  int*,    //number of groups
		  int*);   //number of easvariables (3)

  void nulibtable_single_species_range_energy_(
		  double*, //rho
		  double*, //temp
		  double*, //Ye
		  int*,    //species number
		  double*, //eas_energy (2D array)
		  int*,    //number of groups
		  int*);   //number of easvariables(3)

  void nulibtable_epannihil_single_species_range_energy_(
		  double* temp,  // MeV
		  double* eta,   // mu/kT
		  int* lns,      // species number
		  double* phi,   // phi[legendre-p/a index][this_group][anti-group]
		  int* ngroups1,
		  int* ngroups2,
		  int* n_phis);

  void nulibtable_inelastic_single_species_range_energy_(
		  double* temp,  // MeV
		  double* eta,   // mu/kT
		  int* lns,      // species number
		  double* phi,   // phi[legendre index][out group][in group]
		  int* ngroups1, // ng in
		  int* ngroups2, // ng out (should be same as eas_n1)
		  int* n_phis);   // number of legendre terms (=2)

  void nulibtable_reader_(char*,int*,int*,int*,int);
}

/**************/
/* nulib_init */
/**************/
void nulib_init(string filename, int use_scattering_kernels){
  int read_Ielectron = 0;
  int read_epannihil = 0;
  int read_delta = 0;
  if(hdf5_dataset_exists(filename.c_str(),"/scattering_delta")) read_delta = 1;
  if(hdf5_dataset_exists(filename.c_str(),"/inelastic_phi0"))   read_Ielectron = 1;  
  nulibtable_reader_((char*)filename.c_str(), &read_Ielectron, &read_epannihil, &read_delta, filename.length());
   
  cout << "#   logrho range: {" << __nulibtable_MOD_nulibtable_logrho_min << "," << __nulibtable_MOD_nulibtable_logrho_max << "} g/ccm" << endl;
  cout << "#   logT   range: {" << __nulibtable_MOD_nulibtable_logtemp_min << "," << __nulibtable_MOD_nulibtable_logtemp_max << "} MeV" << endl;
  cout << "#   Ye  range: {" << __nulibtable_MOD_nulibtable_ye_min << "," << __nulibtable_MOD_nulibtable_ye_max << "}" << endl;
  cout << "#   E   range: {" << __nulibtable_MOD_nulibtable_ebottom[0] << "," << __nulibtable_MOD_nulibtable_etop[__nulibtable_MOD_nulibtable_number_groups-1] << "} MeV" << endl;
  cout << "#   n_species = " << __nulibtable_MOD_nulibtable_number_species << endl;
  cout << "#   n_rho   = " << __nulibtable_MOD_nulibtable_nrho << endl;
  cout << "#   n_T     = " << __nulibtable_MOD_nulibtable_ntemp << endl;
  cout << "#   n_Ye    = " << __nulibtable_MOD_nulibtable_nye << endl;
  cout << "#   n_E     = " << __nulibtable_MOD_nulibtable_number_groups << endl;
}


#endif
