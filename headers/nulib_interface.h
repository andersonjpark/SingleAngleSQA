#ifndef _NULIB_INTERFACE_H
#define _NULIB_INTERFACE_H

#include "H5Cpp.h"
#include <string>
#include <cstdlib>
#include "mstl.h"

// physical constants
const double clight = 2.99792458e10;
const double hplanck = 1.0545716e-27;
const double MeV_to_ergs = 1.60217646e-6;

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

class EAS{
 public:
  int ns, ng, nv;
  vector<double> storage;

  void resize(int ns_in, int ng_in, int nv_in){
    ns = ns_in;
    ng = ng_in;
    nv = nv_in;
    storage.resize(ns*ng*nv);
  }

  int index(int is,int ig,int iv){
    return is + ig*ns + iv*ns*ng;
  }

  /* void fix_units(){ */
  /*   for(int ig=0; ig<ng; ig++){ */
  /*     double Emid = __nulibtable_MOD_nulibtable_energies[ig]*MeV_to_ergs; */
  /*     double dE = __nulibtable_MOD_nulibtable_ewidths[ig]*MeV_to_ergs; */
  /*     double tmp = hplanck*hplanck*hplanck * clight*clight /  (Emid*Emid*Emid * dE); */
  /*     for(int is=0; is<ns; is++) */
  /* 	storage[index(is,ig,0)] *= tmp; */
  /*   } */
  /* } */

  double emis(int is,int ig){
    return storage[index(is,ig,0)];
  }
  double abs(int is,int ig){
    return storage[index(is,ig,1)];
  }
  double scat(int is,int ig){
    return storage[index(is,ig,2)];
  }
  double Bnu(int is, int ig){
    return emis(is,ig) / abs(is,ig);
  }
};
EAS eas;


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

  eas.resize(__nulibtable_MOD_nulibtable_number_species,
	     __nulibtable_MOD_nulibtable_number_groups,
	     __nulibtable_MOD_nulibtable_number_easvariables);

}

void initialize(vector<vector<MATRIX<complex<double>,NF,NF> > >& fmatrixf,
	      double rho, double T, double Ye){

  // don't do anything if too sparse
  if(log10(rho) <= __nulibtable_MOD_nulibtable_logrho_min)
    abort();

  // T should be MeV
  // nspecies, ngroups, nvars
  nulibtable_range_species_range_energy_(&rho, &T, &Ye, &eas.storage.front(),
					 &__nulibtable_MOD_nulibtable_number_species,
					 &__nulibtable_MOD_nulibtable_number_groups,
					 &__nulibtable_MOD_nulibtable_number_easvariables);
  //eas.fix_units();

  for(int i=0; i<NE; i++){
    for(state m=matter; m<=antimatter; m++)
      for(flavour f1=e; f1<=mu; f1++)
	for(flavour f2=e; f2<=mu; f2++) 
	  fmatrixf[m][i][f1][f2] = 0;

    fmatrixf[    matter][i][e ][e ] = eas.Bnu(0,i);
    fmatrixf[    matter][i][mu][mu] = eas.Bnu(2,i);
    fmatrixf[antimatter][i][e ][e ] = eas.Bnu(1,i);
    fmatrixf[antimatter][i][mu][mu] = eas.Bnu(2,i);
  }
}
void interact(vector<vector<MATRIX<complex<double>,NF,NF> > >& fmatrixf,
	      double rho, double T, double Ye, double dr){

  // don't do anything if too sparse
  if(log10(rho) <= __nulibtable_MOD_nulibtable_logrho_min)
    return;

  // T should be MeV
  // nspecies, ngroups, nvars
  nulibtable_range_species_range_energy_(&rho, &T, &Ye, &eas.storage.front(),
					 &__nulibtable_MOD_nulibtable_number_species,
					 &__nulibtable_MOD_nulibtable_number_groups,
					 &__nulibtable_MOD_nulibtable_number_easvariables);
  //eas.fix_units();


  double tmp = 0;
  for(int i=0; i<NE; i++){
    //cout << i << " " << eas.emis(0,i) << " " << eas.abs(0,i) << endl;
    // scale the off-diagonal components
    tmp = (eas.abs(0,i) + eas.abs(2,i)) * 0.5;
    fmatrixf[    matter][i][e][mu] *= exp(-tmp * dr);
    fmatrixf[    matter][i][mu][e] = fmatrixf[matter][i][e][mu];

    tmp = (eas.abs(1,i) + eas.abs(2,i)) * 0.5;
    fmatrixf[antimatter][i][e][mu] *= exp(-tmp * dr);
    fmatrixf[antimatter][i][mu][e] = fmatrixf[antimatter][i][e][mu];

    // scale the diagonal components
    fmatrixf[    matter][i][e ][e ] += (eas.emis(0,i) - eas.abs(0,i)*fmatrixf[    matter][i][e ][e ]) * dr;
    fmatrixf[    matter][i][mu][mu] += (eas.emis(2,i) - eas.abs(2,i)*fmatrixf[    matter][i][mu][mu]) * dr;
    fmatrixf[antimatter][i][e ][e ] += (eas.emis(1,i) - eas.abs(1,i)*fmatrixf[antimatter][i][e ][e ]) * dr;
    fmatrixf[antimatter][i][mu][mu] += (eas.emis(2,i) - eas.abs(2,1)*fmatrixf[antimatter][i][mu][mu]) * dr;
  }
}

#endif
