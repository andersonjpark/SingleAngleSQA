#ifndef _NULIB_INTERFACE_H
#define _NULIB_INTERFACE_H

#include "H5Cpp.h"
#include <string>
#include <cstdlib>

// physical constants
const double clight = 2.99792458e10; // cm/s
const double hplanck = 1.0545716e-27; // erg.s
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
		  const double*, //rho
		  const double*, //temp
		  const double*, //ye
		  double*, //eas_species_energy (3D array)
		  int*,    //number of species (3,5,6)
		  int*,    //number of groups
		  int*);   //number of easvariables (3)

  void nulibtable_single_species_range_energy_(
		  const double*, //rho
		  const double*, //temp
		  const double*, //Ye
		  int*,    //species number
		  double*, //eas_energy (2D array)
		  int*,    //number of groups
		  int*);   //number of easvariables(3)

  void nulibtable_epannihil_single_species_range_energy_(
		  const double* temp,  // MeV
		  const double* eta,   // mu/kT
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
  vector<double> eas;
  vector<double> escat_kernel0;
  vector<double> escat_kernel1;

  void resize(int ns_in, int ng_in, int nv_in, int read_Ielectron){
    ns = ns_in;
    ng = ng_in;
    nv = nv_in;
    eas.resize(ns*ng*nv);
    if(read_Ielectron){
      escat_kernel0.resize(ns*ng*ng*2);
      escat_kernel1.resize(ns*ng*ng*2);
    }
  }

  int index(int is,int ig,int iv){
    return is + ig*ns + iv*ns*ng;
  }

  int kernel_index(int is,int igin, int igout){
    return is + igin*ns + igout*ns*ng;
  }

  double emis(int is,int ig){ // 1/cm
    double Emid = __nulibtable_MOD_nulibtable_energies[ig]*MeV_to_ergs; // erg
    double dE3 = pow(__nulibtable_MOD_nulibtable_etop[ig],3) - pow(__nulibtable_MOD_nulibtable_ebottom[ig],3);
    dE3 *= pow(MeV_to_ergs,3);
    double tmp = hplanck*hplanck*hplanck * clight*clight /  (Emid*dE3/3.);//Emid*Emid * dE); // 
    double nulib_emis = eas[index(is,ig,0)]; // erg/ccm/s/sr
    return nulib_emis * tmp;
  }
  double abs(int is,int ig){ // 1/cm
    return eas[index(is,ig,1)];
  }
  double scat(int is,int ig){ // 1/cm
    return eas[index(is,ig,2)];
  }
  double delta(int is,int ig){ // 1/cm
    if(nv==4) return eas[index(is,ig,3)];
    else return 0;
  }
  double Bnu(int is, int ig){
    return emis(is,ig) / abs(is,ig);
  }
  double Phi0(int is,int igin, int igout){ // 1/cm
    double Phi0 = (igin==igout ? scat(is,igin) : 0);
    if(escat_kernel0.size() > 0)
      return Phi0 += escat_kernel0[kernel_index(is,igin,igout)];
    return Phi0;
  }
  double Phi1(int is,int igin, int igout){ // 1/cm
    if(escat_kernel1.size() > 0)
      return escat_kernel1[kernel_index(is,igin,igout)];
    else return 0;
  }

  static MATRIX<double,2,2> avg_matrix(double eval, double muval){
    MATRIX<double,2,2> result;
    result[e ][e ] = eval;
    result[mu][mu] = muval;
    result[e ][mu] = result[mu][e] = (eval + muval) / 2.;
    return result;
  }
  static MATRIX<double,2,2> tilde_matrix(double eval, double muval){
    const double sin2thetaW = 0.23122;
    MATRIX<double,2,2> result;
    result[e ][e ] = 0;
    result[mu][mu] = 0;
    result[e ][mu] = result[mu][e] = (eval - muval) / (4.*sin2thetaW);
    return result;
  }
  static MATRIX<complex<double>,2,2> blocking_term0(MATRIX<double,2,2> Phi0matrix, MATRIX<complex<double>,2,2> f, MATRIX<complex<double>,2,2> fp){
    MATRIX<complex<double>,2,2> result;
    for(flavour fa=e; fa<=mu; fa++)
      for(flavour fb=e; fb<=mu; fb++){
	result[fa][fb] = 0;
	for(flavour fc=e; fc<=mu; fc++){
	  result[fa][fb] += 0.5 * (Phi0matrix[fc][fb]*f[fa][fc]*fp[fc][fb] + Phi0matrix[fa][fc]*fp[fa][fc]*f[fc][fb]);
	}
      }
    return result;
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
	     __nulibtable_MOD_nulibtable_number_easvariables,
	     read_Ielectron);

}


#endif
