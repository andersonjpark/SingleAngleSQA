#ifndef _NULIB_INTERFACE_H
#define _NULIB_INTERFACE_H

#include "H5Cpp.h"
#include <string>
#include <cstdlib>
#include <cmath>

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
extern double __nulib_MOD_m_ref, __eosmodule_MOD_eos_rhomin, __eosmodule_MOD_eos_yemin, __eosmodule_MOD_eos_tempmin;
extern const bool __nulib_MOD_do_weak_mag_corrections,
__nulib_MOD_do_ionion_correlation,
__nulib_MOD_do_heavyscat_formfactor,
__nulib_MOD_do_electronpolarization_correction,
__nulib_MOD_do_nc_virial_correction,
__nulib_MOD_do_strange_coupling,
__nulib_MOD_do_transport_opacities,
__nulib_MOD_add_nue_absorption_on_n,
__nulib_MOD_add_anue_absorption_on_p,
__nulib_MOD_add_nue_absorption_on_a,
__nulib_MOD_add_nux_absorption_on_n_and_p,
__nulib_MOD_add_nue_scattering_n,
__nulib_MOD_add_nue_scattering_p,
__nulib_MOD_add_nue_scattering_heavies,
__nulib_MOD_add_nue_scattering_electrons,
__nulib_MOD_add_nue_scattering_alphas,
__nulib_MOD_add_anue_scattering_n,
__nulib_MOD_add_anue_scattering_p,
__nulib_MOD_add_anue_scattering_heavies,
__nulib_MOD_add_anue_scattering_electrons,
__nulib_MOD_add_anue_scattering_alphas,
__nulib_MOD_add_numu_scattering_n,
__nulib_MOD_add_numu_scattering_p,
__nulib_MOD_add_numu_scattering_heavies,
__nulib_MOD_add_numu_scattering_electrons,
__nulib_MOD_add_numu_scattering_alphas,
__nulib_MOD_add_anumu_scattering_n,
__nulib_MOD_add_anumu_scattering_p,
__nulib_MOD_add_anumu_scattering_heavies,
__nulib_MOD_add_anumu_scattering_electrons,
__nulib_MOD_add_anumu_scattering_alphas,
__nulib_MOD_add_nutau_scattering_n,
__nulib_MOD_add_nutau_scattering_p,
__nulib_MOD_add_nutau_scattering_heavies,
__nulib_MOD_add_nutau_scattering_electrons,
__nulib_MOD_add_nutau_scattering_alphas,
__nulib_MOD_add_anutau_scattering_n,
__nulib_MOD_add_anutau_scattering_p,
__nulib_MOD_add_anutau_scattering_heavies,
__nulib_MOD_add_anutau_scattering_electrons,
__nulib_MOD_add_anutau_scattering_alphas,
__nulib_MOD_add_nue_iscattering_electrons,
__nulib_MOD_add_anue_iscattering_electrons,
__nulib_MOD_add_numu_iscattering_electrons,
__nulib_MOD_add_anumu_iscattering_electrons,
__nulib_MOD_add_nutau_iscattering_electrons,
__nulib_MOD_add_anutau_iscattering_electrons,
__nulib_MOD_add_nue_emission_epannihil,
__nulib_MOD_add_anue_emission_epannihil,
__nulib_MOD_add_numu_emission_epannihil,
__nulib_MOD_add_anumu_emission_epannihil,
__nulib_MOD_add_nutau_emission_epannihil,
__nulib_MOD_add_anutau_emission_epannihil,
__nulib_MOD_add_nue_emission_bremsstrahlung,
__nulib_MOD_add_anue_emission_bremsstrahlung,
__nulib_MOD_add_numu_emission_bremsstrahlung,
__nulib_MOD_add_anumu_emission_bremsstrahlung,
__nulib_MOD_add_nutau_emission_bremsstrahlung,
__nulib_MOD_add_anutau_emission_bremsstrahlung,
__nulib_MOD_add_nue_emission_weakinteraction_ecap,
__nulib_MOD_add_anue_emission_weakinteraction_poscap,
__nulib_MOD_apply_kirchoff_to_pair_creation,
__nulib_MOD_add_nue_kernel_epannihil,
__nulib_MOD_add_anue_kernel_epannihil,
__nulib_MOD_add_numu_kernel_epannihil,
__nulib_MOD_add_anumu_kernel_epannihil,
__nulib_MOD_add_nutau_kernel_epannihil,
__nulib_MOD_add_anutau_kernel_epannihil,
__nulib_MOD_add_anutau_kernel_epannihil,
__nulib_MOD_add_nu4pair_kernel,
__nulib_MOD_add_nu4scat_kernel,
__nulib_MOD_add_nue_kernel_bremsstrahlung,
__nulib_MOD_add_anue_kernel_bremsstrahlung,
__nulib_MOD_add_numu_kernel_bremsstrahlung,
__nulib_MOD_add_anumu_kernel_bremsstrahlung,
__nulib_MOD_add_nutau_kernel_bremsstrahlung,
__nulib_MOD_add_anutau_kernel_bremsstrahlung;

extern "C"{
	void __nulib_MOD_initialize_nulib(const int* neutrino_scheme,
			const int* number_species,
			const int* number_groups);
	void set_eos_variables_(double* eos_variables);
	void read_eos_table_(char* filename);
	void total_absorption_opacities_(const int* neutrino_species,
			const double* neutrino_energy,
			double* absorption_opacity,
			const double* eos_variables);
	void total_scattering_opacity_(const int* neutrino_species,
			const double* neutrino_energy,
			double* scattering_opacity,
			double* delta,
			const double* eos_variables);

	// cm^3/s
	double nes_phi0_thompsonbruenn_(const double* nu_energy_in /*MeV*/,
			const double* nu_energy_out /*MeV*/,
			const double* matter_eta, /*mue/T*/
			const double* matter_temperature, /*MeV*/
			const int* neutrino_species);
	double nes_phi1_thompsonbruenn_(const double* nu_energy_in /*MeV*/,
			const double* nu_energy_out /*MeV*/,
			const double* matter_eta, /*mue/T*/
			const double* matter_temperature, /*MeV*/
			const int* neutrino_species);
	double epannihil_phi_bruenn_(const double* nu_energy_x, /*E/kT*/
			const double* nubar_energy_x, /*E/kT*/
			const double* matter_eta, /*mue/kT*/
			const int* neutrino_species,
			const int* pro_or_ann,
			const int* which_l);
	double bremsstrahlung_phi0_hannestad_(const double* nu_energy_x, /*E/kT*/
			const double* nubar_energy_x, /*E/kT*/
			const double* matter_temperature, /*MeV*/
			const double* n_N, /*cm^-3*/
			const int* neutrino_species,
			const int* pro_ann);

	// 1/cm, all arguments in MeV
	double __nulib_MOD_nu4scat_kernel_single(double* k, double* q1, double* q2, double* q3); // cm^5
	double __nulib_MOD_nu4pair_kernel_single(double* k, double* q1, double* q2, double* q3); // cm^5
}


class EAS{
public:

	vector<double> eos_variables;
	double munue_kT, eta;

	array<array<double,NE>,NS_NULIB> absorption_opacity, emissivities, scattering_opacity, delta;

	EAS(){
		eos_variables.resize(__nulib_MOD_total_eos_variables);
		eta = 0;
		munue_kT = 0;

		int neutrino_scheme = 2; // e,ebar,x,xbar
		int number_species = 6;
		__nulib_MOD_initialize_nulib(&neutrino_scheme, &number_species, &NE);

		// output all nulib parameters
		cout << endl;
		cout << "REQUESTED INTERACTIONS" << endl;
		cout << __nulib_MOD_do_weak_mag_corrections << " do_weak_mag_corrections" << endl;
		cout << __nulib_MOD_do_ionion_correlation << " do_ionion_correlation" << endl;
		cout << __nulib_MOD_do_heavyscat_formfactor << " do_heavyscat_formfactor" << endl;
		cout << __nulib_MOD_do_electronpolarization_correction << " do_electronpolarization_correction" << endl;
		cout << __nulib_MOD_do_nc_virial_correction << " do_nc_virial_correction" << endl;
		cout << __nulib_MOD_do_strange_coupling << " do_strange_coupling" << endl;
		cout << __nulib_MOD_do_transport_opacities << " do_transport_opacities" << endl;
		cout << __nulib_MOD_add_nue_absorption_on_n << " add_nue_absorption_on_n" << endl;
		cout << __nulib_MOD_add_anue_absorption_on_p << " add_anue_absorption_on_p" << endl;
		cout << __nulib_MOD_add_nue_absorption_on_a << " add_nue_absorption_on_A" << endl;
		cout << __nulib_MOD_add_nux_absorption_on_n_and_p << " add_nux_absorption_on_n_and_p" << endl;
		cout << __nulib_MOD_add_nue_scattering_n << " add_nue_scattering_n" << endl;
		cout << __nulib_MOD_add_nue_scattering_p << " add_nue_scattering_p" << endl;
		cout << __nulib_MOD_add_nue_scattering_heavies << " add_nue_scattering_heavies" << endl;
		cout << __nulib_MOD_add_nue_scattering_electrons << " add_nue_scattering_electrons" << endl;
		cout << __nulib_MOD_add_nue_scattering_alphas << " add_nue_scattering_alphas" << endl;
		cout << __nulib_MOD_add_anue_scattering_n << " add_anue_scattering_n" << endl;
		cout << __nulib_MOD_add_anue_scattering_p << " add_anue_scattering_p" << endl;
		cout << __nulib_MOD_add_anue_scattering_heavies << " add_anue_scattering_heavies" << endl;
		cout << __nulib_MOD_add_anue_scattering_electrons << " add_anue_scattering_electrons" << endl;
		cout << __nulib_MOD_add_anue_scattering_alphas << " add_anue_scattering_alphas" << endl;
		cout << __nulib_MOD_add_numu_scattering_n << " add_numu_scattering_n" << endl;
		cout << __nulib_MOD_add_numu_scattering_p << " add_numu_scattering_p" << endl;
		cout << __nulib_MOD_add_numu_scattering_heavies << " add_numu_scattering_heavies" << endl;
		cout << __nulib_MOD_add_numu_scattering_electrons << " add_numu_scattering_electrons" << endl;
		cout << __nulib_MOD_add_numu_scattering_alphas << " add_numu_scattering_alphas" << endl;
		cout << __nulib_MOD_add_anumu_scattering_n << " add_anumu_scattering_n" << endl;
		cout << __nulib_MOD_add_anumu_scattering_p << " add_anumu_scattering_p" << endl;
		cout << __nulib_MOD_add_anumu_scattering_heavies << " add_anumu_scattering_heavies" << endl;
		cout << __nulib_MOD_add_anumu_scattering_electrons << " add_anumu_scattering_electrons" << endl;
		cout << __nulib_MOD_add_anumu_scattering_alphas << " add_anumu_scattering_alphas" << endl;
		cout << __nulib_MOD_add_nutau_scattering_n << " add_nutau_scattering_n" << endl;
		cout << __nulib_MOD_add_nutau_scattering_p << " add_nutau_scattering_p" << endl;
		cout << __nulib_MOD_add_nutau_scattering_heavies << " add_nutau_scattering_heavies" << endl;
		cout << __nulib_MOD_add_nutau_scattering_electrons << " add_nutau_scattering_electrons" << endl;
		cout << __nulib_MOD_add_nutau_scattering_alphas << " add_nutau_scattering_alphas" << endl;
		cout << __nulib_MOD_add_anutau_scattering_n << " add_anutau_scattering_n" << endl;
		cout << __nulib_MOD_add_anutau_scattering_p << " add_anutau_scattering_p" << endl;
		cout << __nulib_MOD_add_anutau_scattering_heavies << " add_anutau_scattering_heavies" << endl;
		cout << __nulib_MOD_add_anutau_scattering_electrons << " add_anutau_scattering_electrons" << endl;
		cout << __nulib_MOD_add_anutau_scattering_alphas << " add_anutau_scattering_alphas" << endl;
		cout << __nulib_MOD_add_nue_iscattering_electrons << " add_nue_Iscattering_electrons" << endl;
		cout << __nulib_MOD_add_anue_iscattering_electrons << " add_anue_Iscattering_electrons" << endl;
		cout << __nulib_MOD_add_numu_iscattering_electrons << " add_numu_Iscattering_electrons" << endl;
		cout << __nulib_MOD_add_anumu_iscattering_electrons << " add_anumu_Iscattering_electrons" << endl;
		cout << __nulib_MOD_add_nutau_iscattering_electrons << " add_nutau_Iscattering_electrons" << endl;
		cout << __nulib_MOD_add_anutau_iscattering_electrons << " add_anutau_Iscattering_electrons" << endl;
		cout << __nulib_MOD_add_nue_emission_epannihil << " add_nue_emission_epannihil" << endl;
		cout << __nulib_MOD_add_anue_emission_epannihil << " add_anue_emission_epannihil" << endl;
		cout << __nulib_MOD_add_numu_emission_epannihil << " add_numu_emission_epannihil" << endl;
		cout << __nulib_MOD_add_anumu_emission_epannihil << " add_anumu_emission_epannihil" << endl;
		cout << __nulib_MOD_add_nutau_emission_epannihil << " add_nutau_emission_epannihil" << endl;
		cout << __nulib_MOD_add_anutau_emission_epannihil << " add_anutau_emission_epannihil" << endl;
		cout << __nulib_MOD_add_nue_emission_bremsstrahlung << " add_nue_emission_bremsstrahlung" << endl;
		cout << __nulib_MOD_add_anue_emission_bremsstrahlung << " add_anue_emission_bremsstrahlung" << endl;
		cout << __nulib_MOD_add_numu_emission_bremsstrahlung << " add_numu_emission_bremsstrahlung" << endl;
		cout << __nulib_MOD_add_anumu_emission_bremsstrahlung << " add_anumu_emission_bremsstrahlung" << endl;
		cout << __nulib_MOD_add_nutau_emission_bremsstrahlung << " add_nutau_emission_bremsstrahlung" << endl;
		cout << __nulib_MOD_add_anutau_emission_bremsstrahlung << " add_anutau_emission_bremsstrahlung" << endl;
		cout << __nulib_MOD_add_nue_emission_weakinteraction_ecap << " add_nue_emission_weakinteraction_ecap" << endl;
		cout << __nulib_MOD_add_anue_emission_weakinteraction_poscap << " add_anue_emission_weakinteraction_poscap" << endl;
		cout << __nulib_MOD_apply_kirchoff_to_pair_creation << " apply_kirchoff_to_pair_creation" << endl;
		cout << __nulib_MOD_add_nue_kernel_epannihil << " add_nue_kernel_epannihil" << endl;
		cout << __nulib_MOD_add_anue_kernel_epannihil << " add_anue_kernel_epannihil" << endl;
		cout << __nulib_MOD_add_numu_kernel_epannihil << " add_numu_kernel_epannihil" << endl;
		cout << __nulib_MOD_add_anumu_kernel_epannihil << " add_anumu_kernel_epannihil" << endl;
		cout << __nulib_MOD_add_nutau_kernel_epannihil << " add_nutau_kernel_epannihil" << endl;
		cout << __nulib_MOD_add_anutau_kernel_epannihil << " add_anutau_kernel_epannihil" << endl;
		cout << endl;
	}

	void update(const double rho /* g/ccm */, const double T /*MeV*/, const double Ye){
		// set the EOS variables
		for(int i=0; i<__nulib_MOD_total_eos_variables; i++) eos_variables[i] = 0;
		eos_variables[0] = max(rho,__eosmodule_MOD_eos_rhomin);
		eos_variables[1] = max(T,  __eosmodule_MOD_eos_tempmin);
		eos_variables[2] = max(Ye, __eosmodule_MOD_eos_yemin);
		set_eos_variables_(&eos_variables[0]);
		double mue = eos_variables[10];
		double muhat = eos_variables[13];
		munue_kT = (mue-muhat) / T;
		eta = mue / T;
	}

	double abs(int s, double E /*erg*/) const{ // 1/cm
		double EMeV = E / (1e6*eV);
		int s_nulib = s+1;

		double absopac;
		total_absorption_opacities_(&s_nulib, &EMeV, &absopac, &eos_variables[0]);

		return absopac;
	}

	double fermidirac(int s, double E_kT) const{
		double mu_kT;
		if     (s==0) mu_kT =  munue_kT;
		else if(s==1) mu_kT = -munue_kT;
		else          mu_kT = 0;
		double result = 1./(1. + exp(E_kT-mu_kT));
		return result;
	}

	//===================//
	// SCATTERING KERNEL //
	//===================//
	// returns out-scattering rate from Ein to Eout.
	// 1/cm
	array<double,KMOMENTS> escat_opac(int s, double T /*MeV*/, double Ein /*erg*/){
		array<double,KMOMENTS> opac;
		const double EinMeV  = Ein  / (1e6*eV);
		int s_nulib = s+1;
		double scatopac, delta;
		total_scattering_opacity_(&s_nulib, &EinMeV, &scatopac, &delta, &eos_variables[0]);

		opac[0] = scatopac;
		opac[1] = scatopac *delta;

		assert(opac[0]>=0);
		assert(fabs(opac[1])<=opac[0]);
		for(int i=0; i<KMOMENTS; i++){
			assert(opac[i]==opac[i]);
			assert(fabs(opac[i]) < std::numeric_limits<double>::infinity() );
		}
		return opac;
	}


	array<double,KMOMENTS> Phi_iscat(int s, double T /*MeV*/,
			double Ein /*erg*/, double Eout /*erg*/,
			double Vout /*cm^-3*/){
		array<double,KMOMENTS> Phi;
		const double EinMeV  = Ein  / (1e6*eV);
		const double EoutMeV = Eout / (1e6*eV);
		int s_nulib = s+1;

		// elastic scattering contribution
		Phi[0] = 0;
		Phi[1] = 0;

		// inelastic scattering contribution
		// only use NuLib to calculate rate for Ein>=Eout as per NuLib suggestion
		// convert to rate with swapped arguments using Bruenn85 eq. C9
		if(__nulib_MOD_add_nue_iscattering_electrons){
			assert(__nulib_MOD_add_anue_iscattering_electrons);
			assert(__nulib_MOD_add_numu_iscattering_electrons);
			assert(__nulib_MOD_add_anumu_iscattering_electrons);
			assert(__nulib_MOD_add_nutau_iscattering_electrons);
			assert(__nulib_MOD_add_anutau_iscattering_electrons);
			const double Elo = min(EinMeV, EoutMeV);
			const double Ehi = max(EinMeV, EoutMeV);
			double conv_to_out_rate = (Eout>Ein ? exp((Ehi-Elo)/T) : 1.0);
			Phi[0] += nes_phi0_thompsonbruenn_(&Ehi, &Elo, &eta, &T, &s_nulib) * conv_to_out_rate;
			Phi[1] += nes_phi1_thompsonbruenn_(&Ehi, &Elo, &eta, &T, &s_nulib) * conv_to_out_rate;
		}
		else{
			assert(!__nulib_MOD_add_anue_iscattering_electrons);
			assert(!__nulib_MOD_add_numu_iscattering_electrons);
			assert(!__nulib_MOD_add_anumu_iscattering_electrons);
			assert(!__nulib_MOD_add_nutau_iscattering_electrons);
			assert(!__nulib_MOD_add_anutau_iscattering_electrons);
		}

		assert(Phi[0]>=0);
		assert(fabs(Phi[1])<=Phi[0]);
		for(int i=0; i<KMOMENTS; i++){
			assert(Phi[i]==Phi[i]);
			assert(fabs(Phi[i]) < std::numeric_limits<double>::infinity() );
		}
		return Phi;
	}


	//=============//
	// PAIR KERNEL //
	//=============//
	// only return annihilation kernels
	// cm^3/s
	array<double,KMOMENTS> Phi_pair(int s, double E /*erg*/, double Ebar /*erg*/){
	        const double TMeV = eos_variables[1];
		const double Terg = TMeV*1e6*cgs::units::eV;
		const double X    = E    / Terg;
		const double Xbar = Ebar / Terg;
		const int s_nulib = s+1; // Fortran indexing
		const int pro_or_ann = 2; // annihilation
		const double rho = eos_variables[0];
		const double ndens = rho / __nulib_MOD_m_ref; // total number density
		const double Xn = eos_variables[6];
		const double Xp = eos_variables[7];

		array<double,2> Phi;
		Phi[0] = 0;
		Phi[1] = 0;

		if(__nulib_MOD_add_nue_kernel_epannihil){
			assert(__nulib_MOD_add_anue_kernel_epannihil);
			assert(__nulib_MOD_add_numu_kernel_epannihil);
			assert(__nulib_MOD_add_anumu_kernel_epannihil);
			assert(__nulib_MOD_add_nutau_kernel_epannihil);
			assert(__nulib_MOD_add_anutau_kernel_epannihil);
			for(int mom=0; mom<KMOMENTS; mom++)
			  Phi[mom] += epannihil_phi_bruenn_(&X, &Xbar, &eta, &s_nulib, &pro_or_ann, &mom);
		}
		else{
			assert(!__nulib_MOD_add_anue_kernel_epannihil);
			assert(!__nulib_MOD_add_numu_kernel_epannihil);
			assert(!__nulib_MOD_add_anumu_kernel_epannihil);
			assert(!__nulib_MOD_add_nutau_kernel_epannihil);
			assert(!__nulib_MOD_add_anutau_kernel_epannihil);
		}

		if(__nulib_MOD_add_nue_kernel_bremsstrahlung){
			assert(__nulib_MOD_add_anue_kernel_bremsstrahlung);
		        assert(__nulib_MOD_add_numu_kernel_bremsstrahlung);
			assert(__nulib_MOD_add_anumu_kernel_bremsstrahlung);
			assert(__nulib_MOD_add_nutau_kernel_bremsstrahlung);
			assert(__nulib_MOD_add_anutau_kernel_bremsstrahlung);
			double n_N = -1.;

			// neutron-neutron
			n_N = ndens * Xn;
			assert(n_N>=0);
			assert(n_N<=ndens);
			Phi[0] += bremsstrahlung_phi0_hannestad_(&X, &Xbar, &TMeV, &n_N, &s_nulib, &pro_or_ann);

			// proton-proton
			n_N = ndens * Xp;
			assert(n_N>=0);
			assert(n_N<=ndens);
			Phi[0] += bremsstrahlung_phi0_hannestad_(&X, &Xbar, &TMeV, &n_N, &s_nulib, &pro_or_ann);

			// neutron-proton
			n_N = ndens * sqrt(Xn*Xp);
			assert(n_N>=0);
			assert(n_N<=ndens);
			Phi[0] += bremsstrahlung_phi0_hannestad_(&X, &Xbar, &TMeV, &n_N, &s_nulib, &pro_or_ann) * 28./3.;
		}
		else{
		        assert(!__nulib_MOD_add_anue_kernel_bremsstrahlung);
		        assert(!__nulib_MOD_add_numu_kernel_bremsstrahlung);
			assert(!__nulib_MOD_add_anumu_kernel_bremsstrahlung);
			assert(!__nulib_MOD_add_nutau_kernel_bremsstrahlung);
			assert(!__nulib_MOD_add_anutau_kernel_bremsstrahlung);
		}

		// convert to cm^3/s
		for(int mom=0; mom<KMOMENTS; mom++)
		    Phi[mom] *= 2.0*cgs::constants::GF*cgs::constants::GF * Terg*Terg * clight / (2.*M_PI) / pow(cgs::constants::hbarc,4);// from nulib.F90
		assert(std::abs(Phi[1]) <= Phi[0]);
		
		assert(Phi[0]>=0);
		assert(fabs(Phi[1])<=Phi[0]);
		for(int i=0; i<KMOMENTS; i++){
			assert(Phi[i]==Phi[i]);
			assert(fabs(Phi[i]) < std::numeric_limits<double>::infinity() );
		}
		return Phi;
	}

};
EAS eas;

#endif
