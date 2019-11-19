#ifndef STATE_H
#define STATE_H

#include <algorithm>
#include "misc.h"
#include "mixing_angles.h"
#include "isospin.h"
#include "profile.h"
#include "nulib_interface.h"

class State{
public:
	double Ecom_Elab, Elab_Elabstart;
	double r;
	double rho, T, Ye; // g/ccm, MeV

	// energy grid
	array<double,NE> E, Etop, Ecom, Etopcom; // erg

	// distribution function in the direction of the trajectory
	// value at the last reset
	array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> fmatrixf;

	// Evolution variables for neutrino oscillation
	// Y Describes oscillation since Scumulative was last updated
	// Scumulative is evolution matrix from initial mass basis (s0.UU) to mass basis
	// at the time of the last update. S = WBWB(Y) * Scumulative to current mass basis.
	array<array<array<array<double,NY>,NS>,NE>,NM> Y;
	array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> Scumulative;

	// Intermediate quantities used in calculating potentials and K
	array<array<array<double,NF>,NE>,NM> kk;
	array<array<array<double,NF-1>,NE>,NM> dkk;
	array<MATRIX<complex<double>,NF,NF>,NM> VfSI;
	array<array<array<array<double,NF>,NF>,NE>,NM> AA;
	array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> Sf, SThisStep, VfMSW, dVfMSWdr, UU;
	array<array<array<MATRIX<complex<double>,NF,NF>,NF>,NE>,NM> CC;
	array<array<array<MATRIX<complex<double>,NF,NF>,NS>,NE>,NM> BB,WW;

	// other matrices
	array<array<double,NM>,NE> dphi_dr_interact, dtheta_dr_interact;
	array<array<double,NM>,NE> dphi_dr_osc,      dtheta_dr_osc;

	State(const Profile& profile, const double rmin){
		rho=NAN;
		T=NAN;
		Ye=NAN;
		r=rmin;
		Ecom_Elab=profile.Ecom_Elab(r);
		Elab_Elabstart=NAN;

		// set energy bins to match the profile at rmin in the comoving frame
		cout << endl;
		cout<<"NE="<<NE << "   Ecom0/Elab0=" << Ecom_Elab << endl;
		cout << "mid(com) \t mid(lab) \t top(lab)" << endl;
		this->Ecom = profile.Ecom;
		this->Etopcom = profile.Etopcom;
		for(int i=0;i<NE;i++){
			E[i]    = Ecom[i]    / Ecom_Elab;
			Etop[i] = Etopcom[i] / Ecom_Elab;
			cout << Ecom[i] / (1.e6*cgs::units::eV) << " \t ";
			cout << E[i]    / (1.e6*cgs::units::eV) << " \t ";
			cout << Etop[i] / (1.e6*cgs::units::eV) << endl;
		}
		cout.flush();


		// set Scumulative to identity
		for(int m=0; m<NM; m++){
			for(int ig=0; ig<NE; ig++){
				for(int f1=0; f1<NF; f1++){
					for(int f2=0; f2<NF; f2++){
						Scumulative[m][ig][f1][f2] = 0.;
						Sf[m][ig][f1][f2] = 0.;
					}
					Scumulative[m][ig][f1][f1] = 1.;
					Sf[m][ig][f1][f1] = 1.;
				}
				Y[m][ig] = YIdentity;
			}
		}
	}

	int find_bin(const double x) const{
	        if(x<0) return -1;

		// upper_bound returns first element greater than xval
		// values mark bin tops, so this is what we want
		int ind = std::upper_bound(Etop.begin(), Etop.end(), x) - Etop.begin();
		assert(ind >= 0);
		assert(ind <= NE);
		return ind;
	}

	void update_background(const Profile& profile, const State& s0){

		rho = exp(profile.lnrho(r));
		T = profile.temperature(r);
		Ye = profile.Ye(r);
		Ecom_Elab = profile.Ecom_Elab(r);
		Elab_Elabstart = profile.Elab_Elabstart(r);

		for(int i=0; i<NE; i++){
			E[i]       = s0.E[i]    * Elab_Elabstart;
			Etop[i]    = s0.Etop[i] * Elab_Elabstart;
			Ecom[i]    = E[i]       * Ecom_Elab;
			Etopcom[i] = Etop[i]    * Ecom_Elab;
		}
	}

	void update_potential(const Profile& profile, const State& s0){
		// vacuum potential
		array<array<double,NF>,NE> kV = set_kV(E);
		array<MATRIX<complex<double>,NF,NF>,NM> UV = Evaluate_UV();
		array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> VfVac = Evaluate_VfVac(kV,UV);

		// derivative of vacuum potential (which is proportional to 1/E)
		double VfVac_derivative_fac = -profile.Elab_Elabstart.Derivative(r) / Elab_Elabstart;

		// Matter Potential
		array<MATRIX<complex<double>,NF,NF>,NM> VfMatter, dVfMatterdr;
		double matter_potential=M_SQRT2*cgs::constants::GF/cgs::constants::Mp*rho*Ye*Ecom_Elab;
		VfMatter[matter][e ][e ] = matter_potential;
		VfMatter[matter][mu][mu] = 0;
		VfMatter[matter][e ][mu] = 0;
		VfMatter[matter][mu][e ] = 0;
		VfMatter[antimatter]=-Conjugate(VfMatter[matter]);

		// derivative of matter potential
		double dlogrhodr=profile.lnrho.Derivative(r);
		double dYedr=profile.Ye.Derivative(r);
		dVfMatterdr[matter] = VfMatter[matter] * (dlogrhodr + dYedr/Ye);
		dVfMatterdr[antimatter]=-Conjugate(dVfMatterdr[matter]);

		// SI potential
		VfSI[matter] = MATRIX<complex<double>,NF,NF>();
		VfSI[antimatter] = MATRIX<complex<double>,NF,NF>();
#pragma omp parallel for collapse(2)
		for(int m=matter; m<=antimatter; m++){
			for(int i=0;i<=NE-1;i++){

				// stuff that used to be in K()
				MATRIX<complex<double>,NF,NF> dVfVacdr = VfVac[m][i] * VfVac_derivative_fac;
				VfMSW[m][i] = VfVac[m][i]+VfMatter[m];
				dVfMSWdr[m][i] = dVfMatterdr[m] + dVfVacdr;
				kk[m][i] = k(VfMSW[m][i]);
				dkk[m][i] = deltak(VfMSW[m][i]);
				CC[m][i]  = CofactorMatrices(VfMSW[m][i],kk[m][i]);
				AA[m][i] = MixingMatrixFactors(CC[m][i],s0.CC[m][i],s0.AA[m][i]);
				UU[m][i] = U(dkk[m][i],CC[m][i],AA[m][i]);
				BB[m][i][msw] = B(Y[m][i][msw]);
				BB[m][i][si ] = B(Y[m][i][si ]);
				WW[m][i][msw] = W(Y[m][i][msw]);
				WW[m][i][si ] = W(Y[m][i][si ]);
				SThisStep[m][i] = WW[m][i][msw] * BB[m][i][msw] * WW[m][i][si] * BB[m][i][si];
				Sf[m][i] = UU[m][i] * SThisStep[m][i] * Scumulative[m][i] * Adjoint(s0.UU[m][i]);

				assert( abs( Trace( Scumulative[m][i]*Adjoint(Scumulative[m][i])) - (double)NF) < 1e-5);
				assert( abs( Trace(          Sf[m][i]*Adjoint(         Sf[m][i])) - (double)NF) < 1e-5);
			}
		}

		array<array<array<MATRIX<complex<double>,NF,NF>,NMOMENTS>,NE>,NM> MBackground = oscillated_moments(profile, s0);

		// calculate the self-interaction potential
		for(int m=matter; m<=antimatter; m++){
			for(int i0=0; i0<NE; i0++){
				MATRIX<complex<double>,NF,NF> VfSIE = (MBackground[m][i0][0] - MBackground[m][i0][1]) * sqrt(2.)*cgs::constants::GF;
				VfSI[matter] += (m==matter ? VfSIE : -Conjugate(VfSIE));
			}
		}

		// convert comoving-frame potential to lab-frame potential
		VfSI[matter] *= Ecom_Elab;

		// set antimatter potential
		VfSI[antimatter]=-Conjugate(VfSI[matter]);
	}


	//====================//
	// OSCILLATED MOMENTS //
	//====================//
	// oscillate moments, keeping the moments' comoving-frame energy grid
	array<array<array<MATRIX<complex<double>,NF,NF>,NMOMENTS>,NE>,NM> oscillated_moments(const Profile& profile, const State& s0) const{

		array<array<array<MATRIX<complex<double>,NF,NF>,NMOMENTS>,NE>,NM> MBackground;

#pragma omp parallel for collapse(2)
		for(int m=matter; m<=antimatter; m++){
			for(int i0=0; i0<NE; i0++){
				double V0 = Vphase(i0, s0.Etop);

				// fill in the un-oscillated diagonals
				array<MATRIX<complex<double>,NF,NF>,NMOMENTS> unosc_moment;
				for(flavour f=e; f<=mu; f++){
					unosc_moment[0][f][f] += profile.Dens_unosc[m][i0][f](r);
					unosc_moment[1][f][f] += profile.Flux_unosc[m][i0][f](r);
					unosc_moment[2][f][f] += profile.Pres_unosc[m][i0][f](r);
				}

				double total_overlap_fraction = 0; // should end up being 1
				for(int ilab=0; ilab<NE; ilab++){
					assert( abs(Trace(Sf[m][ilab]*Adjoint(Sf[m][ilab])) - (double)NF) < 1e-5);

					// calculate fraction of bin i0 that overlaps with bin ilab
					double V_overlap = Vphase_overlap_comoving(i0, s0.Etop, ilab, Etop, Ecom_Elab);

					// calculate contribution to the moments due to this overlapping
					// segment of bin i0, oscillating it with Sf from bin ilab
					// oscillate the moments
					if(V_overlap>0){
					  double overlap_fraction = V_overlap / V0;
					  total_overlap_fraction += overlap_fraction;

					  if(ASSUME_ISOTROPY){
					    MBackground[m][i0][0] += fmatrixf[m][ilab] * V_overlap;
					    MBackground[m][i0][1] += MBackground[m][i0][0] * 0;
					    MBackground[m][i0][2] += MBackground[m][i0][0] * 1./3.;
					  }
					  else for(int mom=0; mom<NMOMENTS; mom++)
						 MBackground[m][i0][mom] += Sf[m][ilab]*unosc_moment[mom]*Adjoint(Sf[m][ilab]) * overlap_fraction;
					}
				}

				// add in rest of moment using highest-bin Sf
				assert(total_overlap_fraction < 1.+1e-6);
				assert(total_overlap_fraction >= 0);
				if(total_overlap_fraction < 1.-1e-6){
					assert(s0.Etop[i0] > Etop[NE-1]*Ecom_Elab);
					for(int mom=0; mom<NMOMENTS; mom++)
						MBackground[m][i0][mom] += Sf[m][NE-1]*unosc_moment[mom]*Adjoint(Sf[m][NE-1]) * (1.-total_overlap_fraction);
				}

				// check that it's reasonable
				for(flavour f1=e; f1<=mu; f1++){
				  for(flavour f2=e; f2<=mu; f2++){
				    assert(abs(MBackground[m][i0][0][f1][f1])/V0 <= 1.);
				  }
				  assert(real(MBackground[m][i0][0][f1][f1])/V0 >= 0.);
				  assert(real(MBackground[m][i0][0][f1][f1])/V0 <= 1.);
				}
			}
		}

		return MBackground;
	}

	void accumulate_S(double dr, const State& sReset, const State& s0){
#pragma omp parallel for collapse(2)
		for(int m=matter;m<=antimatter;m++){
			for(int i=0;i<=NE-1;i++){
				Scumulative[m][i] = SThisStep[m][i] * Scumulative[m][i];
				Sf[m][i] = UU[m][i] * Scumulative[m][i] * Adjoint(s0.UU[m][i]);
				assert( abs( Trace( Scumulative[m][i]*Adjoint(Scumulative[m][i])) - (double)NF) < 1e-5);
				assert( abs( Trace(          Sf[m][i]*Adjoint(         Sf[m][i])) - (double)NF) < 1e-5);

				// convert fmatrix from flavor basis to (reset-point) mass basis
				// evolve fmatrix from reset-point to current-point mass basis
				// convert fmatrix from (current-point) mass basis to flavor basis
				MATRIX<complex<double>,NF,NF> SfThisStep =
						UU[m][i]
							  * SThisStep[m][i]
											 * Adjoint(sReset.UU[m][i]);
				fmatrixf[m][i] = SfThisStep * fmatrixf[m][i] * Adjoint(SfThisStep);

				// reset the evolution matrix to identity
				Y[m][i] = YIdentity;

				// get rate of change of fmatrix from oscillation
				array<double,4> hold = pauli_decompose(sReset.fmatrixf[m][i]);
				array<double,4> hnew = pauli_decompose(       fmatrixf[m][i]);
				double oldmag   = sqrt(hold[0]*hold[0] + hold[1]*hold[1] + hold[2]*hold[2]);
				double newmag   = sqrt(hnew[0]*hnew[0] + hnew[1]*hnew[1] + hnew[2]*hnew[2]);
				double costheta = (hold[0]*hnew[0] + hold[1]*hnew[1] + hold[2]*hnew[2]) / (newmag*oldmag);
				assert(costheta-1. < 1e-10);
				costheta = min(1.,costheta);
				dtheta_dr_osc[i][m] = (acos(hnew[2]/newmag) - acos(hold[2]/oldmag)) / dr;
				dphi_dr_osc[i][m] = (atan2(hnew[1],hnew[0]) - atan2(hold[1],hold[0])) / dr;
			}
		}
	}

	void assert_noNaN(double accuracy){
		for(state m=matter; m<=antimatter; m++){
			for(int i=0; i<NE; i++){

				for(flavour f1=e; f1<=mu; f1++){
					assert(real(fmatrixf[m][i][f1][f1]) <= 1.);
					assert(real(fmatrixf[m][i][f1][f1]) >= 0.);
					assert(abs(imag(fmatrixf[m][i][f1][f1])) < accuracy);
					for(flavour f2=e; f2<=mu; f2++)
						assert(fmatrixf[m][i][f1][f2] == fmatrixf[m][i][f1][f2]);
				}

				array<double,4> isospin = pauli_decompose(fmatrixf[m][i]);
				double fperp2 = isospin[0]*isospin[0] + isospin[1]*isospin[1];
				assert(fabs(isospin[0]) <= isospin[3]);
				assert(fperp2*(1.-1e10) <= pow(min(isospin[3], 1.-isospin[3]),2) - isospin[2]*isospin[2]);

				for(solution x=msw;x<=si;x++)
					for(int j=0;j<=NY-1;j++)
						assert(Y[m][i][x][j] == Y[m][i][x][j]);
			}
		}
	}


	//============//
	// Initialize //
	//============//
	void initialize(const array<array<array<DISCONTINUOUS,NF>,NE>,NM>& D_unosc){
		array<array<double,NF>,NE> kV = set_kV(E);
		array<MATRIX<complex<double>,NF,NF>,NM> UV = Evaluate_UV();
		array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> VfVac = Evaluate_VfVac(kV,UV);
		array<array<MATRIX<complex<double>,NF,NF>,NF>,NE> CV = Evaluate_CV(kV, VfVac);
		array<array<array<double,NF>,NF>,NE> AV = Evaluate_AV(kV,VfVac,UV);
		for(state m=matter; m<=antimatter; m++){
			for(int i=0;i<=NE-1;i++){
				for(int j=0;j<=NF-1;j++){
					if(real(CC[m][i][j][mu][e]*CV[i][j][mu][e]) < 0.)
						AA[m][i][j][e]=-AV[i][j][e];
					else AA[m][i][j][e]=AV[i][j][e];
					AA[m][i][j][mu]=AV[i][j][mu];
				}
				UU[m][i]=U(dkk[m][i],CC[m][i],AA[m][i]);
			}
		}

		// determine eigenvalue ordering
		if(kV[0][1]>kV[0][0])
			cout<<"\n\nNormal hierarchy" << endl;
		else{
			if(kV[0][1]<kV[0][0])
				cout<<"\n\nInverted hierarchy" << endl;
			else{
				cout<<endl<<endl<<"Neither normal or Inverted"<<endl;
				abort();
			}
		}


		// T should be MeV
		cout << "Setting initial data." << endl;
		cout << "rho = " << rho << " g/ccm" << endl;
		cout << "T = " << T << " MeV" << endl;
		cout << "Ye = " << Ye << endl;
		//Ye = max(Ye,__nulibtable_MOD_nulibtable_ye_min);
		/* nulibtable_range_species_range_energy_(&rho, &T, &Ye, &eas.eas.front(), */
		/* 					 &__nulibtable_MOD_nulibtable_number_species, */
		/* 					 &__nulibtable_MOD_nulibtable_number_groups, */
		/* 					 &__nulibtable_MOD_nulibtable_number_easvariables); */

		for(int i=0; i<NE; i++){
			for(state m=matter; m<=antimatter; m++){
				fmatrixf[m][i] = MATRIX<complex<double>,NF,NF>();
				for(flavour f=e; f<=mu; f++){
				  //double D = D_unosc[m][i][f](r);
				  //double V = Vphase(i,Etop);
				  fmatrixf[m][i][f][f] = 1.e-200; //D / V;
				  assert(abs(fmatrixf[m][i][f][f]) >= 0);
				  assert(abs(fmatrixf[m][i][f][f]) <= 1);
				}
			}

			cout << "GROUP " << i << " f = {";
			cout << real(fmatrixf[matter][i][e][e]) << ", ";
			cout << real(fmatrixf[antimatter][i][e][e]) << ", ";
			cout << real(fmatrixf[matter][i][mu][mu]) <<"}" << endl;
		}
	}


};



#endif
