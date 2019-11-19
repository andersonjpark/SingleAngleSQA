#ifndef TIMEDERIVATIVES_H
#define TIMEDERIVATIVES_H

#include "nulib_interface.h"

//============//
// Koscillate //
//============//
array<array<array<array<double,NY>,NS>,NE>,NM> Koscillate(const State& s){

	array<array<array<array<double,NY>,NS>,NE>,NM> K;

#pragma omp parallel for collapse(2)
	for(int m=matter; m<=antimatter; m++){
		for(int i=0;i<=NE-1;i++){

			array<double,NF-1> phase;
			MATRIX<complex<double>,NF,NF> Ha,HB;
			phase[0] = M_2PI*(s.Y[m][i][msw][4]-s.Y[m][i][msw][5]);
			Ha[0][1]=0.;
			for(int j=0;j<=NF-2;j++)
				for(int k=j+1;k<=NF-1;k++)
					for(flavour f=e;f<=mu;f++)
						Ha[j][k]+= conj(s.UU[m][i][f][j])*s.dVfMSWdr[m][i][f][f]*s.UU[m][i][f][k];

			Ha[0][1] *= I*cgs::constants::hbarc/s.dkk[m][i][0]*exp(I*phase[0]);
			Ha[1][0] = conj(Ha[0][1]);

			// HB = -I/cgs::constants::hbarc*Ha*BB;
			HB[0][0]=-I/cgs::constants::hbarc*( Ha[0][1]*s.BB[m][i][msw][1][0] );
			HB[0][1]=-I/cgs::constants::hbarc*( Ha[0][1]*s.BB[m][i][msw][1][1] );

			array<double,4> dvdr;
			dvdr[0]=real(HB[0][1]);
			dvdr[1]=imag(HB[0][1]);
			dvdr[2]=real(HB[0][0]);
			dvdr[3]=imag(HB[0][0]);

			MATRIX<double,3,4> JI = JInverse(s.Y[m][i][msw]);

			array<double,NF> dkkdr = dkdr(s.UU[m][i],s.dVfMSWdr[m][i]);
			array<MATRIX<complex<double>,NF,NF>,NF> dCCdr = CofactorMatricesDerivatives(s.dVfMSWdr[m][i],dkkdr);
			array<double,NF> QQ =  Q(s.UU[m][i],s.dkk[m][i],s.CC[m][i],dCCdr);

			for(int j=0;j<=2;j++){
				K[m][i][msw][j]=0.;
				for(int k=j;k<=3;k++)
					K[m][i][msw][j] += JI[j][k]*dvdr[k];
			}
			K[m][i][msw][3] = 0.;
			K[m][i][msw][4] = (s.kk[m][i][0]+QQ[0])/M_2PI/cgs::constants::hbarc;
			K[m][i][msw][5] = (s.kk[m][i][1]+QQ[1])/M_2PI/cgs::constants::hbarc;


			// *********************
			// SI part of solution *
			// *********************
			MATRIX<complex<double>,NF,NF> UWBW =
					s.UU[m][i]
							* s.WW[m][i][msw]
										 * s.BB[m][i][msw]
													  * s.WW[m][i][si];
			Ha = Adjoint(UWBW)*s.VfSI[m]*UWBW;

			K[m][i][si][4]=real(Ha[0][0])/(M_2PI*cgs::constants::hbarc);
			K[m][i][si][5]=real(Ha[1][1])/(M_2PI*cgs::constants::hbarc);

			HB[0][0]=-I/cgs::constants::hbarc*( Ha[0][1]*s.BB[m][i][si][1][0] );
			HB[0][1]=-I/cgs::constants::hbarc*( Ha[0][1]*s.BB[m][i][si][1][1] );

			//array<double,4> dvdr;
			dvdr[0]=real(HB[0][1]);
			dvdr[1]=imag(HB[0][1]);
			dvdr[2]=real(HB[0][0]);
			dvdr[3]=imag(HB[0][0]);

			//MATRIX<double,3,4>
			JI = JInverse(s.Y[m][i][si]);

			for(int j=0;j<=2;j++){
				K[m][i][si][j]=0.;
				for(int k=j;k<=3;k++) K[m][i][si][j]+=JI[j][k]*dvdr[k];
			}

			K[m][i][si][3]=0.;

			for(solution x=msw;x<=si;x++)
				for(int j=0;j<=NY-1;j++)
					assert(K[m][i][x][j] == K[m][i][x][j]);
		}
	}

	return K;
}// end of K function


//==================//
// HELPER FUNCTIONS //
//==================//
MATRIX<double,2,2> avg_matrix(const double eval, const double muval){
	MATRIX<double,2,2> result;
	result[e ][e ] = eval;
	result[mu][mu] = muval;
	result[e ][mu] = result[mu][e] = (eval + muval) / 2.;
	return result;
}
MATRIX<double,2,2> tilde_matrix(const double eval, const double muval){
	MATRIX<double,2,2> result;
	result[e ][e ] = 0;
	result[mu][mu] = 0;
	result[e ][mu] = result[mu][e] = (eval - muval) / (4.*sin2thetaW);
	return result;
}

MATRIX<complex<double>,NF,NF>
blocking_term(const array<MATRIX<double,NF,NF>,KMOMENTS>& Phi /*cm^3/s/sr*/,
		const MATRIX<complex<double>,NF,NF>& f/*dimensionless*/,
		const array<MATRIX<complex<double>,NF,NF>,NMOMENTS>& Mp /* cm^-3*/){

	MATRIX<complex<double>,NF,NF> result;
	for(flavour fa=e; fa<=mu; fa++)
		for(flavour fb=e; fb<=mu; fb++){
			result[fa][fb] = 0;

			for(flavour fc=e; fc<=mu; fc++){
				result[fa][fb] += 0.5 * Phi[0][fc][fb]*f[fa][fc]*Mp[0][fc][fb];
				result[fa][fb] += 0.5 * Phi[0][fa][fc]*Mp[0][fa][fc]*f[fc][fb];
				result[fa][fb] += 1.5 * Phi[1][fc][fb]*f[fa][fc]*Mp[1][fc][fb];
				result[fa][fb] += 1.5 * Phi[1][fa][fc]*Mp[1][fa][fc]*f[fc][fb];
			}
		}

	return result * 0.5; // 1/s (really M involves vol int. which gets rid of sr units)
}


//===========//
// Kinteract //
//===========//
array<array<MATRIX<complex<double>,NF,NF>,NE>,NM>
Kinteract(const State& s, const State& s0, const Profile& profile){

	// set up the array to be returned
	array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> dfdr;

	// get oscillated background density
	array<array<array<MATRIX<complex<double>,NF,NF>,NMOMENTS>,NE>,NM> MBackground = s.oscillated_moments(profile, s0.Etopcom);

#pragma omp parallel for collapse(2)
	for(int m=matter; m<=antimatter; m++){
		for(int i=0; i<NE; i++){

			// get nulib species indices
			const int se = (m==matter ? 0 : 1);
			const int sx = (m==matter ? 2 : 3);
			const state mbar = (m==matter ? antimatter : matter);

			// intermediate variables
			MATRIX<complex<double>,NF,NF> block, Pi_plus, Pi_minus, tmp;
			array<MATRIX<double,NF,NF>,KMOMENTS> Phi, PhiAvg,  PhiTilde;
			array<double,KMOMENTS> Phi_e, Phi_x;
			double conv_to_in_rate;

			//=========================//
			// absorption and emission //
			//=========================//
			MATRIX<double,NF,NF> kappa_abs_avg = avg_matrix(eas.abs(se,s.Ecom[i]), eas.abs(sx,s.Ecom[i]));
			double E_kT = s.Ecom[i]/(1e6*eV) / s.T;
			dfdr[m][i][e ][e ] += kappa_abs_avg[e ][e ] * eas.fermidirac(se, E_kT);
			dfdr[m][i][mu][mu] += kappa_abs_avg[mu][mu] * eas.fermidirac(sx, E_kT);
			for(flavour f1=e; f1<=mu; f1++)
				for(flavour f2=e; f2<=mu; f2++)
					dfdr[m][i][f1][f2] -= kappa_abs_avg[f1][f2] * s.fmatrixf[m][i][f1][f2];


			// scattering and pair annihilation
			for(int j=0; j<NE; j++){
				double Vj = Vphase(j, s0.Etop); // cm^-3

				//============//
				// scattering //
				//============//
				bool include_elastic = (j == s0.find_bin(s.Ecom[i])); // E[i] is inside current background energy bin j
				Phi_e = eas.Phi_scat(se, s.T, s.Ecom[i], s0.Ecom[j], Vj, include_elastic); // cm^3/s
				Phi_x = eas.Phi_scat(sx, s.T, s.Ecom[i], s0.Ecom[j], Vj, include_elastic); // cm^3/s
				for(int mom=0; mom<KMOMENTS; mom++){
					PhiAvg[mom]   = avg_matrix(  Phi_e[mom], Phi_x[mom]);
					PhiTilde[mom] = tilde_matrix(Phi_e[mom], Phi_x[mom]);
					Phi[mom] = PhiAvg[mom] - PhiTilde[mom];
				}
				block = blocking_term(Phi, s.fmatrixf[m][i], MBackground[m][j]);
				// target energy difference negative of neutrino energy difference
				conv_to_in_rate = exp((s0.Ecom[j]-s.Ecom[i])/(s.T*1e6*cgs::units::eV));
				for(flavour f1=e; f1<=mu; f1++){
					for(flavour f2=e; f2<=mu; f2++){

						complex<double> in_rate = conv_to_in_rate * (
								MBackground[m][j][0][f1][f2] * 0.5*Phi[0][f1][f2] +
								MBackground[m][j][1][f1][f2] * 1.5*Phi[1][f1][f2] -
								block[f1][f2]);

						complex<double> out_rate =
								s.fmatrixf[m][i][f1][f2] * 0.5*PhiAvg[0][f1][f2]*Vj -
								block[f1][f2];

						dfdr[m][i][f1][f2] += (in_rate - out_rate) / (cgs::constants::c);
					}
				}

				//==============//
				// annihilation //
				//==============//
				Phi_e = eas.Phi_pair(se, s.Ecom[i], s0.Ecom[j]); // cm^3/s
				Phi_x = eas.Phi_pair(sx, s.Ecom[i], s0.Ecom[j]); // cm^3/s
				for(int mom=0; mom<KMOMENTS; mom++){
					PhiAvg[mom]   = avg_matrix(  Phi_e[mom], Phi_x[mom]);
					PhiTilde[mom] = tilde_matrix(Phi_e[mom], Phi_x[mom]);
					Phi[mom] = PhiAvg[mom] - PhiTilde[mom];
				}
				block = blocking_term(Phi, s.fmatrixf[m][i], MBackground[mbar][j]);
				conv_to_in_rate = exp(-(s0.Ecom[j]+s.Ecom[i])/(s.T*1e6*cgs::units::eV));
				for(flavour f1=e; f1<=mu; f1++){
					for(flavour f2=e; f2<=mu; f2++){
						complex<double> in_rate = conv_to_in_rate * (
								(f1==f2 ? 1. : 0.)         * Vj   * 0.5*Phi   [0][f1][f2]
								- s.fmatrixf[m][i][f1][f2] * Vj   * 0.5*PhiAvg[0][f1][f2]
								- MBackground[mbar][j][0][f1][f2] * 0.5*Phi   [0][f1][f2]
								- MBackground[mbar][j][1][f1][f2] * 1.5*Phi   [1][f1][f2]
								+ block[f1][f2]
								);

						complex<double> out_rate = block[f1][f2];

						dfdr[m][i][f1][f2] += (in_rate - out_rate) / (cgs::constants::c);
				    }
				}


				//===============//
				// nu4 processes //
				//===============//
				if(__nulib_MOD_add_nu4scat_kernel or __nulib_MOD_add_nu4pair_kernel){
					assert(__nulib_MOD_add_nu4scat_kernel);
					assert(__nulib_MOD_add_nu4pair_kernel);
					double k = s.Ecom[i];
					double kernel = NAN;

					double q1 = s0.Ecom[j];
					for(int j3=0; j3<NE; j3++){
						double q3 = s0.Ecom[j3];
						double V3 = Vphase(j3, s0.Etop); // cm^-3

						double q2 = q1+q3-k;
						int j2 = s0.find_bin(q2);

						if(j2>=0 and j2<NE){
							double V2 = Vphase(j2, s0.Etop); // cm^-3
							MATRIX<complex<double>,NF,NF> fj, fj2, fj3;

							// SCATTERING
							fj  = MBackground[m][j ][0] / Vj;
							fj2 = MBackground[m][j2][0] / V2;
							fj3 = MBackground[m][j3][0] / V3;
							kernel = __nulib_MOD_nu4scat_kernel_single(&k, &q1, &q2, &q3) * Vj*V3; // 1/cm

							tmp = (1.-fj) * fj2;
							Pi_minus += (Trace(tmp) + tmp) * (1.-fj3) * kernel;

							tmp = fj*(1.-fj2);
							Pi_plus  += (Trace(tmp) + tmp) * fj3 * kernel;

							// PAIR PROCESSES
							fj  = MBackground[mbar][j ][0] / Vj;
							fj2 = MBackground[mbar][j2][0] / V2;
							kernel = __nulib_MOD_nu4pair_kernel_single(&k, &q1, &q2, &q3) * Vj*V3; // 1/cm

							tmp = fj2 * (1.-fj);
							Pi_minus += (Trace(tmp) + tmp) * (1.-fj3) * kernel;

							tmp = (1.-fj3) * (1.-fj);
							Pi_minus += (Trace(tmp) + tmp) * fj2 * kernel;

							tmp = (1.-fj2)*fj;
							Pi_plus += (Trace(tmp) + tmp) * fj3 * kernel;

							tmp = fj3*fj;
							Pi_plus += (Trace(tmp) + tmp) * (1.-fj2) * kernel;
						}
					}
					dfdr[m][i] += Pi_plus *(1.-s.fmatrixf[m][i]) + (1.-s.fmatrixf[m][i])*Pi_plus ;
					dfdr[m][i] -= Pi_minus*    s.fmatrixf[m][i]  +     s.fmatrixf[m][i] *Pi_minus;
				} // nu4 processes


			} // other group

			// Make sure dfdr is Hermitian
			dfdr[m][i][mu][e ] = conj(dfdr[m][i][e][mu]);
			assert(dfdr[m][i][e ][e ] == dfdr[m][i][e ][e ]);
			assert(dfdr[m][i][e ][mu] == dfdr[m][i][e ][mu]);
			assert(dfdr[m][i][mu][mu] == dfdr[m][i][mu][mu]);
		} // group
	} // state

	return dfdr;
}



#endif
