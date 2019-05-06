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
	array<array<array<MATRIX<complex<double>,NF,NF>,NMOMENTS>,NE>,NM> MBackground = s.oscillated_moments(profile, s0);

#pragma omp parallel for collapse(2)
	for(int m=matter; m<=antimatter; m++){
		for(int i=0; i<NE; i++){

			// get nulib species indices
			const int se = (m==matter ? 0 : 1);
			const int sx = (m==matter ? 2 : 3);
			const state mbar = (m==matter ? antimatter : matter);

			// intermediate variables
			complex<double> unblock_in, unblock_out;
			MATRIX<complex<double>,NF,NF> block, Pi_plus, Pi_minus;

			// absorption and emission
			MATRIX<double,NF,NF> kappa_abs_avg = avg_matrix(eas.abs(se,s.E[i]), eas.abs(sx,s.E[i]));
			double E_kT = s.E[i]/(1e6*eV) / s.T;
			dfdr[m][i][e ][e ] += kappa_abs_avg[e ][e ] * eas.fermidirac(se, E_kT);
			dfdr[m][i][mu][mu] += kappa_abs_avg[mu][mu] * eas.fermidirac(sx, E_kT);
			for(flavour f1=e; f1<=mu; f1++)
				for(flavour f2=e; f2<=mu; f2++)
					dfdr[m][i][f1][f2] -= kappa_abs_avg[f1][f2] * s.fmatrixf[m][i][f1][f2];

			// scattering and pair annihilation
			double Vi = Vphase(i, s.Etop);
			for(int j=0; j<NE; j++){
				double Vj = Vphase(j, s.Etop);
				array<MATRIX<double,NF,NF>,KMOMENTS> Phi, PhiAvg,  PhiTilde;
				array<double,KMOMENTS> Phi_e, Phi_x;
				double conv_to_in_rate;

				//out-scattering from i to j and reverse
				Phi_e = eas.Phi_scat(se, s.T, s.E[i], s.E[j], Vi, Vj);
				Phi_x = eas.Phi_scat(sx, s.T, s.E[i], s.E[j], Vi, Vj);
				for(int mom=0; mom<KMOMENTS; mom++){
					PhiAvg[mom]   = avg_matrix(  Phi_e[mom], Phi_x[mom]);
					PhiTilde[mom] = tilde_matrix(Phi_e[mom], Phi_x[mom]);
					Phi[mom] = PhiAvg[mom] - PhiTilde[mom];
				}
				block = blocking_term(Phi, s.fmatrixf[m][i], MBackground[m][j]);
				// target energy difference negative of neutrino energy difference
				conv_to_in_rate = exp((s.E[j]-s.E[i])/(s.T*1e6*cgs::units::eV));
				for(flavour f1=e; f1<=mu; f1++)
					for(flavour f2=e; f2<=mu; f2++){

						complex<double> in_rate = conv_to_in_rate * (
								MBackground[m][j][0][f1][f2] * 0.5*Phi[0][f1][f2] +
								MBackground[m][j][1][f1][f2] * 1.5*Phi[1][f1][f2] -
								block[f1][f2]);

						complex<double> out_rate =
								s.fmatrixf[m][i][f1][f2] * 0.5*PhiAvg[0][f1][f2]*Vj -
								block[f1][f2];

						dfdr[m][i][f1][f2] += in_rate - out_rate;
					}

				// annihilation with group j and reverse
				Phi_e = eas.Phi_pair(se, s.T, s.E[i], s.E[j], Vi, Vj);
				Phi_e = eas.Phi_pair(sx, s.T, s.E[i], s.E[j], Vi, Vj);
				for(int mom=0; mom<KMOMENTS; mom++){
					PhiAvg[mom]   = avg_matrix(  Phi_e[mom], Phi_x[mom]);
					PhiTilde[mom] = tilde_matrix(Phi_e[mom], Phi_x[mom]);
					Phi[mom] = PhiAvg[mom] - PhiTilde[mom];
				}
				block = blocking_term(Phi, s.fmatrixf[m][i], MBackground[mbar][j]);
				conv_to_in_rate = exp(-(s.E[j]+s.E[i])/(s.T*1e6*cgs::units::eV));
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

						dfdr[m][i][f1][f2] += in_rate - out_rate;
				    }
				}

			} // other group

			// Make sure dfdr is Hermitian
			dfdr[m][i][mu][e ] = conj(dfdr[m][i][e][mu]);
		} // group
	} // state

	return dfdr;
}



#endif
