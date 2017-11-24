//Yonglin Zhu
//zhuygln@gmail.com
//Calculate isospin vector for neutrino and antineutrino in two flavor model.
//Ref[Phys.Rev. D94 (2016) no.10, 105006,Phys.Rev. D86 (2012) 085015]
////
#include <cmath>
//
//
//
vector<double> isospin(MATRIX<complex<double> > Sf, MATRIX<complex<double> > Sfbar) 
{
//	int NE=E.size();
	vector<double> s(6);
//	vector<double> delta(NE);
	s[0] = real(Sf[e][e])*real(Sf[mu][e])+imag(Sf[e][e])*imag(Sf[mu][e]);
	s[1] = real(Sf[e][e])*imag(Sf[mu][e])-imag(Sf[e][e])*real(Sf[mu][e]);
	s[2] = (norm(Sf[e][e])-(pow(real(Sf[mu][e]),2)+pow(imag(Sf[e][mu]),2)))/2;
	s[3] =-(real(Sfbar[e][e])*real(Sfbar[mu][e])+imag(Sfbar[e][e])*imag(Sfbar[mu][e]));
        s[4] =  real(Sfbar[e][e])*imag(Sfbar[mu][e])-imag(Sfbar[e][e])*real(Sfbar[mu][e]);
        s[5] = -(norm(Sfbar[e][e])-(pow(real(Sfbar[mu][e]),2)+pow(imag(Sfbar[e][mu]),2)))/2;

	return s;
}

