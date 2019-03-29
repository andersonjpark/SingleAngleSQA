//Yonglin Zhu
//zhuygln@gmail.com
//Calculate isospin vector for neutrino and antineutrino in two flavor model.
//Ref[Phys.Rev. D94 (2016) no.10, 105006,Phys.Rev. D86 (2012) 085015]
#ifndef _H_ISOSPIN
#define _H_ISOSPIN

#include <cmath>

vector<double> isospin(MATRIX<complex<double>,NF,NF> Sf, MATRIX<complex<double>,NF,NF> Sfbar){
  //int NE=E.size();
  vector<double> s(6);
  //vector<double> delta(NE);
  s[0] = real(Sf[e][e])*real(Sf[mu][e])+imag(Sf[e][e])*imag(Sf[mu][e]);
  s[1] = real(Sf[e][e])*imag(Sf[mu][e])-imag(Sf[e][e])*real(Sf[mu][e]);
  s[2] = (norm(Sf[e][e])-(pow(real(Sf[mu][e]),2)+pow(imag(Sf[e][mu]),2)))/2;
  s[3] =-(real(Sfbar[e][e])*real(Sfbar[mu][e])+imag(Sfbar[e][e])*imag(Sfbar[mu][e]));
  s[4] =  real(Sfbar[e][e])*imag(Sfbar[mu][e])-imag(Sfbar[e][e])*real(Sfbar[mu][e]);
  s[5] = -(norm(Sfbar[e][e])-(pow(real(Sfbar[mu][e]),2)+pow(imag(Sfbar[e][mu]),2)))/2;
  
  return s;
}

complex<double> pauli[4][2][2] = {
  { {0, 1}, {1, 0} }, // x
  { {0,-I}, {I, 0} }, // y
  { {1, 0}, {0,-1} }, // z
  { {1, 0}, {0, 1} }  // t
};

void pauli_decompose(const MATRIX<complex<double>,2,2>& M, double coefficients[4]){
  //assert(M[1][0] == conj(M[0][1]));
  //assert(imag(M[0][0]) == 0);
  //assert(imag(M[1][1]) == 0);

  coefficients[0] = 0.5 * (real(M[1][0]) + real(M[0][1]));
  coefficients[1] = 0.5 * (imag(M[1][0]) - imag(M[0][1]));
  coefficients[2] = 0.5 * (real(M[0][0]) - real(M[1][1]));
  coefficients[3] = 0.5 * (real(M[0][0]) + real(M[1][1]));
}

template<typename T>
void pauli_reconstruct(const T coefficients[4], MATRIX<complex<double>,2,2>& M){
  for(unsigned i=0; i<2; i++){
    for(unsigned j=0; j<2; j++){
      M[i][j] = 0;
      for(unsigned k=0; k<4; k++)
	M[i][j] += coefficients[k] * pauli[k][i][j];
    }
  }
}

#endif
