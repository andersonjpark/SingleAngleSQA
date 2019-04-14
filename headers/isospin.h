//Yonglin Zhu
//zhuygln@gmail.com
//Calculate isospin vector for neutrino and antineutrino in two flavor model.
//Ref[Phys.Rev. D94 (2016) no.10, 105006,Phys.Rev. D86 (2012) 085015]
#ifndef _H_ISOSPIN
#define _H_ISOSPIN

#include <cmath>

const complex<double> pauli[4][2][2] = {
  { {0, 1}, {1, 0} }, // x
  { {0,-I}, {I, 0} }, // y
  { {1, 0}, {0,-1} }, // z
  { {1, 0}, {0, 1} }  // t
};

array<double,4> pauli_decompose(const MATRIX<complex<double>,2,2>& M){
  array<double,4> coefficients;
  coefficients[0] = 0.5 * (real(M[1][0]) + real(M[0][1]));
  coefficients[1] = 0.5 * (imag(M[1][0]) - imag(M[0][1]));
  coefficients[2] = 0.5 * (real(M[0][0]) - real(M[1][1]));
  coefficients[3] = 0.5 * (real(M[0][0]) + real(M[1][1]));

  return coefficients;
}

template<typename T>
MATRIX<complex<double>,2,2> pauli_reconstruct(const array<T,4> coefficients){
  MATRIX<complex<double>,2,2> M;
  for(unsigned i=0; i<2; i++){
    for(unsigned j=0; j<2; j++){
      M[i][j] = 0;
      for(unsigned k=0; k<4; k++)
	M[i][j] += coefficients[k] * pauli[k][i][j];
    }
  }
  return M;
}

#endif
