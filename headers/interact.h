#ifndef INTERACT_H
#define INTERACT_H

void interact(double dr, State& s,
	      const array<array<array<DISCONTINUOUS,NF>,NE>,NM>& D_unosc){
  
  // save old fmatrix
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> old_fmatrixf = s.fmatrixf;

  // let neutrinos interact
  my_interact(s.fmatrixf, dr, s, D_unosc);
  
  // loop through getting rotation matrices
  double hold[4];
  double hnew[4];
  for(int i=0; i<NE; i++){
    for(state m=matter; m<=antimatter; m++){
      pauli_decompose(old_fmatrixf[m][i], hold);
      pauli_decompose(  s.fmatrixf[m][i], hnew);

      // get the theta and phi contribution
      double oldmag2   = hold[0]*hold[0] + hold[1]*hold[1] + hold[2]*hold[2];
      double newmag2   = hnew[0]*hnew[0] + hnew[1]*hnew[1] + hnew[2]*hnew[2];
      if(oldmag2==0 or newmag2==0){
	continue;
      }
      s.dtheta_dr_interact[i][m] = (acos(hnew[2]/sqrt(newmag2)) - acos(hold[2]/sqrt(oldmag2))) / dr;
      s.dphi_dr_interact[i][m] = (atan2(hnew[1],hnew[0]) - atan2(hold[1],hold[0])) / dr;

      // get the axis of rotation
      double lrot[3];
      lrot[0] =   hold[1]*hnew[2] - hold[2]*hnew[1] ;
      lrot[1] = -(hold[0]*hnew[2] - hold[2]*hnew[0]);
      lrot[2] =   hold[0]*hnew[1] - hold[1]*hnew[0] ;
      double lmag2 = lrot[0]*lrot[0] + lrot[1]*lrot[1] + lrot[2]*lrot[2];
      double lmag = sqrt(lmag2);
      double sinalpha = sqrt(lmag2 / (oldmag2*newmag2));
      if(lmag > 0)
	for(unsigned i=0; i<3; i++) lrot[i] /= sqrt(lmag);
      else continue;

      // get the rotation operator in the flavor basis
      double alpha = asin(sinalpha);
      complex<double> Rcoeff[4];
      for(int i=0; i<3; i++) Rcoeff[i] = -I * sin(alpha/2.) * lrot[i];
      Rcoeff[3] = cos(alpha/2.);
      
      MATRIX<complex<double>,NF,NF> R;
      pauli_reconstruct(Rcoeff, R);

      // apply to Scumulative
      s.Scumulative[m][i] = MATRIX<complex<double>,NF,NF>(Adjoint(s.U0[m][i]) * R * s.U0[m][i] * s.Scumulative[m][i]);
    }
  }
}

#endif
