#ifndef EVOLVE_H
#define EVOLVE_H

void evolve_oscillations(State& s,
			 const State& s0,
			 const double r_end,
			 double& dr,
			 const DISCONTINUOUS& lnrho,
			 const DISCONTINUOUS& temperature,
			 const DISCONTINUOUS& Ye,
			 const array<array<array<DISCONTINUOUS,NF>,NE>,NM>& P_unosc,
			 const double accuracy,
			 const double increase,
			 const array<array<MATRIX<complex<double>,NF,NF>,NE>,NM>& HfV){
  array<array<array<array<array<double,NY>,NS>,NE>,NM>,NRK> dYdr;

  State sBeforeBlock = s;
  
  bool finish = false;
  do{ 
    if(s.r+dr>r_end){
      dr=r_end-s.r;
      finish=true;
    }

    // beginning of RK section
    const State sReset = s;
    bool repeat = false;
    double maxerror = 0;
    do{
      repeat=false;
      for(int k=0;k<=NRK-1;k++){
	s = sReset;
	s.r +=AA[k]* dr;

	for(state m = matter; m <= antimatter; m++)
	  for(int i=0;i<=NE-1;i++)
	    for(solution x=msw;x<=si;x++)
	      for(int j=0;j<=NY-1;j++)
		for(int l=0;l<=k-1;l++)
		  s.Y[m][i][x][j] += BB[k][l] * dYdr[l][m][i][x][j] * dr;

	s.update_potential(lnrho,temperature,Ye,P_unosc,HfV,s0);
	dYdr[k] = K(s);
      }
	  
      // increment all quantities and update C and A arrays
      s = sReset;
      s.r= sReset.r + dr;
      maxerror=0.;
      for(state m=matter;m<=antimatter;m++){
	for(int i=0;i<=NE-1;i++){
	  for(solution x=msw;x<=si;x++){
	    for(int j=0;j<=NY-1;j++){
	      double Yerror = 0.;
	      for(int k=0;k<=NRK-1;k++){
		assert(CC[k] == CC[k]);
		s.Y[m][i][x][j] += CC[k] * dYdr[k][m][i][x][j] * dr;
		Yerror += (CC[k]-DD[k]) * dYdr[k][m][i][x][j] * dr;
	      }
	      maxerror = max(maxerror, fabs(Yerror));
	    }
	  }
	}
      }
      s.update_potential(lnrho,temperature,Ye,P_unosc,HfV,s0);

      // decide whether to accept step, if not adjust step size
      if(maxerror>accuracy){
	dr *= 0.9 * pow(accuracy/maxerror, 1./(NRKOrder-1.));
	if(dr > 4.*s.r*numeric_limits<double>::epsilon()) repeat=true;
      }

      // reset integration variables to those at beginning of step
      if(repeat==true) s = sReset;
      else s.counter++;
      
    }while(repeat==true); // end of RK section

    // adjust step size based on RK error
    // could be moved up to RK section but better left here 
    // in case adjustments are necessary based on new S matrices
    dr = min(dr*pow(accuracy/maxerror,1./max(1,NRKOrder)),increase*dr);
    double drmin = 4.*s.r*numeric_limits<double>::epsilon();
    dr = max(dr,drmin);
  } while(!finish);

  
  // accumulate S and reset variables
  // need sBeforeBlock because we tack on to the end of Scumulative
  s.accumulate_S(dr, sBeforeBlock);
  s.assert_noNaN();
}

#endif
