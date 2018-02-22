void initialize(vector<vector<MATRIX<complex<double>,NF,NF> > >& fmatrixf,
	      double rho, double T, double Ye){
  for(int i=0; i<NE; i++){
    for(state m=matter; m<=antimatter; m++)
      for(flavour f1=e; f1<=mu; f1++)
	for(flavour f2=e; f2<=mu; f2++) 
	  fmatrixf[m][i][f1][f2] = 0;
  }
}

double get_rho(const double r){
  return exp(lnrho(log(r)));
}
double get_drhodr(const double rrho, const double r){
  return rrho*lnrho.Derivative(log(r))/r;
}
double get_Ye(const double r){
  return Ye(r);
}
double get_dYedr(const double r){
  return Ye.Derivative(r);
}

//=======//
// Ebins //
//=======//
void set_Ebins(vector<double>& E){
  const double NEP=8;
  cout<<"NE="<<NE << " NEP="<<NEP;
  E.resize(NE);
  for(int i=0;i<NE;i++){
    unsigned ind = i;
    if(NE==1||NE==2||NE==4) ind = i*NEP/NE+NEP/NE/2;

    E[i] = 1.e6*cgs::units::eV * exp(log(2*1000000) + ind*(log(37.48*1000000) - log(2*1000000))/(NE-1))/1000000;

    cout<<E[i]<<endl;
    cout.flush();
  }
  cout.flush();
}

//===================//
// Vacuum Potentials //
//===================//
double deltaV(const double E){ // erg
  return abs(dm21)*cgs::constants::c4 / (2.*E);
}

void set_kV(vector<vector<double> >& kV){
  assert(NF==2);
  for(int i=0;i<NE;i++){
    kV[i][0] = m1*m1 * cgs::constants::c4 /2./E[i];
    kV[i][1] = kV[i][0] + deltaV(E[i]);
  }
}

//===================//
// Matter Potentials //
//===================//
double Ve(double rho, double Ye){
  return (M_SQRT2*cgs::constants::GF/cgs::constants::Mp)*rho*Ye;
}

double dVedr(double rho, double drhodr, double Ye, double dYedr){
  return (M_SQRT2*cgs::constants::GF/cgs::constants::Mp) * (drhodr*Ye + rho*dYedr );
}

double Vmu(double rho, double Ye){ return 0.;}

double dVmudr(double rho, double drhodr, double Ye, double dYedr){ return 0.;}


//=============================//
// Self-Interaction Potentials //
//=============================//
void getP(const double r,
	  const vector<vector<MATRIX<complex<double>,NF,NF> > > U0, 
	  const vector<vector<MATRIX<complex<double>,NF,NF> > > Scumulative, 
	  vector<vector<MATRIX<complex<double>,NF,NF> > >& pmatrixf0,
	  vector<vector<MATRIX<complex<double>,NF,NF> > >& pmatrixm0){
  
  for(int i=0;i<=NE-1;i++){

    // set all elements to zero
    for(flavour f=e;f<=mu;f++) for(flavour fp=e; fp<=mu; fp++){
	pmatrixf0[matter    ][i][f][fp] = 0;
	pmatrixf0[antimatter][i][f][fp] = 0;
      }

    // read un-oscillated potentials from file
    pmatrixf0[matter    ][i][e ][e ]=eP[i](r);
    pmatrixf0[antimatter][i][e ][e ]=eBarP[i](r);
    pmatrixf0[antimatter][i][mu][mu]=xP[i](r);
    pmatrixf0[matter    ][i][mu][mu]=xP[i](r);

    // oscillate potentials and put in mass basis
    for(state m=matter;m<=antimatter;m++){
      pmatrixm0[m][i] = Scumulative[m][i]
	* Adjoint(U0[m][i])
	* pmatrixf0[m][i]
	* U0[m][i]
	* Adjoint(Scumulative[m][i]);
      pmatrixf0[m][i] = U0[m][i]
	* pmatrixm0[m][i]
	* Adjoint(U0[m][i]);
    }
  }
}

void interact(vector<vector<MATRIX<complex<double>,NF,NF> > >& fmatrixf,
	      double rho, double T, double Ye, double dr){
  return;
  // don't do anything if too sparse
  //if(log10(rho) <= __nulibtable_MOD_nulibtable_logrho_min)
  //  return;

  // T should be MeV
  // nspecies, ngroups, nvars
  //nulibtable_range_species_range_energy_(&rho, &T, &Ye, &eas.storage.front(),
  //					 &__nulibtable_MOD_nulibtable_number_species,
  //					 &__nulibtable_MOD_nulibtable_number_groups,
  //					 &__nulibtable_MOD_nulibtable_number_easvariables);
  //eas.fix_units();


  //double tmp = 0;
  //for(int i=0; i<NE; i++){
    //cout << i << " " << eas.emis(0,i) << " " << eas.abs(0,i) << endl;
    // scale the off-diagonal components
    //tmp = (eas.abs(0,i) + eas.abs(2,i)) * 0.5;
    //fmatrixf[    matter][i][e][e ] *= exp(tmp     * dr);
    //fmatrixf[    matter][i][e][mu] *= exp(tmp*0.5 * dr);
    //fmatrixf[    matter][i][mu][e] = conj(fmatrixf[matter][i][e][mu]);

    // scale the diagonal components
    //fmatrixf[    matter][i][e ][e ] += (eas.emis(0,i) - eas.abs(0,i)*fmatrixf[    matter][i][e ][e ]) * dr;
    //fmatrixf[    matter][i][mu][mu] += (eas.emis(2,i) - eas.abs(2,i)*fmatrixf[    matter][i][mu][mu]) * dr;
    //fmatrixf[antimatter][i][e ][e ] += (eas.emis(1,i) - eas.abs(1,i)*fmatrixf[antimatter][i][e ][e ]) * dr;
    //fmatrixf[antimatter][i][mu][mu] += (eas.emis(2,i) - eas.abs(2,1)*fmatrixf[antimatter][i][mu][mu]) * dr;
  //}
}
