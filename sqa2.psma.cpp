#include "sqa2.psma.h"

//======//
// MAIN //
//======//
int main(int argc, char *argv[]){
  try{ 
    int in=1;
    string inputfilename,potential_directory;
    ofstream fout,foutC,foutP,foutS, foutf;
    string outputfilename,rhofilename, Yefilename, vfilename, spectrapath, nulibfilename, temperaturefilename;
    string outputfilenamestem;
    string nt, note;
    
    inputfilename=string(argv[in]);
    ifstream fin(inputfilename.c_str());
    
    // load the nulib table
    fin>>nulibfilename;
    fin>>potential_directory;
    cout << nulibfilename << endl;
    nulib_init(nulibfilename, 0);

    fin>>rhofilename;
    fin>>Yefilename;
    fin>>temperaturefilename;
    //fin>>spectrapath;
    fin>>outputfilename;
    outputfilenamestem = outputfilename+"/";
    
    double rmin, rmax;
    fin>>rmin>>rmax; // cm
    
    fin>>NE>>Emin>>Emax; // MeV
    
    //fin>>Rnu; // cm
    //fin>>t;   // s
    
    m1=0.;
    fin>>dm21;
    fin>>theta12V;
    
    alphaV[0]=0.;
    alphaV[1]=0.;
    
    betaV[0]=0.;
    
    double accuracy;
    fin>>accuracy;
    fin>>note;
    int iNT,fNT;
    string id;
    // fin>>iNT;
    // fin>>fNT;
    // fin>>id;//number of tracer
    fin>>nt;
    
    cout<<"\n\n*********************************************************\n";
    cout<<"\nrho\t"<<rhofilename;
    cout<<"\nYe\t"<<Yefilename;
    cout<<"\nT\t"<<temperaturefilename;
    cout<<"\noutput\t"<<outputfilename;
    cout<<"\nrmin\t"<<rmin<<"\trmax\t"<<rmax;
    //cout<<"\nRnu\t"<<Rnu<<"\nt\t"<<t;
    
    cout << "\n\nNE\t" << NE << "\tEmin\t" << Emin << "\tEmax\t" << Emax;
    
    cout<<"\n\nm1\t"<<m1<<"\tdm21^2\t"<<dm21;
    cout<<"\ntheta12V\t"<<theta12V;
    cout<<"\nalpha1V\t"<<alphaV[0]<<"\talpha2V\t"<<alphaV[1];
    cout<<"\nbeta1V\t"<<betaV[0];
    
    cout<<"\naccuracy\t"<<accuracy<<"\n";
    cout.flush();
    
    // load rho and Ye data
    rho.Open(rhofilename,'#');
    Ye.Open(Yefilename,'#');
    //temperature.Open(Yefilename,'#'); // UNCOMMENT
    rmin=max(rmin,max(rho.XMin(),Ye.XMin()) );
    //rmin = max(rmin, temperature.XMin() ); //UNCOMMENT
    rmax=min(rmax,min(rho.XMax(),Ye.XMax()) );
    //rmax = min(rmax, temperature.XMax() ); //UNCOMMENT

    lnrho=rho;
    lnrho.TransformX(log);
    lnrho.TransformY(log);
    
    // load and compute spectral data    
    eP.resize(NE);
    eBarP.resize(NE);
    xP.resize(NE);
    
    // load and compute spectral data
    for(int i=0;i<=NE-1;i++){
      eP   [i].Open(potential_directory+"/v_potential1_"+patch::to_string(i+1)+"_"+patch::to_string(nt)+note+".txt",'#');
      eBarP[i].Open(potential_directory+"/v_potential2_"+patch::to_string(i+1)+"_"+patch::to_string(nt)+note+".txt",'#');
      xP   [i].Open(potential_directory+"/v_potential3_"+patch::to_string(i+1)+"_"+patch::to_string(nt)+note+".txt",'#');
    }

    // output filestreams: the arrays of ofstreams cannot use the vector container - bug in g++
    foutS.open((outputfilename+"/"
		+  "S"+patch::to_string(rmin)
		+  "-"+patch::to_string(rmax)
		+ "km"+patch::to_string(accuracy)
		+"ACC"+patch::to_string(NE)
		+ "NE"+patch::to_string(nt)
		+ note+".txt").c_str());
    foutS.precision(12);
    foutS.flush();
    foutP.open((outputfilename+"/"
		+string("2p")+patch::to_string(rmin)
		+         "-"+patch::to_string(rmax)
		+        "km"+patch::to_string(accuracy)
		+       "ACC"+patch::to_string(NE)
		+        "NE"+patch::to_string(nt)
		+ note+".txt").c_str());
    foutP.precision(12);
    foutP.flush();
    foutf.open((outputfilename+"/"+"f.dat").c_str());
    foutf.precision(12);
    foutf << "# 1:r ";
    for(int i=0; i<NE; i++)
      for(state m=matter; m<=antimatter; m++)
	for(flavour f1=e; f1<=mu; f1++)
	  for(flavour f2=e; f2<=mu; f2++) {
	    int istart = 2*( f2 + f1*2 + m*2*2 + i*2*2*2) + 2;
	    foutf << istart   << ":ie"<<i<<"m"<<m<<"f"<<f1<<f2<<"R\t";
	    foutf << istart+1 << ":ie"<<i<<"m"<<m<<"f"<<f1<<f2<<"I\t";
    }
    foutf << endl;
    foutf.flush();
    
    fout.open((outputfilename+"/"
	       +      patch::to_string(rmin)
	       +  "-"+patch::to_string(rmax)
	       + "km"+patch::to_string(accuracy)
	       +"ACC"+patch::to_string(NE)
	       + "NE"+patch::to_string(nt)
	       + note+".txt").c_str());
    fout.precision(12);
    fout.flush();
    string comma(","), colon(":"), dat(".dat");;
    stringstream filename, filenamestart;
    filenamestart<<outputfilenamestem;
    filenamestart<<colon<<dm21<<comma<<theta12V;
    filenamestart<<colon<<NE<<comma<<Emin<<comma<<Emax;

    ofstream foutPvsr[NE],foutFvsr[NE];
    ofstream fPvsE, fFvsE;
    
    for(int i=0;i<=NE-1;i++){
      filename << filenamestart.str() 
	       << string(":E=") 
	       << ((NE-1.-i)*Emin+i*Emax)/(NE-1.) 
	       << string("MeV.P.dat");
      foutPvsr[i].open((filename.str()).c_str());
      foutPvsr[i].precision(12);
      filename.str("");
      
      filename << filenamestart.str() 
	       << string(":E=") 
	       << ((NE-1.-i)*Emin+i*Emax)/(NE-1.)
	       <<string("MeV.F.dat");
      foutFvsr[i].open((filename.str()).c_str());
      foutFvsr[i].precision(12);
      filename.str("");
    }
    
    // unit conversion to cgs
    Emin *= 1.*mega*cgs::units::eV;
    Emax *= 1.*mega*cgs::units::eV;
    m1   *= 1.*cgs::units::eV/cgs::constants::c2;
    dm21 *= 1.*cgs::units::eV*cgs::units::eV/cgs::constants::c4;
    theta12V *= M_PI/180.;
    c12V = cos(theta12V);
    s12V = sin(theta12V);
    
    // *************************************************
    // set up global variables defined in parameters.h *
    // *************************************************

    // vectors of energies and vacuum eigenvalues
    E = vector<double>(NE);
    kV = vector<vector<double> >(NE,vector<double>(NF));
    vector<double> energybin(NE);
    energybin=Ebins(NE);
    cout<<"NE"<<NE;
    cout.flush();
    for (int i = 0;i<NE;i++){
      cout<<energybin[i]<<endl;
      cout.flush();
    }
    for(int i=0;i<=NE-1;i++){
      if(NE==1||NE==2|NE==4) E[i]=energybin[i*NEP/NE+NEP/NE/2];
      else E[i]=energybin[i];
      
      kV[i][0] = m1*m1 * cgs::constants::c4 /2./E[i];
      kV[i][1] = (m1*m1 + dm21) * cgs::constants::c4 /2./E[i];
    }
    
    // determine eigenvalue ordering
    ordering[0]=0; ordering[1]=1;
    vector<double> tempkV(kV[0]);
    Sort(tempkV,ordering,ascending);
    
    if(kV[0][1]>kV[0][0]){
      a1=-1;
      a2=+1;
      cout<<"\n\nNormal hierarchy";
    }
    else{ 
      if(kV[0][1]<kV[0][0]){
	a1=+1; 
	a2=-1; 
	cout<<"\n\nInverted hierarchy";
      }
      else{ 
	cout<<endl<<endl<<"Neither normal or Inverted"<<endl; 
	abort();
      }
    }
    
    // vaccum mixing matrices and Hamiltonians
    Evaluate_UV();
    
    HfV[matter] = vector<MATRIX<complex<double>,NF,NF> >(NE);
    HfV[antimatter] = vector<MATRIX<complex<double>,NF,NF> >(NE);
    Evaluate_HfV();
    
    // cofactor matrices in vacuum
    CV=vector<vector<MATRIX<complex<double>,NF,NF> > >(NE,vector<MATRIX<complex<double>,NF,NF> >(NF));
    Evaluate_CV();
    
    // mixing matrix element prefactors in vacuum
    AV=vector<vector<vector<double> > >(NE,vector<vector<double> >(NF,vector<double>(NF)));
    Evaluate_AV();
    
    // **************************************
    // quantities evaluated at inital point *
    // **************************************
    
    // MSW potential matrix
    double rrho = exp(lnrho(log(rmin)));
    double YYe=Ye(rmin);
    
    MATRIX<complex<double>,NF,NF> VfMSW0, Hf0;
    vector<double> k0, deltak0;
    
    VfMSW0[e][e]=Ve(rrho,YYe);
    VfMSW0[mu][mu]=Vmu(rrho,YYe);
    
    // cofactor matrices at initial point - will be recycled as cofactor matrices at beginning of every step
    vector<vector<vector<MATRIX<complex<double>,NF,NF> > > > 
      C0(NM,vector<vector<MATRIX<complex<double>,NF,NF> > >(NE,vector<MATRIX<complex<double>,NF,NF> >(NF)));

    // mixing matrix element prefactors at initial point - will be recycled like C0
    vector<vector<vector<vector<double> > > > 
      A0(NM,vector<vector<vector<double> > >(NE,vector<vector<double> >(NF,vector<double>(NF))));
    
    // mixing angles to MSW basis at initial point
    U0[matter] = vector<MATRIX<complex<double>,NF,NF> >(NE);
    U0[antimatter] = vector<MATRIX<complex<double>,NF,NF> >(NE);
    
    for(int i=0;i<=NE-1;i++){
      Hf0=HfV[matter][i]+VfMSW0;
      k0=k(Hf0);
      deltak0=deltak(Hf0);
      C0[matter][i]=CofactorMatrices(Hf0,k0);
      
      for(int j=0;j<=NF-1;j++){
	if(real(C0[matter][i][j][mu][e]*CV[i][j][mu][e]) < 0.)
	  A0[matter][i][j][e]=-AV[i][j][e];
	else A0[matter][i][j][e]=AV[i][j][e];
	A0[matter][i][j][mu]=AV[i][j][mu];
      }
      U0[matter][i]=U(deltak0,C0[matter][i],A0[matter][i]);
      
      Hf0=HfV[antimatter][i]-VfMSW0;
      k0=kbar(Hf0);
      deltak0=deltakbar(Hf0);
      C0[antimatter][i]=CofactorMatrices(Hf0,k0);
      for(int j=0;j<=NF-1;j++){
	if(real(C0[antimatter][i][j][mu][e]*CV[i][j][mu][e]) < 0.)
	  A0[antimatter][i][j][e]=-AV[i][j][e];
	else A0[antimatter][i][j][e]=AV[i][j][e];
	A0[antimatter][i][j][mu]=AV[i][j][mu];
      }
      U0[antimatter][i]=Conjugate(U(deltak0,C0[antimatter][i],A0[antimatter][i]));
    }
    
    vector<vector<vector<double> > > P0 (NM,vector<vector<double> >(NF,vector <double>(NE)));

    // density matrices at initial point, rhomatrixm0 - but not rhomatrixf0
    // will be updated whenever discontinuities are crossed and/or S is reset
    pmatrixf0[matter]=pmatrixf0[antimatter]=vector<MATRIX<complex<double>,NF,NF> >(NE);
    pmatrixm0[matter]=pmatrixm0[antimatter]=vector<MATRIX<complex<double>,NF,NF> >(NE);
    fmatrixf[matter]=fmatrixf[antimatter]=vector<MATRIX<complex<double>,NF,NF> >(NE);
    fmatrixm[matter]=fmatrixm[antimatter]=vector<MATRIX<complex<double>,NF,NF> >(NE);

    // yzhu14 density/potential matrices art rmin
    getP(rmin);
    double rho0 = rho(rmin);
    double T0 = 5.0;//temperature(rmin);
    double ye0 = Ye(rmin);
    initialize(fmatrixf,rho0,T0,ye0);

    // ***************************************
    // quantities needed for the calculation *
    // ***************************************
    double r,r0,dr,drmin;
    int ND;
    vector<double> rs;
    
    double maxerror,increase=3.;
    bool repeat, finish, resetflag, output;
    int counterout,step;
    
    // comment out if not following as a function of r
    fin>>step;
    
    // self-interaction integration factor
    NSI=M_SQRT2*cgs::constants::GF*(Emax-Emin)/(NE-1.);
    
    // ***************************************
    // variables followed as a function of r *
    // ***************************************
    
    vector<vector<vector<vector<double> > > > 
      Y(NM,vector<vector<vector<double> > >(NE,vector<vector<double> >(NS,vector<double>(NY))));
    vector<vector<vector<vector<double> > > > 
      dY(NM,vector<vector<vector<double> > >(NE,vector<vector<double> >(NS,vector<double>(NY))));
    vector<vector<vector<vector<double> > > > 
      Y0(NM,vector<vector<vector<double> > >(NE,vector<vector<double> >(NS,vector<double>(NY))));
    vector<vector<vector<vector<double> > > > 
      Yerror(NM,vector<vector<vector<double> > >(NE,vector<vector<double> >(NS,vector<double>(NY))));
    
    // cofactor matrices
    vector<vector<vector<MATRIX<complex<double>,NF,NF> > > > C=C0;;
    // mixing matrix prefactors
    vector<vector<vector<vector<double> > > > A=A0;
    
    // accumulated S matrices from prior integration domains
    vector<vector<MATRIX<complex<double>,NF,NF> > > 
      Scumulative(NM,vector<MATRIX<complex<double>,NF,NF> >(NE,UnitMatrix<complex<double> >(NF)));
    
    // ************************
    // Runge-Kutta quantities *
    // ************************
    int NRK,NRKOrder;
    const double *AA=NULL,**BB=NULL,*CC=NULL,*DD=NULL;
    RungeKuttaCashKarpParameters(NRK,NRKOrder,AA,BB,CC,DD);
    
    vector<vector<vector<vector<vector<double> > > > > 
      Ks(NRK,vector<vector<vector<vector<double> > > >
	 (NM,vector<vector<vector<double> > >(NE,vector<vector<double> >(NS,vector<double>(NY)))));
    
    // temporaries
    MATRIX<complex<double>,NF,NF> SSMSW,SSSI,SThisStep;
    
    // *********************
    // integration domians *
    // *********************
    
    rho.FindDomains();
    
    rs.push_back(rmin);
    for(int d=1;d<=rho.NDomains();d++){
      r=rho.Discontinuity(d);
      if(r>rmin && r<rmax) rs.push_back(r);
    }
    rs.push_back(rmax);
    
    sort(rs.begin(),rs.end());
    ND=rs.size()-1;
    
    // **********************
    // start of calculation *
    // **********************
    
    for(int d=0;d<=ND-1;d++){
      if(d==0) rmin=rs.front();
      else rmin=rs[d]+1.*cgs::units::cm;
      
      if(d==ND-1) rmax=rs.back();
      else rmax=rs[d+1]-1.*cgs::units::cm;
      
      cout<<"\n"<<d<<"\t"<<rmin<<"\t"<<rmax << endl;
      cout << endl;
      cout << "r(km)  dr(cm)" << endl;
      cout.flush();
      
      // oscillate pmatrix for background potential
      for(int i=0;i<=NE-1;i++){
	pmatrixm0[matter][i] = Scumulative[matter][i]
	  * Adjoint(U0[matter][i])
	  * pmatrixf0[matter][i]
	  * U0[matter][i]
	  * Adjoint(Scumulative[matter][i]);
	
	pmatrixm0[antimatter][i] = Scumulative[antimatter][i]
	  * Adjoint(U0[antimatter][i])
	  * pmatrixf0[antimatter][i]
	  * U0[antimatter][i]
	  * Adjoint(Scumulative[antimatter][i]);
      }

      // *****************************************
      // initialize at beginning of every domain *
      // *****************************************
      r=rmin;
      dr=1e-3*cgs::units::cm;
      drmin=4.*r*numeric_limits<double>::epsilon();
      
      for(state m=matter;m<=antimatter;m++)
	for(int i=0;i<=NE-1;i++){
	  Y[m][i][msw][0] = Y[m][i][si][0] = M_PI/2.;
	  Y[m][i][msw][1] = Y[m][i][si][1] = M_PI/2.;
	  Y[m][i][msw][2] = Y[m][i][si][2] = 0.;
	  Y[m][i][msw][3] = Y[m][i][si][3] = 1.; // The determinant of the S matrix
	  Y[m][i][msw][4] = Y[m][i][si][4] = 0.;
	  Y[m][i][msw][5] = Y[m][i][si][5] = 0.;
	}
      
      // *************************************************
      // comment out if not following as a function of r *
      // *************************************************
	
      finish=output=false;
      counterout=1;
      Outputvsr(fout,foutP,foutf,foutPvsr,foutFvsr,r,Y,C,A,Scumulative);
	
      // ***********************
      // start the loop over r *
      // ***********************
      do{ 
	double intkm = int(r/1e5)*1e5;
	if(r - intkm <= dr){
	  cout << r/1e5 << " " << dr << endl;
	  cout.flush();
	}

	if(r+dr>rmax){
	  dr=rmax-r;
	  finish=true;
	  output=true;
	}
	  
	r0=r;
	Y0=Y;
	C0=C;
	A0=A;
	  
	// beginning of RK section
	do{ 
	  repeat=false;
	  // first step: assumes derivatives are evaluated at r
	  getP(r);

	  for(int i=0;i<=NE-1;i++)
	    for(state m=matter;m<=antimatter;m++)
		pmatrixm0[m][i] = Scumulative[m][i]
		  * Adjoint(U0[m][i])
		  * pmatrixf0[m][i]
		  * U0[m][i]
		  * Adjoint( Scumulative[m][i] );

	  // pmatrixmm=pmatrixm0[matter];
	  // pmatrixmmbar=pmatrixm0[antimatter];
	  K(r,dr,Y,C,A,Ks[0]);
	    
	  // second step
	  r=r0+AA[1]*dr;
	  for(state m=matter;m<=antimatter;m++)
	    for(int i=0;i<=NE-1;i++)
	      for(solution x=msw;x<=si;x++)
		for(int j=0;j<=NY-1;j++)
		  Y[m][i][x][j] += BB[1][0] * Ks[0][m][i][x][j];

	  getP(r);     
	  for(int i=0;i<=NE-1;i++)
	    for(state m=matter;m<=antimatter;m++)
		pmatrixm0[m][i] = Scumulative[m][i]
		  * Adjoint(U0[m][i])
		  * pmatrixf0[m][i]
		  * U0[m][i]
		  * Adjoint( Scumulative[m][i] );

	  //pmatrixmm=pmatrixm0[matter];
	  //pmatrixmmbar=pmatrixm0[antimatter];
	  K(r,dr,Y,C,A,Ks[1]);
	    
	  // remaining steps
	  for(int k=2;k<=NRK-1;k++){
	    r=r0+AA[k]*dr;
	    Y=Y0;

	    for(state m = matter; m <= antimatter; m++)
	      for(int i=0;i<=NE-1;i++)
		for(solution x=msw;x<=si;x++)
		  for(int j=0;j<=NY-1;j++)
		    for(int l=0;l<=k-1;l++)
		      Y[m][i][x][j] += BB[k][l] * Ks[l][m][i][x][j];
	    
	    getP(r);     
	    for(int i=0;i<=NE-1;i++)
	      for(state m=matter;m<=antimatter;m++)
		  pmatrixm0[m][i] = Scumulative[m][i]
		    * Adjoint(U0[m][i])
		    * pmatrixf0[m][i]
		    * U0[m][i]
		    * Adjoint( Scumulative[m][i] );

	    K(r,dr,Y,C,A,Ks[k]);
	  }
	  
	  // increment all quantities and update C and A arrays
	  r=r0+dr;
	  for(state m=matter;m<=antimatter;m++){
	    for(int i=0;i<=NE-1;i++){
	      for(solution x=msw;x<=si;x++){
		for(int j=0;j<=NY-1;j++){
		  Y[m][i][x][j] = Y0[m][i][x][j];
		  Yerror[m][i][x][j] = 0.;
		  for(int k=0;k<=NRK-1;k++){
		    Y[m][i][x][j] += CC[k] * Ks[k][m][i][x][j];
		    Yerror[m][i][x][j] += (CC[k]-DD[k]) * Ks[k][m][i][x][j];
		  }
		}
	      }
	    }
	  }
	  
	  C=UpdateC(r,Ye(r));
	  A=UpdateA(C,C0,A0);
	    
	  // find largest error
	  maxerror=0.;
	  for(state m=matter;m<=antimatter;m++)
	    for(int i=0;i<=NE-1;i++)
	      for(solution x=msw;x<=si;x++)
		for(int j=0;j<=NY-1;j++)
		  maxerror = max( maxerror, fabs(Yerror[m][i][x][j]) );
	    
	  // decide whether to accept step, if not adjust step size
	  if(maxerror>accuracy){
	    dr *= 0.9 * pow(accuracy/maxerror, 1./(NRKOrder-1.));
	    if(dr > drmin) repeat=true;
	  }

	  // reset integration variables to those at beginning of step
	  if(repeat==true){
	    r=r0;
	    Y=Y0;
	    C=C0;
	    A=A0;
	  }
	  
	}while(repeat==true); // end of RK section

	// interact with the matter
	interact(fmatrixf, rho0/*rho(r)*/, T0/*temperature(r)*/, ye0/*Ye(r)*/, dr);

	// accumulate S and reset variables
	for(state m=matter;m<=antimatter;m++){
	  for(int i=0;i<=NE-1;i++){
	    SSMSW = W(Y[m][i][msw])*B(Y[m][i][msw]);
	    SSSI  = W(Y[m][i][si ])*B(Y[m][i][si ]);
	    SThisStep = SSMSW*SSSI;	
	    Scumulative[m][i]=MATRIX<complex<double>,NF,NF>(SThisStep*Scumulative[m][i] );
	      
	    // convert fmatrix from flavor basis to mass basis
	    // oscillate fmatrix in mass basis
	    // convert back to flavor basis.
	    // don't need to modify pmatrix since it's re-read at each timestep
	    fmatrixm[m][i] = SThisStep
	      * Adjoint( U0[m][i] ) 
	      * fmatrixf[m][i] 
	      * U0[m][i]
	      * Adjoint( SThisStep );
	    fmatrixf[m][i] = U0[m][i]
	      * fmatrixm[m][i] 
	      * Adjoint(  U0[m][i] );
	    
	    // reset the evolution matrix to identity
	    Y[m][i][msw][0]=Y[m][i][msw][1]=M_PI/2.;
	    Y[m][i][msw][2]=0.;
	    Y[m][i][msw][3]=1.;
	    Y[m][i][msw][4]=Y[m][i][msw][5]=0.;
	    
	    Y[m][i][si][0]=Y[m][i][si][1]=M_PI/2.;
	    Y[m][i][si][2]=0.;
	    Y[m][i][si][3]=1.;
	    Y[m][i][si][4]=Y[m][i][si][5]=0.;
	  }
	}

	// comment out if not following as a function of r
	if(counterout==step){
	  output=true;
	  counterout=1;
	}
	else counterout++;
	
	if(output==true || finish==true){
	  Outputvsr(fout,foutP,foutf,foutPvsr,foutFvsr,r,Y,C,A,Scumulative);
	  output=false;
	}

	// adjust step size based on RK error
	// could be moved up to RK section but better left here 
	// in case adjustments are necessary based on new S matrices
	dr = min(dr*pow(accuracy/maxerror,1./max(1,NRKOrder)),increase*dr);
	drmin = 4.*r*numeric_limits<double>::epsilon();
	dr = max(dr,drmin);

      } while(finish==false);

      // if this is not the last domain then carry the S matrix across the discontinuities
      if(d<=ND-2){
	double rminus=rs[d+1]-1.*cgs::units::cm;
	double rplus=rs[d+1]+1.*cgs::units::cm;
	Scumulative = UpdateSm(rminus,rplus,Ye(rminus),Ye(rplus),Y,C,A,Scumulative);
	C0=C;
	C=UpdateC(rplus,Ye(rplus));
	A=UpdateA(C,C0,A0);
      }
      else{ // output at the end of the code
	Outputvsr(fout,foutP,foutf,foutPvsr,foutFvsr,rs[d+1],Y,C,A,Scumulative);
	//OutputvsE(fPvsE,fFvsE,rs[d+1],Y,C,A,Scumulative);
      }

    }// end of r loop

    for(int i=0;i<=NE-1;i++){
      foutPvsr[i].close();
      foutFvsr[i].close();
    }
    fPvsE.close();
    fFvsE.close();

  }
  catch(OUT_OF_RANGE<unsigned int> OOR){ OOR.Message();}
  catch(OUT_OF_RANGE<int> OOR){ OOR.Message();}
  catch(OUT_OF_RANGE<double> OOR){ OOR.Message();}
  catch(EMPTY E){ E.Message();}
  catch(FUNCTION_ERROR FE){ FE.Message();}
  catch(BASIC_ERROR BE){ BE.Message();}
  catch(...){ UNKNOWN_ERROR("main");}

  cout<<"\nFinished\n\a"; cout.flush();

  return 0;
}
