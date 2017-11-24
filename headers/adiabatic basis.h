

MATRIX<complex<double>,NF,NF> W(vector<double> Y);

vector<double> Q(MATRIX<complex<double>,NF,NF> Hf,MATRIX<complex<double>,NF,NF> dHfdr,MATRIX<complex<double>,NF,NF> U,vector<double> k,vector<double> deltak,vector<double> dkdr);
vector<double> Q(MATRIX<complex<double>,NF,NF> U,vector<double> deltak,vector<MATRIX<complex<double>,NF,NF> > C,vector<MATRIX<complex<double>,NF,NF> > dCdr);

MATRIX<complex<double>,NF,NF> Gamma(MATRIX<complex<double>,NF,NF> dHfdr,MATRIX<complex<double>,NF,NF> U,vector<double> k,vector<double> Q);

// *********************************************************************
// *********************************************************************

MATRIX<complex<double>,NF,NF> W(vector<double> Y)
       { MATRIX<complex<double>,NF,NF> w; 
         w[0][0]=exp(-I*M_2PI*Y[4]); w[1][1]=exp(-I*M_2PI*Y[5]);
         return w;
        }

// ********************************************************************************

vector<double> Q(MATRIX<complex<double>,NF,NF> Hf,MATRIX<complex<double>,NF,NF> dHfdr,MATRIX<complex<double>,NF,NF> U,vector<double> k,vector<double> deltak,vector<double> dkdr)
           { vector<double> Q(NF);

             flavour beta;
             double Xi;

             for(int i=0;i<=NF-1;i++)
                { beta=e;
                  if(norm(U[mu][i])>norm(U[beta][i])){ beta=mu;}

                  Xi=0.5/(norm(U[beta][i])*pow(deltak[0],2.)) *cgs::constants::hbarc; 

                  if(beta==e){ Q[i]=real(-I*Xi*( C<mu,e>(Hf,k[i])*dCdr<e,mu>(dHfdr,dkdr[i]) -C<e,mu>(Hf,k[i])*dCdr<mu,e>(dHfdr,dkdr[i]) ));}
                  if(beta==mu){ Q[i]=real(-I*Xi*( C<e,mu>(Hf,k[i])*dCdr<mu,e>(dHfdr,dkdr[i]) -C<mu,e>(Hf,k[i])*dCdr<e,mu>(dHfdr,dkdr[i]) ));}
                 }

             return Q;
            }

vector<double> Q(MATRIX<complex<double>,NF,NF> U,vector<double> deltak,vector<MATRIX<complex<double>,NF,NF> > C,vector<MATRIX<complex<double>,NF,NF> > dCdr)
           { vector<double> Q(NF);

             flavour beta;
             double Xi;

             for(int i=0;i<=NF-1;i++)
                { beta=e;
                  if(norm(U[mu][i])>norm(U[beta][i])){ beta=mu;}

                  Xi=0.5/(norm(U[beta][i])*pow(deltak[0],2.)) *cgs::constants::hbarc;

                  if(beta==e){ Q[i]=real(-I*Xi*( C[i][mu][e]*dCdr[i][e][mu] -C[i][e][mu]*dCdr[i][mu][e] ));}
                  if(beta==mu){ Q[i]=real(-I*Xi*( C[i][e][mu]*dCdr[i][mu][e] -C[i][mu][e]*dCdr[i][e][mu] ));}
                 }

             return Q;
            }

// *****************************************************************************************

MATRIX<complex<double>,NF,NF> Gamma(MATRIX<complex<double>,NF,NF> dHfdr,MATRIX<complex<double>,NF,NF> U,vector<double> k,vector<double> Q)
           { MATRIX<complex<double>,NF,NF> Gamma;
             MATRIX<complex<double>,NF,NF> UdagdHdxU=Adjoint(U)*dHfdr*U;

             for(int i=0;i<=NF-1;i++)
                { for(int j=1;j<=NF-1;j++)
                     { Gamma[i][j]=M_2PI*exp(I*(alphaV[i]-alphaV[j])) / (k[i]-k[j]) / (k[i]-k[j]+Q[i]-Q[j]) * UdagdHdxU[i][j] *cgs::constants::hbarc;
                       Gamma[j][i]=conj(Gamma[i][j]);
                      }
                 }

             return Gamma;
            }





