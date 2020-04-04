#ifndef PROFILE_H
#define PROFILE_H

#include "H5Cpp.h"

class Profile{
 public:
  DISCONTINUOUS rho; // ln(g/ccm)
  DISCONTINUOUS Ye;
  DISCONTINUOUS temperature; // MeV
  DISCONTINUOUS Ecom_Elab;
  DISCONTINUOUS Elab_Elabstart;
  array<array<array<DISCONTINUOUS,NF>,NE>,NM> dens_unosc;
  array<array<array<DISCONTINUOUS,NF>,NE>,NM> fluxfac_unosc;
  array<array<array<DISCONTINUOUS,NF>,NE>,NM> eddfac_unosc;
  array<double,NE> Ecom, Etopcom; // erg
  double rstart, Elabstart_Ecomstart;

  Profile(string inputfile, double rhostart, bool do_SR, bool do_GR){
    H5::H5File file(inputfile, H5F_ACC_RDONLY );

    // get the dimenions of the dataset
    hsize_t dims[3];
    int ndims = file.openDataSet("Ndens(1|ccm)").getSpace().getSimpleExtentDims(dims);
    assert(ndims == 3);
    const hsize_t ns = dims[0];
    const hsize_t ng = dims[1];
    const hsize_t nr = dims[2];
    cout<<ng<<endl;
    cout<<NE<<endl;
    assert(ng==NE);
    assert(ns==2*NF);
    vector<double> x(nr), data(nr);
  
    // energy grid
    file.openDataSet("Ecom(erg)"   ).read(&Ecom,    H5::PredType::NATIVE_DOUBLE);
    file.openDataSet("Etopcom(erg)").read(&Etopcom, H5::PredType::NATIVE_DOUBLE);

    // load rho and Ye data
    file.openDataSet("ct(cm)"    ).read(&x[0],    H5::PredType::NATIVE_DOUBLE);

    if(do_SR) file.openDataSet("Ecom_Elab" ).read(&data[0], H5::PredType::NATIVE_DOUBLE);
    else for(size_t i=0; i<data.size(); i++) data[i] = 1.;
    Ecom_Elab.SetData(x, data);

    if(do_GR) file.openDataSet("Elab_Elab0").read(&data[0], H5::PredType::NATIVE_DOUBLE);
    else for(size_t i=0; i<data.size(); i++) data[i] = 1.;
    Elab_Elabstart.SetData(x, data);

    file.openDataSet("Ye"        ).read(&data[0], H5::PredType::NATIVE_DOUBLE);
    Ye.SetData(x, data);

    file.openDataSet("T(MeV)"    ).read(&data[0], H5::PredType::NATIVE_DOUBLE);
    temperature.SetData(x, data);

    file.openDataSet("rho(g|ccm)").read(&data[0], H5::PredType::NATIVE_DOUBLE);
    rho.SetData(x, data);
    
    // get the starting radius
    rstart = rho.XMin();
    for(size_t i=0; i<rho.data.size(); i++){
      if(rho.data[i] >= rhostart){
	cout << i << " " << rho.data[i] << " " << rhostart << endl;
	double dr = rho.x[i+1] - rho.x[i];
	double drho = rho.data[i+1] - rho.data[i];
	if(abs(drho)>0 and i<rho.x.size()-1){
	  rstart = rho.x[i] + (rhostart-rho.data[i]) * dr/drho;
	}
	else rstart = rho.x[i];
      }
    }
    if(rstart==rho.XMin()) rstart = max(rho.XMin(),0.); // have reference point be at 0 if starting from boundary

    // normalize the lab-frame neutrino energy relative to the start of the calculation
    assert(rstart >= rho.XMin());
    assert(rstart <= rho.XMax());
    Elabstart_Ecomstart = 1./Ecom_Elab(rstart);
    double startval = Elab_Elabstart( rstart );
    for(size_t i=0; i<Elab_Elabstart.data.size(); i++)
      Elab_Elabstart.data[i] /= startval;

    // load and compute spectral data
    double Ndens[ns][ng][nr];
    double Fdens[ns][ng][nr];
    double Pdens[ns][ng][nr];

    file.openDataSet("Ndens(1|ccm)").read(Ndens, H5::PredType::NATIVE_DOUBLE);
    file.openDataSet("Fdens(1|ccm)").read(Fdens, H5::PredType::NATIVE_DOUBLE);
    file.openDataSet("Pdens(1|ccm)").read(Pdens, H5::PredType::NATIVE_DOUBLE);

    for(int m=matter; m<=antimatter; m++){
    	for(int i=0; i<NE; i++){
    		for(int f=e; f<=mu; f++){
    			int s = m + 2*f; // species index. 0-e 1-ebar 2-x 3-xbar

    			for(size_t ir=0; ir<nr; ir++) data[ir] = Ndens[s][i][ir];
			dens_unosc[m][i][f].SetData(x,data);

    			for(size_t ir=0; ir<nr; ir++) data[ir] = Fdens[s][i][ir]/Ndens[s][i][ir];
    			fluxfac_unosc[m][i][f].SetData(x,data);

    			for(size_t ir=0; ir<nr; ir++) data[ir] = Pdens[s][i][ir]/Ndens[s][i][ir];
    			eddfac_unosc[m][i][f].SetData(x,data);

    		}
    	}
    }

  }

};


#endif
