#ifndef PROFILE_H
#define PROFILE_H

#include "H5Cpp.h"

class Profile{
 public:
  DISCONTINUOUS lnrho; // ln(g/ccm)
  DISCONTINUOUS Ye;
  DISCONTINUOUS temperature; // MeV
  DISCONTINUOUS Ecom_Elab;
  DISCONTINUOUS Elab_Elab0;
  array<array<array<DISCONTINUOUS,NF>,NE>,NM> Dens_unosc;
  array<array<array<DISCONTINUOUS,NF>,NE>,NM> Flux_unosc;
  array<array<array<DISCONTINUOUS,NF>,NE>,NM> Pres_unosc;
  array<double,NE> Ecom, Etopcom; // erg
  double rstart;

  Profile(string inputfile, double rmin, bool do_SR, bool do_GR){
    H5::H5File file(inputfile, H5F_ACC_RDONLY );

    // get the dimenions of the dataset
    hsize_t dims[3];
    int ndims = file.openDataSet("Ndens(1|ccm)").getSpace().getSimpleExtentDims(dims);
    assert(ndims == 3);
    const hsize_t ns = dims[0];
    const hsize_t ng = dims[1];
    const hsize_t nr = dims[2];
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
    Elab_Elab0.SetData(x, data);

    file.openDataSet("Ye"        ).read(&data[0], H5::PredType::NATIVE_DOUBLE);
    Ye.SetData(x, data);

    file.openDataSet("T(MeV)"    ).read(&data[0], H5::PredType::NATIVE_DOUBLE);
    temperature.SetData(x, data);

    file.openDataSet("rho(g|ccm)").read(&data[0], H5::PredType::NATIVE_DOUBLE);
    lnrho.SetData(x, data);
    lnrho = lnrho.copy_logy();
    
    // normalize the lab-frame neutrino energy relative to the start of the calculation
    rstart = max(rmin, lnrho.XMin());
    assert(rstart >= lnrho.XMin());
    assert(rstart <= lnrho.XMax());
    double startval = Elab_Elab0( rstart );
    for(size_t i=0; i<Elab_Elab0.data.size(); i++)
      Elab_Elab0.data[i] /= startval;

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
    			Dens_unosc[m][i][f].SetData(x,data);

    			for(size_t ir=0; ir<nr; ir++) data[ir] = Fdens[s][i][ir];
    			Flux_unosc[m][i][f].SetData(x,data);

    			for(size_t ir=0; ir<nr; ir++) data[ir] = Pdens[s][i][ir];
    			Pres_unosc[m][i][f].SetData(x,data);

    		}
    	}
    }

  }

};


#endif
