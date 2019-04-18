#ifndef OUTPUT_H
#define OUTPUT_H
#include <string>
#include "H5Cpp.h"
#include "State.h"
using namespace std;

void load_input_data(string inputfile,
		     DISCONTINUOUS& lnrho,
		     DISCONTINUOUS& Ye,
		     DISCONTINUOUS& temperature,
		     array<array<array<DISCONTINUOUS,NF>,NE>,NM>& P_unosc,
		     array<array<array<DISCONTINUOUS,NF>,NE>,NM>& D_unosc){
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
  
  // load rho and Ye data
  file.openDataSet("ct(cm)"    ).read(&x[0],    H5::PredType::NATIVE_DOUBLE);
  file.openDataSet("Ye"        ).read(&data[0], H5::PredType::NATIVE_DOUBLE);
  Ye.SetData(x, data);
  file.openDataSet("T(MeV)"    ).read(&data[0], H5::PredType::NATIVE_DOUBLE);
  temperature.SetData(x, data);
  file.openDataSet("rho(g|ccm)").read(&data[0], H5::PredType::NATIVE_DOUBLE);
  lnrho.SetData(x, data);
  lnrho = lnrho.copy_logy();
    
  // load and compute spectral data
  double Ndens[ns][ng][nr];
  double flux_factor[ns][ng][nr];
  double eddington_factor[ns][ng][nr];
  double pot_coeff = sqrt(2)*cgs::constants::GF;
  file.openDataSet("Ndens(1|ccm)"    ).read(Ndens,           H5::PredType::NATIVE_DOUBLE);
  file.openDataSet("flux_factor"     ).read(flux_factor,     H5::PredType::NATIVE_DOUBLE);
  file.openDataSet("eddington_factor").read(eddington_factor,H5::PredType::NATIVE_DOUBLE);
  for(int m=matter; m<=antimatter; m++){
    for(int i=0; i<NE; i++){
      for(int f=e; f<=mu; f++){
	int s = m + 2*f; // species index. 0-e 1-ebar 2-x 3-xbar
	
	for(size_t ir=0; ir<nr; ir++) data[ir] = Ndens[s][i][ir];
	D_unosc[m][i][f].SetData(x,data);
	for(size_t ir=0; ir<nr; ir++) data[ir] *= (1.-flux_factor[s][i][ir]) * pot_coeff;
	P_unosc[m][i][f].SetData(x,data);
      }
    }
  }

}

class FilePointers{
 public:
  hid_t file;
  hid_t dset_f, dset_S, dset_U, dset_r, dset_dr_osc, dset_dr_int, dset_dr_block,
    dset_rho, dset_Ye, dset_T, dset_VfSI, dset_dtdrO, dset_dtdrI, dset_dpdrO, dset_dpdrI;
  const hsize_t       dims[6]   = {0,             NM, NE, NF, NF, 2};
  const hsize_t   max_dims[6]   = {H5S_UNLIMITED, NM, NE, NF, NF, 2};
  const hsize_t chunk_dims[6]   = {1,             NM, NE, NF, NF, 2};
  const hsize_t     dims_V[5]   = {0,             NM, NF, NF, 2};
  const hsize_t max_dims_V[5]   = {H5S_UNLIMITED, NM, NF, NF, 2};
  const hsize_t chunk_dims_V[5] = {1,             NM, NF, NF, 2};
  const hsize_t     dims3[3]    = {0,             NM, NE};
  const hsize_t max_dims3[3]    = {H5S_UNLIMITED, NM, NE};
  const hsize_t chunk_dims3[3]  = {1,             NM, NE};
};

//============//
// setup_file //
//============//
FilePointers setup_HDF5_file(const array<double,NE>& E, const array<double,NE>& Vphase){
  FilePointers fp;
  
  fp.file = H5Fcreate("output.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  hid_t file_space;
  hsize_t ndims;
  hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_layout(plist, H5D_CHUNKED);

  // FMATRIXF //
  ndims = 6;
  file_space = H5Screate_simple(ndims, fp.dims, fp.max_dims);
  H5Pset_chunk(plist, ndims, fp.chunk_dims);
  H5Dcreate(fp.file, "fmatrixf", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Dcreate(fp.file, "S",        H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Dcreate(fp.file, "U",        H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  fp.dset_f = H5Dopen(fp.file, "fmatrixf", H5P_DEFAULT);
  fp.dset_S = H5Dopen(fp.file, "S", H5P_DEFAULT);
  fp.dset_U = H5Dopen(fp.file, "U", H5P_DEFAULT);

  // VfSI //
  ndims = 5;
  file_space = H5Screate_simple(ndims, fp.dims_V, fp.max_dims_V);
  H5Pset_chunk(plist, ndims, fp.chunk_dims_V);
  H5Dcreate(fp.file, "VfSI(erg)",H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  fp.dset_VfSI = H5Dopen(fp.file, "VfSI(erg)", H5P_DEFAULT);

  // dangle_dr //
  ndims = 3;
  file_space = H5Screate_simple(ndims, fp.dims3, fp.max_dims3);
  H5Pset_chunk(plist, ndims, fp.chunk_dims3);
  H5Dcreate(fp.file, "dtheta_dr_osc(rad|cm)",H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Dcreate(fp.file, "dtheta_dr_int(rad|cm)",H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Dcreate(fp.file, "dphi_dr_osc(rad|cm)",H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Dcreate(fp.file, "dphi_dr_int(rad|cm)",H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  fp.dset_dtdrO = H5Dopen(fp.file, "dtheta_dr_osc(rad|cm)", H5P_DEFAULT);
  fp.dset_dtdrI = H5Dopen(fp.file, "dtheta_dr_int(rad|cm)", H5P_DEFAULT);
  fp.dset_dpdrO = H5Dopen(fp.file, "dphi_dr_osc(rad|cm)", H5P_DEFAULT);
  fp.dset_dpdrI = H5Dopen(fp.file, "dphi_dr_int(rad|cm)", H5P_DEFAULT);

  
  // RADIUS/TIME //
  ndims = 1;
  file_space = H5Screate_simple(ndims, fp.dims, fp.max_dims);
  H5Pset_chunk(plist, ndims, fp.chunk_dims);
  H5Dcreate(fp.file, "r(cm)", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Dcreate(fp.file, "dr_block(cm)", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Dcreate(fp.file, "dr_osc(cm)", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Dcreate(fp.file, "dr_int(cm)", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Dcreate(fp.file, "rho(g|ccm)", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Dcreate(fp.file, "Ye", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Dcreate(fp.file, "T(MeV)", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  fp.dset_r        = H5Dopen(fp.file, "r(cm)", H5P_DEFAULT);
  fp.dset_dr_osc   = H5Dopen(fp.file, "dr_osc(cm)", H5P_DEFAULT);
  fp.dset_dr_int   = H5Dopen(fp.file, "dr_int(cm)", H5P_DEFAULT);
  fp.dset_dr_block = H5Dopen(fp.file, "dr_block(cm)", H5P_DEFAULT);
  fp.dset_rho      = H5Dopen(fp.file, "rho(g|ccm)", H5P_DEFAULT);
  fp.dset_Ye       = H5Dopen(fp.file, "Ye", H5P_DEFAULT);
  fp.dset_T        = H5Dopen(fp.file, "T(MeV)", H5P_DEFAULT);
  
  
  // energy grid //
  hid_t mem_space =   H5Screate_simple(ndims, &fp.dims[2], NULL);
  file_space = H5Screate_simple(ndims, &fp.dims[2], &fp.dims[2]);
  hid_t dset_Egrid = H5Dcreate(fp.file, "Egrid(erg)", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Dwrite(dset_Egrid, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &E[0]);
  dset_Egrid = H5Dcreate(fp.file, "Vphase(1|ccm)", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Dwrite(dset_Egrid, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &Vphase[0]);

  // clear resources
  H5Sclose(file_space);
  H5Pclose(plist);
  H5Fflush(fp.file,H5F_SCOPE_LOCAL);

  // preliminary stdout
  cout << "r(cm) \t rho(g/ccm) \t Ye \t T(MeV) \t dr_osc(cm) \t dr_int(cm) \t dr_block(cm)" << endl;
  
  return fp;
}

//============//
// write_data //
//============//
void write_data_HDF5(FilePointers& fp, const State& s, double dr_osc, double dr_int, double dr_block){
  // output to stdout
  cout << s.r/1e5 << "\t";
  cout << s.rho << "\t";
  cout << s.Ye << "\t";
  cout << s.T << "\t";
  cout << dr_osc << "\t";
  cout << dr_int << "\t";
  cout << dr_block << "\t";
  /* cout << impact << endl; */
  cout << endl;
  cout.flush();

  hid_t mem_space, file_space;
  hsize_t ndims;

  // create the memory space
  ndims = 6;
  mem_space = H5Screate_simple(ndims, fp.chunk_dims, NULL);
  file_space = H5Dget_space (fp.dset_f);
  hsize_t dims[ndims];
  H5Sget_simple_extent_dims(file_space, dims, NULL);
  hsize_t start[ndims] = {dims[0], 0, 0, 0, 0, 0};
  dims[0]++;

  // fmatrixf
  H5Dset_extent(fp.dset_f, dims);
  file_space = H5Dget_space(fp.dset_f);
  H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, fp.chunk_dims, NULL);
  H5Dwrite(fp.dset_f, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &s.fmatrixf);

  // S
  H5Dset_extent(fp.dset_S, dims);
  file_space = H5Dget_space(fp.dset_S);
  H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, fp.chunk_dims, NULL);
  H5Dwrite(fp.dset_S, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &s.Sf);

  // U
  H5Dset_extent(fp.dset_U, dims);
  file_space = H5Dget_space(fp.dset_U);
  H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, fp.chunk_dims, NULL);
  H5Dwrite(fp.dset_U, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &s.UU);

  // VfSI
  ndims = 5;
  mem_space = H5Screate_simple(ndims, fp.chunk_dims_V, NULL);
  file_space = H5Dget_space (fp.dset_VfSI);
  H5Sget_simple_extent_dims(file_space, dims, NULL);
  dims[0]++;
  H5Dset_extent(fp.dset_VfSI, dims);
  file_space = H5Dget_space (fp.dset_VfSI);
  H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, fp.chunk_dims_V, NULL);
  H5Dwrite(fp.dset_VfSI, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &s.VfSI);

  // dangle_dr
  ndims = 3;
  mem_space = H5Screate_simple(ndims, fp.chunk_dims3, NULL);
  file_space = H5Dget_space (fp.dset_dtdrO);
  H5Sget_simple_extent_dims(file_space, dims, NULL);
  dims[0]++;

  H5Dset_extent(fp.dset_dtdrO, dims);
  file_space = H5Dget_space (fp.dset_dtdrO);
  H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, fp.chunk_dims3, NULL);
  H5Dwrite(fp.dset_dtdrO, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &s.dtheta_dr_osc);
  
  H5Dset_extent(fp.dset_dtdrI, dims);
  file_space = H5Dget_space (fp.dset_dtdrI);
  H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, fp.chunk_dims3, NULL);
  H5Dwrite(fp.dset_dtdrI, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &s.dtheta_dr_interact);
  
  H5Dset_extent(fp.dset_dpdrO, dims);
  file_space = H5Dget_space (fp.dset_dpdrO);
  H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, fp.chunk_dims3, NULL);
  H5Dwrite(fp.dset_dpdrO, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &s.dphi_dr_osc);
  
  H5Dset_extent(fp.dset_dpdrI, dims);
  file_space = H5Dget_space (fp.dset_dpdrI);
  H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, fp.chunk_dims3, NULL);
  H5Dwrite(fp.dset_dpdrI, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &s.dphi_dr_interact);
  
  // 1D stuff
  ndims = 1;
  mem_space = H5Screate_simple(ndims, fp.chunk_dims, NULL);
  
  // r
  H5Dset_extent(fp.dset_r, dims);
  file_space = H5Dget_space(fp.dset_r);
  H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, fp.chunk_dims, NULL);
  H5Dwrite(fp.dset_r, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &s.r);

  // rho
  H5Dset_extent(fp.dset_rho, dims);
  file_space = H5Dget_space(fp.dset_rho);
  H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, fp.chunk_dims, NULL);
  H5Dwrite(fp.dset_rho, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &s.rho);

  // Ye
  H5Dset_extent(fp.dset_Ye, dims);
  file_space = H5Dget_space(fp.dset_Ye);
  H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, fp.chunk_dims, NULL);
  H5Dwrite(fp.dset_Ye, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &s.Ye);

  // T
  H5Dset_extent(fp.dset_T, dims);
  file_space = H5Dget_space(fp.dset_T);
  H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, fp.chunk_dims, NULL);
  H5Dwrite(fp.dset_T, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &s.T);

  // free resources
  H5Sclose(file_space);
  H5Sclose(mem_space);
  H5Fflush(fp.file,H5F_SCOPE_LOCAL);
}



#endif
