#ifndef OUTPUT_H
#define OUTPUT_H
#include <string>
#include "hdf5.h"
#include "State.h"
using namespace std;

void load_input_data(string input_directory,
		     DISCONTINUOUS& lnrho,
		     DISCONTINUOUS& Ye,
		     DISCONTINUOUS& temperature,
		     array<array<array<DISCONTINUOUS,NF>,NE>,NM>& P_unosc,
		     array<array<array<DISCONTINUOUS,NF>,NE>,NM>& D_unosc){

  // load rho and Ye data
  lnrho.Open(input_directory+"/rho.txt",'#');
  Ye.Open(input_directory+"/Ye.txt",'#');
  temperature.Open(input_directory+"/temp.txt",'#');
  lnrho = lnrho.copy_logy();
    
  // load and compute spectral data
  for(int i=0; i<NE; i++){
    P_unosc[matter    ][i][e ].Open(input_directory+"/potential_s1_g"+to_string(i+1)+"_.txt",'#');
    P_unosc[antimatter][i][e ].Open(input_directory+"/potential_s2_g"+to_string(i+1)+"_.txt",'#');
    P_unosc[matter    ][i][mu].Open(input_directory+"/potential_s3_g"+to_string(i+1)+"_.txt",'#');
    P_unosc[antimatter][i][mu].Open(input_directory+"/potential_s3_g"+to_string(i+1)+"_.txt",'#');
    D_unosc[matter    ][i][e ].Open(input_directory+"/density_s1_g"+to_string(i+1)+"_.txt",'#');
    D_unosc[antimatter][i][e ].Open(input_directory+"/density_s2_g"+to_string(i+1)+"_.txt",'#');
    D_unosc[matter    ][i][mu].Open(input_directory+"/density_s3_g"+to_string(i+1)+"_.txt",'#');
    D_unosc[antimatter][i][mu].Open(input_directory+"/density_s3_g"+to_string(i+1)+"_.txt",'#');
  }

}

void init_output(string outputfilename,
		 ofstream &fout,
		 ofstream &foutP,
		 ofstream &foutf,
		 ofstream &foutdangledr){
    foutdangledr.open((outputfilename+"/dangledr.dat").c_str());
    foutdangledr.precision(12);
    foutdangledr << "# 1:r ";
    for(state m=matter; m<=antimatter; m++)
      for(int i=0; i<NE; i++)
	foutdangledr << 2+0*NE*2 + i + m*NE << ":OscThetaie"<<i<<"m"<<m<<"\t";
    for(state m=matter; m<=antimatter; m++)
      for(int i=0; i<NE; i++)
	foutdangledr << 2+1*NE*2 + i + m*NE << ":OscPhiie"<<i<<"m"<<m<<"\t";
    for(state m=matter; m<=antimatter; m++)
      for(int i=0; i<NE; i++)
	foutdangledr << 2+2*NE*2 + i + m*NE << ":InteractThetaie"<<i<<"m"<<m<<"\t";
    for(state m=matter; m<=antimatter; m++)
      for(int i=0; i<NE; i++)
	foutdangledr << 2+3*NE*2 + i + m*NE << ":InteractPhiie"<<i<<"m"<<m<<"\t";
    foutdangledr << endl;
    foutdangledr.flush();

}


void Outputvsr(ofstream &foutdangledr,
	       const State& s,
	       const array<array<array<DISCONTINUOUS,NF>,NE>,NM>& P_unosc){
  foutdangledr << s.r << "\t";
  for(state m=matter; m<=antimatter; m++)
    for(int i=0; i<NE; i++)
      foutdangledr << s.dtheta_dr_osc[i][m] << "\t";
  for(state m=matter; m<=antimatter; m++)
    for(int i=0; i<NE; i++)
      foutdangledr << s.dphi_dr_osc[i][m] << "\t";
  for(state m=matter; m<=antimatter; m++)
    for(int i=0; i<NE; i++)
      foutdangledr << s.dtheta_dr_interact[i][m] << "\t";
  for(state m=matter; m<=antimatter; m++)
    for(int i=0; i<NE; i++)
      foutdangledr << s.dphi_dr_interact[i][m] << "\t";
  foutdangledr << endl;
  foutdangledr.flush();
}

class FilePointers{
 public:
  hid_t file;
  hid_t dset_f, dset_S, dset_U, dset_r, dset_dr_osc, dset_dr_int, dset_dr_block,
    dset_rho, dset_Ye, dset_T, dset_VfSI;
  const hsize_t       dims[6] = {0,             NM, NE, NF, NF, 2};
  const hsize_t   max_dims[6] = {H5S_UNLIMITED, NM, NE, NF, NF, 2};
  const hsize_t chunk_dims[6] = {1,             NM, NE, NF, NF, 2};
  const hsize_t     dims_V[5] = {0,             NM, NF, NF, 2};
  const hsize_t max_dims_V[5] = {H5S_UNLIMITED, NM, NF, NF, 2};
  const hsize_t chunk_dims_V[5] = {1,           NM, NF, NF, 2};
};

//============//
// setup_file //
//============//
FilePointers setup_HDF5_file(const array<double,NE>& E){
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
  for(int i=0; i<NE; i++) cout << E[i] << " ";
  cout << endl;

  // clear resources
  H5Sclose(file_space);
  H5Pclose(plist);
  H5Fflush(fp.file,H5F_SCOPE_LOCAL);
  
  return fp;
}

//============//
// write_data //
//============//
void write_data_HDF5(FilePointers& fp, const State& s){
  // output to stdout
  double n=0, nbar=0;
  double coeff = 4.*M_PI / pow(cgs::constants::c,3);
  cout << s.counter << "\t";
  cout << s.r/1e5 << "\t";
  /* cout << s.dr_osc/cgs::constants::c << "\t"; */
  /* cout << s.dr_int/cgs::constants::c << "\t"; */
  /* cout << s.dr_block/cgs::constants::c << "\t"; */
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
