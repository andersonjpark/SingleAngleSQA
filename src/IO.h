#ifndef OUTPUT_H
#define OUTPUT_H
#include <string>
#include "H5Cpp.h"
#include "State.h"
#include <iomanip>
using namespace std;


class FilePointers{
 public:
  hid_t file;
  hid_t dset_f, dset_S, dset_U, dset_r, dset_dr_osc, dset_dr_int, dset_dr_block,
    dset_rho, dset_Ye, dset_T, dset_VfSI, dset_dtdrO, dset_dtdrI, dset_dpdrO, dset_dpdrI, dset_Elab_Elabstart, dset_Ecom_Elab;
  const hsize_t        dims0[6] = {0,             NM, NE, NF, NF, 2};
  const hsize_t     max_dims[6] = {H5S_UNLIMITED, NM, NE, NF, NF, 2};
  const hsize_t   chunk_dims[6] = {1,             NM, NE, NF, NF, 2};
  const hsize_t      dims0_V[5] = {0,             NM, NF, NF, 2};
  const hsize_t   max_dims_V[5] = {H5S_UNLIMITED, NM, NF, NF, 2};
  const hsize_t chunk_dims_V[5] = {1,             NM, NF, NF, 2};
  const hsize_t       dims03[3] = {0,             NM, NE};
  const hsize_t    max_dims3[3] = {H5S_UNLIMITED, NM, NE};
  const hsize_t  chunk_dims3[3] = {1,             NM, NE};
};

//============//
// setup_file //
//============//
FilePointers setup_HDF5_file(const string& outputfilename, const array<double,NE>& E, const array<double,NE>& Etop){
  FilePointers fp;
  
  fp.file = H5Fcreate(outputfilename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  hid_t file_space;
  hsize_t ndims;
  hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_layout(plist, H5D_CHUNKED);

  // FMATRIXF //
  ndims = 6;
  file_space = H5Screate_simple(ndims, fp.dims0, fp.max_dims);
  H5Pset_chunk(plist, ndims, fp.chunk_dims);
  H5Dcreate(fp.file, "fmatrixf", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Dcreate(fp.file, "S",        H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Dcreate(fp.file, "U",        H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  fp.dset_f = H5Dopen(fp.file, "fmatrixf", H5P_DEFAULT);
  fp.dset_S = H5Dopen(fp.file, "S", H5P_DEFAULT);
  fp.dset_U = H5Dopen(fp.file, "U", H5P_DEFAULT);

  // VfSI //
  ndims = 5;
  file_space = H5Screate_simple(ndims, fp.dims0_V, fp.max_dims_V);
  H5Pset_chunk(plist, ndims, fp.chunk_dims_V);
  H5Dcreate(fp.file, "VfSI(erg)",H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  fp.dset_VfSI = H5Dopen(fp.file, "VfSI(erg)", H5P_DEFAULT);

  // dangle_dr //
  ndims = 3;
  file_space = H5Screate_simple(ndims, fp.dims03, fp.max_dims3);
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
  file_space = H5Screate_simple(ndims, fp.dims0, fp.max_dims);
  H5Pset_chunk(plist, ndims, fp.chunk_dims);
  H5Dcreate(fp.file, "r(cm)", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Dcreate(fp.file, "dr_block(cm)", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Dcreate(fp.file, "dr_osc(cm)", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Dcreate(fp.file, "dr_int(cm)", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Dcreate(fp.file, "rho(g|ccm)", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Dcreate(fp.file, "Ye", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Dcreate(fp.file, "T(MeV)", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Dcreate(fp.file, "Elab_Elabstart", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Dcreate(fp.file, "Ecom_Elab", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  fp.dset_r        = H5Dopen(fp.file, "r(cm)", H5P_DEFAULT);
  fp.dset_dr_osc   = H5Dopen(fp.file, "dr_osc(cm)", H5P_DEFAULT);
  fp.dset_dr_int   = H5Dopen(fp.file, "dr_int(cm)", H5P_DEFAULT);
  fp.dset_dr_block = H5Dopen(fp.file, "dr_block(cm)", H5P_DEFAULT);
  fp.dset_rho      = H5Dopen(fp.file, "rho(g|ccm)", H5P_DEFAULT);
  fp.dset_Ye       = H5Dopen(fp.file, "Ye", H5P_DEFAULT);
  fp.dset_T        = H5Dopen(fp.file, "T(MeV)", H5P_DEFAULT);
  fp.dset_Elab_Elabstart = H5Dopen(fp.file, "Elab_Elabstart", H5P_DEFAULT);
  fp.dset_Ecom_Elab = H5Dopen(fp.file, "Ecom_Elab", H5P_DEFAULT);
  
  
  // energy grid //
  hid_t mem_space =   H5Screate_simple(ndims, &fp.dims0[2], NULL);
  file_space = H5Screate_simple(ndims, &fp.dims0[2], &fp.dims0[2]);
  hid_t dset_E0 = H5Dcreate(fp.file, "E0(erg)", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Dwrite(dset_E0, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &E[0]);
  hid_t dset_Etop0 = H5Dcreate(fp.file, "Etop0(erg)", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Dwrite(dset_Etop0, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &Etop[0]);

  // clear resources
  H5Sclose(file_space);
  H5Pclose(plist);
  H5Fflush(fp.file,H5F_SCOPE_LOCAL);

  return fp;
}

//============//
// write_data //
//============//
void write_data_HDF5(FilePointers& fp, const State& s, double dr_osc, double dr_int, double dr_block, bool do_header){
  // output to stdout
  int colwidth = 11;
  int precision = 4;

  // integrate neutrino distributions
  array<array<double,NF>,NM> ndens, ndensE;
  for(int i=0; i<NE; i++){
    double vol = Vphase(i,s.Etop);
    for(state m=matter; m<=antimatter; m++){
      for(flavour f=e; f<=mu; f++){
	ndens[m][f]  += real(s.fmatrixf[m][i][f][f]) * vol;
	ndensE[m][f] += real(s.fmatrixf[m][i][f][f]) * vol * s.E[i]/(1e6*cgs::units::eV);
      }
    }
  }
  
  // preliminary stdout
  if(do_header){
    cout << setw(colwidth) << "r(km)";
    cout << setw(colwidth) << "rho(g/ccm)";
    cout << setw(colwidth) << "Ye";
    cout << setw(colwidth) << "T(MeV)";
    //cout << setw(colwidth) << "eta(MeV)";
    //cout << setw(colwidth) << "munue_kT";
    cout << setw(colwidth) << "grv_rdshft";
    cout << setw(colwidth) << "vel_rdshft";
    cout << setw(colwidth) << "dr_osc(cm)";
    cout << setw(colwidth) << "dr_int(cm)";
    cout << setw(colwidth) << "dr_blk(cm)";
    cout << " |";
    for(int m=matter; m<=antimatter; m++){
      for(int f=e; f<=mu; f++){
	string name = string("m")+to_string(m)+string("f")+to_string(f)+string("ndens");
	cout << setw(colwidth) << name;
      }
    }
    for(int m=matter; m<=antimatter; m++){
      for(int f=e; f<=mu; f++){
	string name = string("m")+to_string(m)+string("f")+to_string(f)+string("avgE");
	cout << setw(colwidth) << name;
      }
    }
    cout << endl;
  }
  
  cout << setw(colwidth) << scientific << setprecision(precision) << s.r/1e5;
  cout << setw(colwidth) << scientific << setprecision(precision) << s.rho;
  cout << setw(colwidth) << scientific << setprecision(precision) << s.Ye;
  cout << setw(colwidth) << scientific << setprecision(precision) << s.T;
  //cout << setw(colwidth) << scientific << setprecision(precision) << eas.eta;
  //cout << setw(colwidth) << scientific << setprecision(precision) << eas.munue_kT;
  cout << setw(colwidth) << scientific << setprecision(precision) << s.Elab_Elabstart;
  cout << setw(colwidth) << scientific << setprecision(precision) << s.Ecom_Elab;
  cout << setw(colwidth) << scientific << setprecision(precision) << dr_osc;
  cout << setw(colwidth) << scientific << setprecision(precision) << dr_int;
  cout << setw(colwidth) << scientific << setprecision(precision) << dr_block;
  /* cout << impact << endl; */
  cout << " |";
  for(int m=matter; m<=antimatter; m++){
    for(int f=e; f<=mu; f++){
      cout << setw(colwidth) << scientific << setprecision(precision) << ndens[m][f];
    }
  }
  for(int m=matter; m<=antimatter; m++){
    for(int f=e; f<=mu; f++){
      cout << setw(colwidth) << scientific << setprecision(precision) << ndensE[m][f]/ndens[m][f];
    }
  }
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

  // Elab_Elabstart
  H5Dset_extent(fp.dset_Elab_Elabstart, dims);
  file_space = H5Dget_space(fp.dset_Elab_Elabstart);
  H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, fp.chunk_dims, NULL);
  H5Dwrite(fp.dset_Elab_Elabstart, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &s.Elab_Elabstart);

  // Ecom_Elab
  H5Dset_extent(fp.dset_Ecom_Elab, dims);
  file_space = H5Dget_space(fp.dset_Ecom_Elab);
  H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, fp.chunk_dims, NULL);
  H5Dwrite(fp.dset_Ecom_Elab, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &s.Ecom_Elab);

  // free resources
  H5Sclose(file_space);
  H5Sclose(mem_space);
  H5Fflush(fp.file,H5F_SCOPE_LOCAL);
}



#endif
