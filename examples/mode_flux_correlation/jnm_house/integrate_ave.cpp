#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <math.h>

using namespace std;

int main(int narg, char **arg){
  int i,j,t,Natoms,L,Nmodes,window;
  double AI,kB,temperature,evps2Js_pwr2,temp,**cor_data_ave,**Area_Q_cor,dt,prefactor,sumtemp,**run_ave;

  char cwd[1024];
  char filename_in[64];
  char filename_out[64];
  char address_in[256];
  char address_out[256];

  double eV2J = 1.60217646e-19;
  double A2m = 1e-10;
  double ps2s = 1e-12;

  if (getcwd(cwd, sizeof(cwd)) != NULL)

  AI = 19.23008919*16.65374565*pow(10,-20);
  evps2Js_pwr2 = ((eV2J/A2m)*A2m/ps2s)*((eV2J/A2m)*A2m/ps2s);	// eV2J/A2m : -dU/dr = Force & A2m/ps2s : velocity
  kB = 1.3806 * pow(10,-23);					// J/K
  temperature = 300;						// K

  Nmodes = 1+1;

  dt = 5.0*pow(10,-15);     //s timestep between each datapoint in Q.txt
  L = 100001; // number of lines in cordata.txt
  window = 10 * 1000; // Remeber: window*dt should be 50ps

  prefactor = ( evps2Js_pwr2*dt*1./(kB*temperature*temperature*AI) )/1e6; // I divided by 1e6 to get units MW/m^2-K
  prefactor = 1.0; // for mode flux correlation.

  cor_data_ave = new double* [L];
  for (i=0;i<L;i++) cor_data_ave[i] = new double [Nmodes];

  Area_Q_cor = new double* [L];
  for (i=0;i<L;i++) Area_Q_cor[i] = new double [Nmodes];

  run_ave = new double* [L-window];
  for (i=0;i<L-window;i++) run_ave[i] = new double [Nmodes];

  for (t=0;t<L;t++) {
    for (i=0;i<Nmodes;i++){
      cor_data_ave[t][i] = 0.0;
      Area_Q_cor[t][i] = 0.0;
    }
  }

  for (t=0;t<L-window;t++)
    for (i=0;i<Nmodes;i++)
      run_ave[t][i] = 0.0;

  int ensemble_number = atoi(arg[1]); // which file to open

  sprintf (filename_in, "cordata%d.txt", ensemble_number);
  sprintf (address_in, "%s/%s",cwd,filename_in);
  ifstream readfile(address_in);
  readfile.precision(20);

  for (t=0;t<L;t++)
    for (j=0;j<Nmodes;j++)
      readfile>>cor_data_ave[t][j];

  // Area_calculation (integration of cordata)
  for (j=0;j<Nmodes;j++){
    sumtemp = 0;
    for (t=0;t<L-1;t++){
      sumtemp += prefactor*(cor_data_ave[t][j] + cor_data_ave[t+1][j])/2.;
      Area_Q_cor[t][j] = sumtemp;
    }
  }

  // Run ave calculation
  for (j=0;j<Nmodes;j++){
    sumtemp = 0;
    for (t=0;t<L-1;t++){
      sumtemp += Area_Q_cor[t][j];
      if (t >= window){
        run_ave[t-window][j] = sumtemp/window;
        sumtemp -= Area_Q_cor[t-window][j];
      }
    }
  }

  // Writing everything in files
  sprintf (filename_out, "cordata_CCs%d.txt", ensemble_number);
  sprintf (address_out, "%s/%s",cwd,filename_out);
  ofstream writefile_cor_ave(address_out);

  sprintf (filename_out, "area_evol_CCs%d.txt", ensemble_number);
  sprintf (address_out, "%s/%s",cwd,filename_out);
  ofstream writefile_area_ave(address_out);

  sprintf (filename_out, "running_ave%d.txt", ensemble_number);
  sprintf (address_out, "%s/%s",cwd,filename_out);
  ofstream writefile_running_ave(address_out);

  for (t=0;t<L-1;t++){
    for (j=0;j<Nmodes;j++){
      // cor
      writefile_cor_ave<<t*dt*pow(10,12);
      for (j=0;j<Nmodes;j++)
        writefile_cor_ave<<" "<<(cor_data_ave[t][j] + cor_data_ave[t+1][j])/2.;
      writefile_cor_ave<<endl;
 
      // area_ave
      writefile_area_ave<<t*dt*pow(10,12);
      for (j=0;j<Nmodes;j++)
        writefile_area_ave<<" "<<Area_Q_cor[t][j];
      writefile_area_ave<<endl;

      // write running_ave
      if (t >= window){
        writefile_running_ave<<t*dt*pow(10,12);
        for (j=0;j<Nmodes;j++)
          writefile_running_ave<<" "<<run_ave[t-window][j];
        writefile_running_ave<<endl;
      }
    }
  }

  writefile_cor_ave.close();
  writefile_area_ave.close();
  writefile_running_ave.close();
  return 78;
}

