#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <math.h>

using namespace std;

int main(){
  int i,j,t,Natoms,L,Nmodes;
  double AI,kB,temperature,evps2Js_pwr2,temp,**cor_data_ave,*Area_Q_cor,dt,prefactor;

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

  prefactor = ( evps2Js_pwr2*dt*1./(kB*temperature*temperature*AI) )/1e6; // I divided by 1e6 to get units MW/m^2-K
  prefactor = 1.0; // for mode flux correlation.

  cor_data_ave = new double* [L];
  for (i=0;i<L;i++) cor_data_ave[i] = new double [Nmodes];

  Area_Q_cor = new double [Nmodes];

 for (t=0;t<L;t++) {
   for (i=0;i<Nmodes;i++){
      cor_data_ave[t][i] = 0.0;
    }
  }

  sprintf (filename_in, "cordata.txt");
  sprintf (address_in, "%s/%s",cwd,filename_in);
  ifstream readfile(address_in);
  readfile.precision(20);

  for (t=0;t<L;t++){
    for (j=0;j<Nmodes;j++){
      readfile>>cor_data_ave[t][j];
    }
  }

  sprintf (address_out, "%s/%s",cwd,"cordata_CCs.txt");         //////////////////////
  ofstream writefile_cor_ave(address_out);         //////////////////////

  sprintf (address_out, "%s/%s",cwd,"area_evol_CCs.txt");         //////////////////////
  ofstream writefile_area_ave(address_out);         //////////////////////

  for (j=0;j<Nmodes;j++){
    Area_Q_cor[j] = 0.0;
  }

  //for (t=0;t<L;t++){
  for (t=0;t<L-1;t++){
    for (j=0;j<Nmodes;j++){ 
      //Area_Q_cor[j] += cor_data_ave[t][j];
      //if (t==0) Area_Q_cor[j] += cor_data_ave[t][j];
      //else Area_Q_cor[j] += (cor_data_ave[t-1][j] + cor_data_ave[t][j])/2.;
      //if (t==(L-1)) Area_Q_cor[j] += cor_data_ave[t][j];
      //else Area_Q_cor[j] += (cor_data_ave[t][j] + cor_data_ave[t+1][j])/2.;
      Area_Q_cor[j] += (cor_data_ave[t][j] + cor_data_ave[t+1][j])/2.;
    }

    writefile_cor_ave<<t*dt*pow(10,12);
    for (j=0;j<Nmodes;j++){
      writefile_cor_ave<<" "<<(cor_data_ave[t][j] + cor_data_ave[t+1][j])/2.;
    }
    writefile_cor_ave<<endl;

    writefile_area_ave<<t*dt*pow(10,12);
    for (j=0;j<Nmodes;j++){
      writefile_area_ave<<" "<<prefactor*Area_Q_cor[j];
    }
    writefile_area_ave<<endl;
  }

  writefile_cor_ave.close();
  writefile_area_ave.close();
  return 78;
}

