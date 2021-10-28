#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <math.h>
#include <fftw3.h>

using namespace std;

int main(int narg, char **arg)
{
  int N,t,Noutput,Neq,L,n,Nmodes;
  double sum,aveHeat,**corArray,temp;
  double *in_Q1,*out_cor1,**in_Q;
  fftw_complex *out_Q1,*sp11;
  fftw_plan my_plan1,my_plan_i1;

  char cwd[1024];
  char filename_in[64];
  char filename_out[64];
  char address_in[256];
  char address_out[256];

  if (getcwd(cwd, sizeof(cwd)) != NULL)

  Neq = 0;
  L = 400000;
  N = L-Neq;
  Nmodes = 5+1;
  Noutput = 100001;

  corArray = (double**) malloc(Noutput*sizeof(double*));
  for (int t = 0; t < Noutput; t++)
     corArray[t] = (double*) malloc(Nmodes*sizeof(double));

  in_Q = (double**) malloc(N*sizeof(double*));
  for (int t = 0; t < N; t++)
     in_Q[t] = (double*) malloc(Nmodes*sizeof(double));

  in_Q1 = (double*) malloc(N*sizeof(double));
  out_cor1 = (double*) malloc(N*sizeof(double));

  out_Q1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  sp11 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);

  my_plan1 = fftw_plan_dft_r2c_1d(N, in_Q1, out_Q1, FFTW_ESTIMATE);
  my_plan_i1 = fftw_plan_dft_c2r_1d(N, sp11, out_cor1, FFTW_ESTIMATE);

  int ensemble_number = atoi(arg[1]); // which file to open
  //char infile[64];
  //sprintf (infile, "jnm%d.dat", ensemble_number;

  sprintf (filename_in, "jnm%d.dat", ensemble_number);
  sprintf (address_in, "%s/%s",cwd,filename_in);
  ifstream readfile(address_in);
  readfile.precision(20);
  
  for (t=0;t<Neq;t++){
    for (n=0;n<Nmodes;n++){
      readfile>>temp;
    }
  }

  for (t=0;t<N;t++){
    for (n=0;n<Nmodes;n++){
      readfile>>in_Q[t][n];
      //if (t<10000000){
      //  printf("%e ", in_Q[t][n]);
      //}
    }
    //if (t<10000000) printf("\n");
  }

  for (n=0;n<Nmodes;n++){ // Loop over all the modes
    for (t=0;t<N;t++) in_Q1[t] = in_Q[t][n];

    sum = 0; // subtract the nonzero mean
    for (t = 0 ; t < N ; t++) sum += in_Q1[t];
    aveHeat = sum/(N*1.);
    for (t = 0 ; t < N ; t++) in_Q1[t] -= aveHeat;

    fftw_execute(my_plan1); 

    for (t=0;t<N;t++){
      sp11[t][0] = out_Q1[t][0]*out_Q1[t][0] + out_Q1[t][1]*out_Q1[t][1];
      sp11[t][1] = -out_Q1[t][0]*out_Q1[t][1] + out_Q1[t][1]*out_Q1[t][0];
    }  
    
    fftw_execute(my_plan_i1); 

    for (t=0;t<Noutput;t++){  
      corArray[t][n] = out_cor1[t];
    }
  } // end loop over all modes

  // Write the correlation file
  sprintf (filename_out, "cordata%d.txt", ensemble_number);
  sprintf (address_out, "%s/%s",cwd,filename_out);
  ofstream outfile(address_out);
  outfile.precision(20);
  
  
  for (t=0;t<Noutput;t++){
    for (n=0;n<Nmodes;n++){
      //printf("%d\n", n);
      outfile<<corArray[t][n]/(1.*N*N)<<" ";}
    outfile<<endl;
  }

  // Close everything
  readfile.close();
  fftw_destroy_plan(my_plan1);
  fftw_destroy_plan(my_plan_i1);
  fftw_free(in_Q1);
  fftw_free(out_Q1);
  fftw_free(sp11);
  fftw_free(out_cor1);
  outfile.close();
  return 1;
}
