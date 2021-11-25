/*
This class post processes simulation data to analyze the modes.
Current functions include:

*/

#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <random>
#include "mpi.h"
#include <fftw3.h>

#include "postproc.h"
#include "mem.h"
#include "in.h"
//#include "config.h"

// LAMMPS include files
#include "lammps.h"
#include "input.h"
#include "atom.h"
#include "library.h"
#include "memory.h"

#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "pair.h"
#include "domain.h"

#include <ctime>

using namespace std;

using namespace MC_NS;

Postproc::Postproc(MC *mc) : Ptrs(mc) {
    //fh_debug = fopen("DEBUG_PP","w");
    //fh_fc2 = fopen("FC2_ASR","w");
    
    rank = mc->rank;

}

Postproc::~Postproc() 
{

    //fclose(fh_debug);
    //fclose(fh_fc2);

    //if (order >= 3) mem->deallocate(fc3);
    //if (order >= 4) mem->deallocate(fc4);
    //mem->deallocate(emat);
    
  if (indices_bool){
    mem->deallocate(indices);
  }

};

/*
Initialize
*/

void Postproc::initialize()
{

  //if (rank==0) printf(" Initializing post processing module.\n");
  natoms = lmp->atom->natoms;
  nmodes = 3*natoms;
  if (task==1){
  
    if (rank==0) printf("   Task 1: Calculating DOS overlap.\n");
    // Read INDICES file
    ifstream fh_ind;
    fh_ind.open("INDICES");

    if (fh_ind.is_open()) {
        //printf("Unable to open INDICES file, won't output any mode quantities.\n");
      indices_bool = true;
      string line;
      nind = 0;
      while (getline(fh_ind, line))
      {
          nind=nind+1;
      }

      if (rank==0) printf(" Found %d mode indices in INDICES.\n", nind);

      mem->allocate(indices,nind);

      fh_ind.close();

      ifstream fh_ind;
      fh_ind.open("INDICES");

      for (int n=0; n<nind; n++){
        fh_ind >> indices[n];
      }

      fh_ind.close();
     
    }
    else {
      indices_bool = false;
      if (rank==0) printf("Unable to open INDICES file, won't output any mode quantities.\n");
      nind = 0.0;
    }    
  }
}

void Postproc::task1()
{

  if (rank==0) printf(" Settings:\n");
  if (rank==0) printf("   Ensemble dirname: %s\n", ensemble_dirname.c_str());
  if (rank==0) printf("   ntimesteps: %d\n", ntimesteps);
  if (rank==0) printf("   Data sampled every %f ps\n", sampling_interval);
  if (rank==0) printf("   Number of ensembles to average PS overlap for: %d\n", nens);
  
  //Post process the input settings to determine frequency axis of power spectrums.
  double sampling_frequency = 1.0/sampling_interval; // frequency of data collection
  double end_time = ntimesteps*sampling_interval;
  if (rank==0) printf(" End time: %e\n", end_time);
  
  // There will be (N/2)+1 output values from the FFT (real2complex does it symmetrically), where N = ntimesteps
  // So each index from 0 to (N/2)+1 represents a frequency. 
  // Need to divide these indices by the number of timesteps. 
  int nfreq = (ntimesteps/2)+1;
  if (rank==0) printf(" Number of freqency points: %d\n", nfreq);
  double *freq;
  mem->allocate(freq, (ntimesteps/2)+1);
  for (int f=0; f<(ntimesteps/2)+1; f++){
    freq[f] = f/(end_time);
  }
  // The maximum frequency will be ( (ntimesteps/2)/end_time )
  if (rank==0) printf(" Max frequency we can sample: %f\n", (ntimesteps/2)/end_time);
  
  /*--------------------------------------------------------------
              Begin declaring and allocating variables and arrays
    -------------------------------------------------------------*/
  if (rank==0) printf(" Declaring/allocating variables and arrays.\n");
  double **xm;
  double **vm;
  double *time;
  //double *ps1; // power spectrum of mode amplitude (xn in Qnm = Xn*Vm)
  //double *ps2; // power spectrum of mode velocity (vm in Qnm = Xn*Vm)
  mem->allocate(xm,ntimesteps,nind);
  mem->allocate(vm,ntimesteps,nind);
  mem->allocate(time,ntimesteps);
  //mem->allocate(ps1, ntimesteps);
  //mem->allocate(ps2, ntimesteps);
  // Allocate arrays for storing amplitude PS (xm_ps) and velocity PS (vm_ps)
  double **xm_ps;
  double **vm_ps;
  mem->allocate(xm_ps, nfreq,nind);
  mem->allocate(vm_ps, nfreq,nind);
  for (int f=0; f<nfreq; f++){
    for (int n=0; n<nind; n++){
      xm_ps[f][n] = 0.0;
      vm_ps[f][n] = 0.0;
    }
  }
  // Allocate array for calculating PS overlaps
  double *product, *amplitude_ps, *velocity_ps, *overlaps_p, *overlaps;
  double numerator, denominator;
  mem->allocate(product, nfreq);
  mem->allocate(amplitude_ps, nfreq);
  mem->allocate(velocity_ps, nfreq);
  mem->allocate(overlaps_p, nind*nind);
  mem->allocate(overlaps, nind*nind);
  for (int n=0; n<nind; n++){
    for (int m=0;m<nind;m++){
      overlaps[n*nind+m]=0.0;
      overlaps_p[n*nind+m]=0.0;
    }
  }
  // Declare the FFT output
  // out[i][0] and out[i][1] are the real and imaginary parts of a complex number.
  fftw_complex *out;
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ntimesteps);
  fftw_complex *out2;
  out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ntimesteps);
  // Declare the FFT input and plan
  double *in;
  mem->allocate(in, ntimesteps);
  double *in2;
  mem->allocate(in2, ntimesteps);
  // plan for Fourier transform, save the result in 'out'
  fftw_plan p; // plan for xn
  fftw_plan p2; // plan for vm
  p = fftw_plan_dft_r2c_1d(ntimesteps, in, out, FFTW_ESTIMATE);  
  p2 = fftw_plan_dft_r2c_1d(ntimesteps, in2, out2, FFTW_ESTIMATE);  
  // Split modes across procs, for parallelizing the PS overlap part.
  int *nepp;
  mem->allocate(nepp, mc->nprocs); // number elements (MCC3s) per proc
  for (int p=0; p<mc->nprocs; p++){
      nepp[p] = nind/mc->nprocs;
  }
  // divide up the remainder
  for (int p=0; p<(nind % mc->nprocs); p++){
      nepp[p] += 1;
  }
  int start_indx = 0;
  for (int p=0; p<rank; p++){
      start_indx += nepp[p];
  }
  int end_indx = 0; //napp[0]-1;
  for (int p=0; p<rank+1; p++){
      end_indx += nepp[p];
  }
  end_indx=end_indx+1-1;
  //end_indx = end_indx+1-1;
  //printf("rank startindx endindx: %d %d %d\n", rank,start_indx, end_indx);
  if (rank==0){
      printf(" Splitting modes on procs like:\n");
      for (int p=0; p<mc->nprocs; p++){
          printf("  %d Modes on proc %d.\n", nepp[p],p);
      }
  }
  /*--------------------------------------------------------------
              Done declaring and allocating variables and arrays
    -------------------------------------------------------------*/
    
  
  int ens = 1; // current ensemble
  char filename_xm[1000];
  char filename_vm[1000];
  
  // Loop over all ensembles.
  
  for (int ens=1; ens<nens+1; ens++){
  
    sprintf(filename_xm, "%s%d/xm.dat", ensemble_dirname.c_str(), ens);
    //printf(" Opening %s\n", filename_xm);
    sprintf(filename_vm, "%s%d/vm.dat", ensemble_dirname.c_str(), ens);
    if (rank==0) printf(" Opening %s and %s\n", filename_xm,filename_vm);
    //std::cout << filename_xm << std::endl;
    
    ifstream fh_xm;
    fh_xm.open(filename_xm);
    ifstream fh_vm;
    fh_vm.open(filename_vm);
   
    //nmodes=20;
    double junk;
    if (fh_xm.is_open() && fh_vm.is_open()) {

      for (int t=0; t<ntimesteps; t++){
        //printf("%d\n", t);
        fh_xm>>time[t];
        fh_vm>>junk;
        for (int n=0; n<nind; n++){
          //printf("  %d\n", n);
          fh_xm>>xm[t][n];
          fh_vm>>vm[t][n];
        }
      }

      fh_xm.close();
      fh_vm.close();
    }
    else {
      if (rank==0) printf("Unable to open %s or %s.\n", filename_xm, filename_vm);
    }
    
    // Check values
    /*
    for (int t=0; t<ntimesteps; t++){
      printf("%f ", time[t]);
      for (int n=0; n<nmodes; n++){
        printf("%f ", xm[t][n]);
      }
      printf("\n");
    }
    */
    
    // Make first time to be zero.
    double start_time = time[0];
    for (int t=0; t<ntimesteps; t++){
      time[t] = time[t] - start_time;
    }
    
    // Subtract the mean from the data.
    for (int n=0; n<nind; n++){
      double mean=0.0;
      double mean2=0.0;
      for (int t=0; t<ntimesteps; t++){
        mean += xm[t][n]/ntimesteps;
        mean2 += vm[t][n]/ntimesteps;
      }
      for (int t=0; t<ntimesteps; t++){
        xm[t][n] = xm[t][n]-mean;
        vm[t][n] = vm[t][n]-mean2;
      }
    }
    
    // Loop over mode amplitudes and velocities, calculate their power spectrum, and store it.
    //int modeindx = 11;
    //for (int n=modeindx; n<modeindx+1; n++){
    if (rank==0) printf(" Calculating power spectrums for ensemble %d.\n", ens);
    for (int n=0; n<nind; n++){
    
      //printf(" n = %d\n", n);
    
      // Build the input array for FFT calculation.
      for (int t=0; t<ntimesteps; t++){
        in[t] = xm[t][n];
        in2[t] = vm[t][n];
        //printf("%f\n", in[t]);
      }
      
      // Calculate FFT
      fftw_execute(p); // amplitude
      fftw_execute(p2); // velocity
      
      // Calculate the power spectrum
      //printf("freq | powerspectrum\n");
      //fh = fopen("ps.dat","w");
      for (int f = 0; f < nfreq; f++){
        //ps1[f] = (out[f][0]*out[f][0] + out[f][1]*out[f][1])/(ntimesteps*ntimesteps);
        //ps1[i] = 1.0;
        //printf("freq: %3d %+9.5f %+9.5f I\n", i, out[i][0], out[i][1]);
        //fprintf(fh, "%e %e\n", freq[f],ps1[f]);
        xm_ps[f][n] = (out[f][0]*out[f][0] + out[f][1]*out[f][1])/(ntimesteps*ntimesteps);
        vm_ps[f][n] = (out2[f][0]*out2[f][0] + out2[f][1]*out2[f][1])/(ntimesteps*ntimesteps);
      }
      //fclose(fh);
      
    } // for int n
    
    // Now loop over all possible pairs of xm[n] and vm[m] to get PS overlap.
    
    // Testing integrate function
    /*
    double *test;
    double *x;
    int length = 10000;
    double finaltime = 10.0;
    double dx = finaltime/length;
    mem->allocate(test,length);
    mem->allocate(x,length);
    for (int i=0; i<length; i++){
      test[i] = (2.0*i*i*dx*dx);
      x[i] = (1.0*i*dx);
    }
    double integral = integrate(test,x, length);
    printf(" integral from 0 to %f: %f\n", finaltime,integral);
    mem->deallocate(test);
    mem->deallocate(x);
    */
    
    // Calculate PS overlaps
    if (rank==0) printf(" Calculating PS overlaps for ensemble %d.\n", ens);
    double integral_n;
    for (int n=start_indx; n<end_indx; n++){
      if (rank==0) printf(" n = %d\n", n);
      for (int f=0; f<nfreq; f++){
        amplitude_ps[f] = xm_ps[f][n];
      }
      integral_n = integrate(amplitude_ps, freq, nfreq);
      for (int m=0; m<nind; m++){
      
        // Calculate product of xm_ps[n] and vm_ps[m]
        for (int f=0; f<nfreq; f++){
          product[f] = xm_ps[f][n]*vm_ps[f][m];
          //amplitude_ps[f] = xm_ps[f][n];
          velocity_ps[f] = vm_ps[f][m];
        }
        // Integrate the product, xm_ps and vm_ps
        numerator = integrate(product, freq, nfreq);
        denominator = integral_n*integrate(velocity_ps, freq, nfreq);
        overlaps_p[n*nind+m] += (numerator/denominator)/nens; // contribute to the ensemble averaged PS overlap
      }
    }
  
  } // for (int ens=1; ens<nens+1; ens++)
  
  // Reduce all overlaps on each proc.
  MPI_Allreduce(overlaps_p,overlaps,nind*nind,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  
  // Write the overlaps to a file
  if (rank==0){
    FILE * fh;
    fh = fopen("overlaps.dat", "w");
    for (int n=0; n<nind; n++){
      for (int m=0; m<nind; m++){
        if (n!=m){
          fprintf(fh, "%d %d %.5e\n", indices[n],indices[m],overlaps[n*nind+m]);
        }
      }
    }
    
    fclose(fh);
  }
  
  /*
  mem->deallocate(amplitude_ps);
  mem->deallocate(velocity_ps);
  mem->deallocate(overlaps);
  mem->deallocate(product);
  */
  
  // SIMPLE EXAMPLE TO SHOW THAT FFT WORKS -------------------------------------
  /*
  int N=ntimesteps;
  double *time_test;
  mem->allocate(time_test,N);
  
  double in[N];
  for (int t=0; t<N; t++){
    time_test[t] = t*sampling_interval;
  }
  
  fftw_plan p, q;
  int i;
  // prepare a cosine wave with two frequencies
  for (i = 0; i < N; i++) {
    //in[i][0] = cos(3 * 2*M_PI*i/N);
    //in[i][1] = 0;
    in[i] = cos(2.1 * 2*M_PI*time_test[i]) + sin(4 * 2*M_PI*time_test[i]);
  }
  
  // forward Fourier transform, save the result in 'out'
  p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);
  //p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  // Calculate the power spectrum
  printf("freq | powerspectrum\n");
  for (i = 0; i < (N/2)+1; i++){
    ps1[i] = out[i][0]*out[i][0] + out[i][1]*out[i][1];
    //ps1[i] = 1.0;
    //printf("freq: %3d %+9.5f %+9.5f I\n", i, out[i][0], out[i][1]);
    printf("%f   %e\n", freq[i],ps1[i]);
  }
  fftw_destroy_plan(p);
  
  //mem->deallocate(ps1);
  mem->deallocate(time_test);
  */
  // END OF SIMPLE EXAMPLE -------------------------------------
  
  
  fftw_free(out);
  fftw_free(out2);
  fftw_cleanup();
  fftw_destroy_plan(p);
  
  mem->deallocate(in);
  mem->deallocate(in2);
  mem->deallocate(xm);
  mem->deallocate(vm);
  mem->deallocate(time);
  mem->deallocate(nepp);
  //mem->deallocate(ps1);
  //mem->deallocate(ps2);
  mem->deallocate(freq);
  mem->deallocate(xm_ps);
  mem->deallocate(vm_ps);
  mem->deallocate(amplitude_ps);
  mem->deallocate(velocity_ps);
  mem->deallocate(overlaps_p);
  mem->deallocate(overlaps);
  mem->deallocate(product);
  

}

double Postproc::integrate(double *function, double *x, int length)
{

  double integral = 0.0;
  for (int i=0; i<length-1; i++){
    //printf("%f %f\n", x[i+1], x[i]);
    integral += (x[i+1]-x[i])*0.5*(function[i+1]+function[i]);
  }
  
  
  return integral;

}

void Postproc::calc1()
{

  if (rank==0) printf(" Calculating autocorrelation of mode action pairs.\n");
  double **test;
  double **test2;
  
  long long int n1 = 400000;
  long long int n2 = 3000;
  
  mem->allocate(test,n1,n2);
  mem->allocate(test2,n1,n2);
  size_t size_of_one = sizeof(test[0]);
  printf(" %ld\n", size_of_one);

  long long int n = n1*n2*size_of_one;
  std::cout << n << std::endl;
  n = n/1e6;
  printf(" MB: %lld\n", n);
  
  mem->deallocate(test);
  mem->deallocate(test2);
}

void Postproc::readEmat()
{

    if (rank==0) printf(" Reading EMAT.\n");
    if (rank==0) printf(" Reading eigenvector matrix for %d atoms.\n",natoms);

    ifstream readfile;

    mem->allocate(emat,natoms*3,natoms*3);

    for (int i = 0; i < natoms*3; i++) {
        for (int j = 0; j < natoms*3; j++) {
            emat[i][j] = 0.0;
        }
    }

    readfile.open("EMAT");
    //readfile2.open("ev_real.txt");
  
    if (!readfile.is_open()) {
        cout<<"Unable to open the file!"<<endl;
        exit(1);
    }

    //printf("natoms: %d\n",  natoms);
    for (int i=0;i<3*natoms;i++){
        for (int j=0;j<3*natoms;j++){
            readfile>>emat[i][j];
        }
    }
}
