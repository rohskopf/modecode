/*
This class prints out files that help visualize mode motions.
*/

#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <random>
#include "mpi.h"

#include "visualize.h"
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

Visualize::Visualize(MC *mc) : Ptrs(mc) {
    //fh_debug = fopen("DEBUG_VIS","w");
    //fh_fc2 = fopen("FC2_ASR","w");
    //pr_call=false;
    //esp_call=false;
}

Visualize::~Visualize() 
{

    //fclose(fh_debug);
    
    //fclose(fh_fc2);

    //if (order >= 3) mem->deallocate(fc3);
    //if (order >= 4) mem->deallocate(fc4);

    mem->deallocate(emat);
    mem->deallocate(freq);
    mem->deallocate(xm);
    mem->deallocate(vm);
    mem->deallocate(xm0);
    mem->deallocate(vm0);
    mem->deallocate(gv);

};

/*
Initialize variables.
*/

void Visualize::initialize()
{
    natoms = lmp->atom->natoms;
}

/*
Calculate initial state for a high heat flux superposition state (HHFSPS).
*/

void Visualize::calcInitialState()
{

    printf(" Computing initial state for HHFSPS of mode %d.\n", n_indx); 
    printf(" Temperature = %f K\n", temperature);

    mem->allocate(xm, 3*natoms);
    mem->allocate(vm, 3*natoms);
    mem->allocate(xm0, 3*natoms);
    mem->allocate(vm0, 3*natoms);
    for (int n=0; n<3*natoms; n++){
        xm[n] = 0.0;
        vm[n] = 0.0;
        xm0[n] = 0.0;
        vm0[n] = 0.0;
    }

    // Initialize the state (mode amplitudes and velocities).
    // The mode given by n_indx is initialized with a positive amplitude.
    // The rest of the modes are intitialized with positive velocities.
    xm0[n_indx] = sqrt(kb*temperature)/(2.0*pi*freq[n_indx]*1e12); // Units sqrt(kg)*m
    double qnm;
    for (int n=3; n<3*natoms; n++){
        if (n!=n_indx){
            //printf(" %d %d %e\n", gv[n].n1, gv[n].n2, gv[n].val);
            vm0[n] = sqrt(kb*temperature); // Units sqrt(kg)*m/s
            if (gv[n].val < 0.0) vm0[n] = vm0[n]*(-1.0);
        }
        
        // Calculate qnm just to check it's positive.
        //qnm = gv[n].val*xm0[n_indx]*vm0[n];
        //fprintf(fh_debug, "%d %e %e %e %e\n", n, gv[n].val, xm0[n_indx], vm0[n], qnm);
    }
    
 

}

/*
Calculate time dependence of mode amplitudes and velocities, and print the atomic displacements weighted
by magnitude of GVs.
Realize that the time for a single oscillation to come back is 2*pi/frequency.
So the max time we should use is 2*pi divided by the minimum non-zero frequency.
*/

void Visualize::calcTimeDependence()
{

    printf(" Calculating time dependence of mode amplitudes and velocities.\n");

    // Find maximum time.
    double minfreq = 1.0; // Just choose some number that we know isn't the minimum.
    for (int n=3; n<3*natoms; n++){
        if ( (freq[n] < minfreq) && (freq[n] > 0.0) ){
            minfreq = freq[n];
        }
    }

    printf(" Minimum nonzero frequency: %e THz\n", minfreq);
    double maxtime = (2.0*3.1415926535)/minfreq;
    //double timestep = 10.0; // ps. We don't need a sufficient timestep for MD, this is just visualizing.
    printf(" Need %e ps of time.\n", maxtime);
    int ntimesteps = round(maxtime/timestep);
    ntimesteps = 100;
    printf(" Need %d timesteps.\n", ntimesteps);

    // Now loop through number of timesteps and calculate amplitudes and velocities, and convert to
    // atomic displacements that we can visualize.
    double time; // time in ps
    /*
    for (int t=0; t<ntimesteps; t++){
        time = t*timestep;
        printf("%f\n", time);
    }
    */
    
    // Find largest GV within +/- 10 modes from n_indx.
    largest_gv=abs(gv[n_indx-10].val);
    for (int n=n_indx-9; n<n_indx+10; n++){
        if (abs(gv[n].val) > largest_gv){
            largest_gv = abs(gv[n].val);
        }
    }
    printf(" Largest GV within +/- 10 modes of n=%d: %e\n", n_indx, largest_gv);

    // Calculate time-dependent amplitudes and velocities.
    fh_disp = fopen("disp.xyz", "w");
    mass = lmp->atom->mass;
    type = lmp->atom->type;
    tag = lmp->atom->tag;
    x = lmp->atom->x;
    double qnm;
    double qn;
    double qn_sum;
    for (int t=0; t<ntimesteps; t++){
        qn = 0.0;
        time = t*timestep;
        fprintf(fh_disp, "%d\n", natoms);
        fprintf(fh_disp, "Atoms. Timestep: %d\n", t);
        for (int n=3; n<3*natoms; n++){
            calcAmplitude(n,time);
            calcVelocity(n,time);
            qnm = gv[n].val*xm[n_indx]*vm[n];
            qn += qnm;
        }
        qn_sum += qn/ntimesteps;
        //fprintf(fh_debug, "%f %e\n", time, qn*6.242e+18);
        //fprintf(fh_debug, "%f %e %e %e\n", time, xm[21], xm[22], xm[23]);
        
        // Convert mode amplitudes to atomic displacements.
        convertMode2Cartesian();
        
    }
    printf(" qn_sum: %e\n", qn_sum);
    fclose(fh_disp);
    
}

/*
Calculate mode amplitude based on current time.
Note that frequency is in THz, and time is in ps, so the product is unitless.
*/

void Visualize::calcAmplitude(int n, double time){

    xm[n] = xm0[n]*cos(2.0*pi*freq[n]*time) + (vm0[n]/(2.0*pi*freq[n]*1e12))*sin(2.0*pi*freq[n]*time); // 1e12 because vm0 is sqrt(kg)*m/s and freq is THz.

}

/*
Calculate mode velocity based on current time.
Note that frequency is in THz, and time is in ps, so the product is unitless.
*/

void Visualize::calcVelocity(int n, double time){

    vm[n] = vm0[n]*cos(2.0*pi*freq[n]*time) - xm0[n]*2.0*pi*freq[n]*1e12*sin(2.0*pi*freq[n]*time); // 1e12 because vm0 is sqrt(kg)*m/s and freq is THz.

}

/*
Convert mode amplitudes to atomic displacements.
*/

void Visualize::convertMode2Cartesian(){

    //printf("%e %e\n", abs(gv[21+1].val), largest_gv);

    double ux,uy,uz;
    for (int i=0; i<natoms; i++){
        double mi = mass[type[i]]*(1.0/(6.02214076e23*1e3)); // atom mass in kg
        //printf("%e\n", mi);
        ux = 0.0;
        uy = 0.0;
        uz = 0.0;
        for (int n=3; n<3*natoms; n++){
        
            // Scale displacement by magnitude of GV.
            if (n==n_indx){
                ux += (1.0/sqrt(mi))*emat[3*(tag[i]-1)+0][n]*xm[n]*largest_gv;
                uy += (1.0/sqrt(mi))*emat[3*(tag[i]-1)+1][n]*xm[n]*largest_gv;
                uz += (1.0/sqrt(mi))*emat[3*(tag[i]-1)+2][n]*xm[n]*largest_gv;
            }
        
            else{
                ux += (1.0/sqrt(mi))*emat[3*(tag[i]-1)+0][n]*xm[n]*abs(gv[n].val);
                uy += (1.0/sqrt(mi))*emat[3*(tag[i]-1)+1][n]*xm[n]*abs(gv[n].val);
                uz += (1.0/sqrt(mi))*emat[3*(tag[i]-1)+2][n]*xm[n]*abs(gv[n].val);
            }
            //fprintf(fh_debug, "%d %d %e %e %e\n", i, n, emat[3*(tag[i]-1)+0][n], emat[3*(tag[i]-1)+1][n], emat[3*(tag[i]-1)+2][n]);
            //fprintf(fh_debug, "%d %d %e\n", i, n, xm[n]);

            
        }       
        // Convert to angstrom and scale to better visualize
        ux = ux*1e10*scale_factor;
        uy = uy*1e10*scale_factor;
        uz = uz*1e10*scale_factor;
        //fprintf(fh_debug, "%d %e %e %e\n", i, ux, uy, uz);
        fprintf(fh_disp, "%d %e %e %e\n", type[i], x[i][0]+ux, x[i][1]+uy, x[i][2]+uz);
  }

}

/*
Read eigenvector matrix. 
To index "a" Cartesian component of eigenvector on atom "i" in mode "n1", do:
    emat[a+i*3][n1]
*/

void Visualize::readEmat()
{
    printf(" Reading eigenvector matrix for %d atoms.\n",natoms);

    ifstream readfile;

    mem->allocate(emat,natoms*3,natoms*3);

    for (int i = 0; i < natoms*3; i++) {
        for (int j = 0; j < natoms*3; j++) {
            emat[i][j] = 0.0;
        }
    }

    readfile.open("../EMAT");
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

    /* Read frequencies */
    printf(" Reading frequencies.\n");
    mem->allocate(freq,3*natoms);
    ifstream readfile3;

    for (int i = 0; i < natoms*3; i++) {
        freq[i]=0.0;
    }

    readfile3.open("FREQUENCIES");
    //readfile2.open("ev_real.txt");

    if (!readfile3.is_open()) {
        printf("Unable to open FREQUENCIES.\n");
        exit(1);
    }

    //printf("natoms: %d\n",  natoms);
    for (int i=0;i<3*natoms;i++){
      readfile3>>freq[i];
    }
}

/*
Read GVs.
*/
void Visualize::readGV()
{

    printf(" Reading GVs.\n");

    // Read number of GVs
    int n1,n2;
    double val;
    ifstream fh("../GV");
    string line;
    ngv = 0;
    while (getline(fh, line))
    {
        stringstream ss(line);
        ss >> n1 >> n2 >> val;
        if (n1==n_indx) ngv++;
    }
    //ngv = ngv-1; // subtract 1 for first line
    //printf(" Found %d FC3s.\n", nfc3);
    fh.close();

    printf(" Found %d GVs for mode %d.\n", ngv, n_indx);

    mem->allocate(gv,ngv);

    // Store GVs
    int counter = 0;
    fh.open("../GV");
    getline(fh,line);
    while (getline(fh, line))
    {

        stringstream ss(line);
        ss >> n1 >> n2 >> val;
        if (n1==n_indx){
            gv[counter].n1=n1;
            gv[counter].n2=n2;
            gv[counter].val=val;
            //if (counter==0) printf(" %f\n", val);
            
            //printf("%f\n", fc3[counter].val);
            //printf(" %d\n", counter);
            counter++;

        }
        
    }

    fh.close();

}
