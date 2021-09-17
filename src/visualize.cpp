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
    //fh_debug = fopen("DEBUG_ASR","w");
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
    for (int n=0; n<3*natoms; n++){
        xm[n] = 0.0;
        vm[n] = 0.0;
    }

    // Initialize based on temperature.
    xm[n_indx] = sqrt(kb*temperature)/freq[n_indx]; // Units sqrt(kg)*m
    for (int n=0; n<3*natoms; n++){
        if (n!=n_indx){
            //printf(" %d %d %e\n", gv[n].n1, gv[n].n2, gv[n].val);
            vm[n] = sqrt(kb*temperature); // Units sqrt(kg)*m/s
            if (gv[n].val < 0.0) vm[n] = vm[n]*(-1.0);
        }
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
