#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <random>
#include "mpi.h"

#include "asr.h"
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

Asr::Asr(MC *mc) : Ptrs(mc) {
    fh_debug = fopen("DEBUG_ASR","w");
    fh_fc2 = fopen("FC2_ASR","w");
}

Asr::~Asr() 
{

    fclose(fh_debug);
    fclose(fh_fc2);

    if (order >= 3) mem->deallocate(fc3);
    if (order >= 4) mem->deallocate(fc4);

};

/*
Read FC files and store values
*/

void Asr::go(int order_in)
{
    
    order = order_in;

    natoms = lmp->atom->natoms;

    /* Read IFCs */
    readFcs();

    /* FC2 ASR */
    //fprintf(fh_fc2,"i a j b fc\n");
    double single; // this is the single term that is determined from all other terms
    int ni,nj;
    int ii,jj; // ii and jj represent an (i,a) or (j,b) pair
    printf(" natoms=%d\n", natoms);
    printf(" nfc2=%d\n", nfc2);
    fc2_struct *fc2ii;
    mem->allocate(fc2ii,9*natoms); // 8*natoms fc2ii values
    int counter=0;
    for (int i=0; i<natoms; i++){
        for (int a=0; a<3; a++){
            for (int b=0; b<3; b++){
                fc2ii[counter].i=i;
                fc2ii[counter].j=i;
                fc2ii[counter].a=a;
                fc2ii[counter].b=b;
                fc2ii[counter].val=0.0;
                counter++;
            }
        }
    }

    int ii_indx;
    int jj_indx;
    for (int i=0; i<natoms; i++){
        for (int a=0; a<3; a++){
            for (int b=0; b<3; b++){
                single = 0.0;
                for (int w=0; w<nfc2; w++){
                    if (fc2[w].i==i && fc2[w].a==a && fc2[w].b==b && (fc2[w].j!=i)){
                        //printf("%f\n", fc2[w].val);
                        ii = 3*fc2[w].i+a;
                        jj = 3*fc2[w].j+b;
                        //printf("%d %d\n", ii,jj);
                        
                        if (jj>ii){
                            // need to find fc2ii associated with fc2[w].j
                            ii_indx = 9*fc2[w].i+3*a+b; 
                            fc2ii[ii_indx].val -= fc2[w].val;
                            jj_indx = 9*fc2[w].j+3*a+b; 
                            fc2ii[jj_indx].val -= fc2[w].val;
                            //single -= fc2[w].val;
                        }
                        if (ii==jj) printf(" WTF?\n");
                    }
                }
                //fprintf(fh_fc2,"%d %d %d %d %.12e\n", i+1,a+1,i+1,b+1,single);
            }
        }
    }

    // Count number of self interactions above the tol
    int nfc2ii_tol=0;
    for (int w=0; w<9*natoms; w++){
        if (abs(fc2ii[w].val)>mc->in->tol2) nfc2ii_tol++;
    }
    printf(" nfc2ii_tol: %d\n", nfc2ii_tol);

    // Find total number of fc2 above tol
    //nfc2 += nfc2ii_tol;

    // Write self interactions if they are above the tol
    fprintf(fh_fc2,"%d\n", nfc2+nfc2ii_tol);
    for (int w=0; w<9*natoms; w++){
        if (abs(fc2ii[w].val)>mc->in->tol2) fprintf(fh_fc2, "%d %d %d %d %.12e\n", fc2ii[w].i+1,fc2ii[w].a+1,fc2ii[w].j+1,fc2ii[w].b+1,fc2ii[w].val);
    }
    // Write all other FCs
    for (int w=0; w<nfc2; w++){
        fprintf(fh_fc2, "%d %d %d %d %.12e\n", fc2[w].i+1,fc2[w].a+1,fc2[w].j+1,fc2[w].b+1,fc2[w].val);
    }

}

void Asr::readFcs()
{

    if (order >= 2){
        

        /* Read and store FC2s */
        
        ifstream fh("FC2");
        string line;

        int i,j;
        int a,b;
        double val;
        int counter = 0;

        getline(fh,line);
        stringstream ss2(line);
        ss2 >> nfc2;
        printf(" %d FC2s.\n", nfc2);
        mem->allocate(fc2,nfc2);
        while (getline(fh, line))
        {

            stringstream ss(line);
            ss >> i >> a >> j >> b >> val;
            
            fc2[counter].i=i-1;
            fc2[counter].a=a-1;
            fc2[counter].j=j-1;
            fc2[counter].b=b-1;
            fc2[counter].val=val;
            //if (counter==0) printf(" %f\n", val);
            //std::cout << line << std::endl;
            //printf("%f\n", fc3[counter].val);
            //printf(" %d\n", counter);
            counter++;
            
        }

        fh.close();
        
    }



    if (order >= 3){
        /* Read number of FC3s */
        ifstream fh("FC3_FD");
        string line;
        nfc3 = 0;
        while (getline(fh, line))
        {
            nfc3++;
        }
        nfc3 = nfc3-1; // subtract 1 for first line
        printf(" Found %d FC3s.\n", nfc3);
        fh.close();

        mem->allocate(fc3,nfc3);

        /* Read and store FC3s */
        int i,j,k;
        int a,b,c;
        double val;
        int counter = 0;
        fh.open("FC3_FD");
        getline(fh,line);
        while (getline(fh, line))
        {

            stringstream ss(line);
            ss >> i >> a >> j >> b >> k >> c >> val;
            
            fc3[counter].i=i-1;
            fc3[counter].a=a-1;
            fc3[counter].j=j-1;
            fc3[counter].b=b-1;
            fc3[counter].k=k-1;
            fc3[counter].c=c-1;
            fc3[counter].val=val;
            //if (counter==0) printf(" %f\n", val);
            
            //printf("%f\n", fc3[counter].val);
            //printf(" %d\n", counter);
            counter++;
            
        }

        fh.close();
    }

    if (order >= 4){
        /* Read number of FC3s */
        ifstream fh("FC4_FD");
        string line;
        nfc4 = 0;
        while (getline(fh, line))
        {
            nfc4++;
        }
        nfc4 = nfc4-1; // subtract 1 for first line
        printf(" Found %d FC4s.\n", nfc4);
        fh.close();

        mem->allocate(fc4,nfc4);

        /* Read and store FC4s */
        int i,j,k,l;
        int a,b,c,d;
        double val;
        int counter = 0;
        fh.open("FC4_FD");
        getline(fh,line);
        while (getline(fh, line))
        {

            stringstream ss(line);
            ss >> i >> a >> j >> b >> k >> c >> l >> d >> val;
            fc4[counter].i=i-1;
            fc4[counter].a=a-1;
            fc4[counter].j=j-1;
            fc4[counter].b=b-1;
            fc4[counter].k=k-1;
            fc4[counter].c=c-1;
            fc4[counter].l=l-1;
            fc4[counter].d=d-1;
            fc4[counter].val=val;
            //printf("%d %d %d %d %d %d %d %d\n",i,j,k,l,a,b,c,d);
            counter++;
            //printf("%f\n", val);
        }

        fh.close();

    }


}
