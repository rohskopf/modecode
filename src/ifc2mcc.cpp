#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <random>
#include "mpi.h"

#include "ifc2mcc.h"
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

using namespace LAMMPS_NS;

using namespace MC_NS;

Ifc2mcc::Ifc2mcc(MC *mc) : Ptrs(mc) {
    fh_debug = fopen("DEBUG_IFC2MCC","w");

    rank = mc->rank;

    //char debug[64];
    //sprintf (debug, "debug_ifc2mcc/DEBUG_IFC2MCC%d", rank);
    //fh_debug = fopen(debug,"w");
}

Ifc2mcc::~Ifc2mcc() 
{

    fclose(fh_debug);
    if (task==0){
        if (order >= 3) mem->deallocate(fc3);
        if (order >= 4) mem->deallocate(fc4);
        mem->deallocate(emat);
    }
    else if (task==1){
        mem->deallocate(mcc3);

    }

};

/*
Read MCC file and store values
*/

void Ifc2mcc::go(double tol)
{
    order = mc->order;

    if (rank==0) printf(" Converting IFCs to MCCs, order %d.\n", order);

    natoms = lmp->atom->natoms; // and natoms = 8
    if (rank==0) printf(" natoms: %d\n", natoms);

    /* Read FCs */
    readFcs();

    /* Read eigenvectors */
    readEmat();


    /*
    # Calculate each MCC
    fh = open('MCC3_NEW', 'w')

    for n1 in range(3,3*natoms):
        print(" n1: %d" % (n1))
        for n2 in range(n1,3*natoms):
            for n3 in range(n2,3*natoms):
                value = 0.
                for i in range(0,natoms):
                    for j in range(0,natoms):
                        for k in range(0,natoms):
                            for a in range(0,3):
                                for b in range(0,3):
                                    for c in range(0,3):
                                        value = value + (ifc3[i][j][k][a][b][c]*emat[a+i*3][n1]*emat[b+j*3][n2]*emat[c+k*3][n3])/np.sqrt(m*m*m)
                fh.write("%d %d %d %e\n" % (n1,n2,n3,value))

    fh.close()
    */

    /* Loop through all modes, atoms, directions, and calculate MCC. */

    //printf(" %e\n", lmp->atom->x[1][0]);
    int *type = lmp->atom->type;
    //printf(" %d\n", type[0]);
    double *mass = lmp->atom->mass;

    /*
    mem->allocate(mass,lmp->atom->ntypes);
    // Convert mass to SI units
    for (int t=0;t<lmp->atom->ntypes; t++){
        printf(" %d\n", t);
        mass[t] = lmp->atom->mass[t]/(6.0221409e+23*1e3);
    }
    */

    //printf(" %d\n", lmp->atom->ntypes);
    //printf(" %e\n", mass[type[0]]);

    unsigned long long int n = natoms*3 - 3; // dont count translational modes
    //int nmcc2 = nChoosek(n+2-1,n-1);
    double mass_factor;

    int ii,jj,kk,ll; // dynamical matrix indices e.g. ii = 3*i+a

    int i,j;
    int a,b;
        
    //printf("rank: %d\n", rank);
    if (order==2){
        if (rank==0){ // Do MCC2 calculations on one proc, it's easier to have in a single file.

            FILE * fh_mcc2;
            fh_mcc2 = fopen("MCC2","w");
            //double mass = 4.6637066e-26; // Assume every atom has silicon mass for now.
            printf(" natoms: %d\n", natoms);
            //int nmcc2 = nChoosek(n+2-1,n-1);
            int nmcc2 = 3*natoms-3;
            printf(" %d MCC2s.\n", nmcc2);

            fprintf(fh_mcc2, "%d\n", nmcc2);
            mass_factor = 1.0/(6.0221409e+23*1e3*6.0221409e+23*1e3); // convert 2 masses to kg in denominator

            double fc2_value;
            double mcc2_value;
            for (int n1=3; n1<3*natoms; n1++){
                //printf(" n1: %d\n", n1);
                for (int n2=n1; n2<3*natoms; n2++){
                        if (n1==n2){
                        mcc2_value = 0.0;
                        for (int f=0; f<nfc2; f++){ // Loop through nonzero FC3s
                            i = fc2[f].i;
                            j = fc2[f].j;
                            a = fc2[f].a;
                            b = fc2[f].b;
                            fc2_value = fc2[f].val*(2.179874099E-18)*1.89e+10*1.89e+10; // Convert from Ryd/Bohr^2 to J/m^2
                            //printf(" %d %d %d %d %e\n", i,a,j,b,fc2_value);
                            //printf(" %f\n", emat[a+i*3][n1]);
                            //printf(" %f\n", emat[b+j*3][n2]);
                            //printf(" %f\n", mass[type[i]]);
                            //printf(" %f\n", mass[type[j]]);
                            mcc2_value += (fc2_value*emat[a+i*3][n1]*emat[b+j*3][n2])/(sqrt(mass[type[i]]*mass[type[j]]*mass_factor));
                            //printf("ASDF\n");

                            ii=3*i+a;
                            jj=3*j+b;

                            //fprintf(fh_debug, "n1, n2: %d %d\n", n1,n2);
                            //fprintf(fh_debug, "  i j a b: %d %d %d %d\n", i,j,a,b);
                            //fprintf(fh_debug, "  fc2: %e\n", fc2_value);
                            //fprintf(fh_debug, "  term: %e\n", (fc2_value*emat[a+i*3][n1]*emat[b+j*3][n2])/(sqrt(mass*mass)));
                        }
                        fprintf(fh_mcc2, "%d %d %.20e\n", n1,n2,mcc2_value);
                    }
                }
            }

            fclose(fh_mcc2);

        } // if rank == 0
    } // if order == 2

    if (order==3){
        FILE * fh_mcc3;

        char fh_mcc3_filename[64];
        sprintf (fh_mcc3_filename, "mcc3/MCC3_%d", rank);
        fh_mcc3 = fopen(fh_mcc3_filename,"w");
        //double mass1 = 4.6637066e-26; // Assume every atom has silicon mass for now.
        //natoms = 8; // and natoms = 8

        n = natoms*3 - 3; // dont count translational modes
        //n = natoms*3; // including translational modes
        //n = 2000;
        printf(" n: %d\n", n);
        unsigned long long int nmcc3_smaller = nChoosek(n+3-1,n-1);
        //printf(" %d MCC3s with nchoosek.\n", nmcc3_smaller);
        unsigned long long int nmcc3 = n*n*n;
        if (rank==0){ 
            //printf(" %d MCC3s.\n", nmcc3);
            std::cout << "nmcc3: " << nmcc3 << std::endl;
            //std::cout << "nmcc3_smaller: " << nmcc3_smaller << std::endl;
            // Calculate number of bytes required
            std::cout << 8*nmcc3/1e9 << "GB of doubles." << std::endl;
            std::cout << nmcc3*(4+20)/1e9 << "GB of file storage." << std::endl;
            printf(" --- Using nChooseK will reduce by a factor of 6. ---\n");
        }

        // Split modes across procs

        int *nepp;
        mem->allocate(nepp, mc->nprocs); // number elements (MCC3s) per proc
        for (int p=0; p<mc->nprocs; p++){
            nepp[p] = n/mc->nprocs;
        }

        // divide up the remainder
        for (int p=0; p<(n % mc->nprocs); p++){
            nepp[p] += 1;
        }

        int start_indx = 3;
        //int start_indx=0; // including translational modes
        for (int p=0; p<rank; p++){
            start_indx += nepp[p];
        }
        int end_indx = 0; //napp[0]-1;
        for (int p=0; p<rank+1; p++){
            end_indx += nepp[p];
        }
        end_indx=end_indx+4-1; // including translational modes
        //end_indx = end_indx+1-1;

        if (rank==0){
            printf(" Splitting modes on procs like:\n");
            for (int p=0; p<mc->nprocs; p++){
                printf("  %d Modes on proc %d.\n", nepp[p],p);
            }
        }
        //printf("  Proc %d: %d-%d.\n", rank,start_indx,end_indx);


        fprintf(fh_mcc3, "%d\n", nepp[rank]*n*n); // Number of MCC3s on this proc

        mass_factor = 1.0/(6.0221409e+23*1e3*6.0221409e+23*1e3*6.0221409e+23*1e3); // convert 3 masses to kg in denominator
        //mass_factor=1.0;
        int k;
        int c;
        double fc;
        double mcc3_value;
        double sqrtmasses;
        int counter=0;
        //for (int n1=3; n1<3*natoms; n1++){
        for (int n1=start_indx; n1<end_indx; n1++){
        //for (int n1=40; n1<41; n1++){
            if (rank==0) printf(" n1: %d\n", n1);
            for (int n2=3; n2<3*natoms; n2++){

                for (int n3=3; n3<3*natoms;n3++){
                    mcc3_value = 0.0;

                    for (int f=0; f<nfc3; f++){ // Loop through nonzero FC3s

                        i = fc3[f].i;
                        j = fc3[f].j;
                        k = fc3[f].k;
                        a = fc3[f].a;
                        b = fc3[f].b;
                        c = fc3[f].c;

                        fc = fc3[f].val*(2.179874099E-18)*1.89e+10*1.89e+10*1.89e+10;
                        //printf(" %d %d %d %d %d %d %f\n", i,j,k,a,b,c,fc3_value);
                                    
                        ii=3*i+a;
                        jj=3*j+b;
                        kk=3*k+c;
                        
    
                        sqrtmasses = sqrt(mass[type[i]]*mass[type[j]]*mass[type[k]]*mass_factor);
                        //printf("%d %d %d\n",ii,jj,kk);
                        //fprintf(fh_debug, "%d %d %d %d %d %d\n",i,a,j,b,k,c);

                        mcc3_value += (fc*emat[ii][n1]*emat[jj][n2]*emat[kk][n3])/sqrtmasses;
                        //fprintf(fh_mcc3, "%d %d %d %d %d %d %e\n", i,a,j,b,k,c,fc3_value);

                     
                        

                    }
                    if (abs(mcc3_value)>tol){
                        counter++;
                        //fprintf(fh_mcc3, "%d %d %d %.20e\n", n1,n2,n3,mcc3_value);
                    }
                    if (abs(mcc3_value)>tol){ // Originally tol was 1e44.
                        fprintf(fh_mcc3, "%d %d %d %.20e\n", n1,n2,n3,mcc3_value);
                    }
                }
            }
        }

        //printf(" Found %d MCC3s > 1e48.\n", counter);


        fclose(fh_mcc3);
        mem->deallocate(nepp);
    }

    if (order ==4){

        FILE * fh_mcc4;
        char fh_mcc4_filename[64];
        sprintf (fh_mcc4_filename, "mcc4/MCC4_%d", rank);
        fh_mcc4 = fopen(fh_mcc4_filename,"w");

        n = natoms*3 - 3; // dont count translational modes
        //int nmcc3 = nChoosek(n+3-1,n-1);
        int nmcc4 = n*n*n*n;

        if (rank==0) printf(" %d MCC4s.\n", nmcc4);

        int *nepp;
        mem->allocate(nepp, mc->nprocs); // number elements (MCC3s) per proc
        for (int p=0; p<mc->nprocs; p++){
            nepp[p] = n/mc->nprocs;
        }

        // divide up the remainder
        for (int p=0; p<(n % mc->nprocs); p++){
            nepp[p] += 1;
        }

        int start_indx = 3;
        for (int p=0; p<rank; p++){
            start_indx += nepp[p];
        }
        int end_indx = 0; //napp[0]-1;
        for (int p=0; p<rank+1; p++){
            end_indx += nepp[p];
        }
        end_indx=end_indx+4-1;

        if (rank==0){
            printf(" Splitting modes on procs like:\n");
            for (int p=0; p<mc->nprocs; p++){
                printf("  %d Modes on proc %d.\n", nepp[p],p);
            }
        }
        //printf("  Proc %d: %d-%d.\n", rank,start_indx,end_indx);


        fprintf(fh_mcc4, "%d\n", nepp[rank]*n*n*n); // Number of MCC4s on this proc

        mass_factor = 1.0/(6.0221409e+23*1e3*6.0221409e+23*1e3*6.0221409e+23*1e3*6.0221409e+23*1e3); // convert 4 masses to kg in denominator
        int i,j,k,l;
        int a,b,c,d;
        double fc;
        double val;
        double sqrtmasses;
        int nmodes = 3*natoms;
        //nmodes = 4;
        //for (int n1=3; n1<nmodes; n1++){
        for (int n1=start_indx; n1<end_indx; n1++){
            if (rank==0) printf(" n1: %d\n", n1);
            for (int n2=3; n2<nmodes; n2++){
                for (int n3=3; n3<nmodes;n3++){
                    for (int n4=3; n4<nmodes; n4++){
                        val = 0.0;
                        for (int f=0; f<nfc4; f++){ // Loop through nonzero FC3s
                            ii = fc4[f].i;
                            jj = fc4[f].j;
                            kk = fc4[f].k;
                            ll = fc4[f].l;
                            a = fc4[f].a;
                            b = fc4[f].b;
                            c = fc4[f].c;
                            d = fc4[f].d;
                            fc = fc4[f].val*(2.179874099E-18)*1.89e+10*1.89e+10*1.89e+10*1.89e+10;
                            //printf(" %d %d %d %d %d %d %d %d %f\n", i,j,k,l,a,b,c,d,fc4_value);

                            i=3*ii+a;
                            j=3*jj+b;
                            k=3*kk+c;
                            l=3*ll+d;


                            sqrtmasses = (sqrt(mass[type[ii]]*mass[type[jj]]*mass[type[kk]]*mass[type[ll]]*mass_factor));
                            val += (fc*emat[i][n1]*emat[j][n2]*emat[k][n3]*emat[l][n4])/sqrtmasses;
                            

                        }
                        fprintf(fh_mcc4, "%d %d %d %d %.20e\n", n1,n2,n3,n4,val);
                    }
                }
            }
        }

        fclose(fh_mcc4);
    }
    

}

/*
Extract MCC3s for some desired modes, and write them in a MCC3 file. 
*/
void Ifc2mcc::extract()
{

    FILE * fh_mcc3;
    fh_mcc3 = fopen("MCC3","w");
    int nfiles = 1000; // Number of files to scan through, ID starting at zero.
    int nmodes = 2400; //-3; // Number of modes.
    int nmodes_i = highmode-lowmode+1; // Number of modes "i" we are concerned with in extract MCC3s for. 
    //fprintf(fh_mcc3, "%d\n", nmodes_i*nmodes*nmodes);
    for (int f=0; f<nfiles; f++){

        printf(" File: %d\n", f);

        /* Read file and look for relevant mode. */

        char filename[64];
        sprintf (filename, "mcc3/MCC3_%d", f);
        //fh_debug = fopen(debug,"w");

        ifstream fh(filename);
        string line;

        int i,j,k;
        double mcc;
        int counter = 0;

        getline(fh,line); // Ignore first line.
        //stringstream ss2(line);
        //ss2 >> nfc3;
        //if (rank==0) printf(" %d FC3s.\n", nfc3);
        //mem->allocate(fc3,nfc3);
        while (getline(fh, line))
        {

            stringstream ss(line);
            ss >> i >> j >> k >> mcc;
            if (counter % 10000 == 0){
                printf("  i: %d\n", i);
            }

            if (lowmode <= i && i <= highmode){
                if (abs(mcc)>1e47){ // originally tol was 1e47
                    fprintf(fh_mcc3, "%d %d %d %.15e\n", i,j,k,mcc);
                }
            }
            
            counter++;
            
        }

        fh.close();

        

    }
    fclose(fh_mcc3);

}

/*
Extract MCC3s for some desired modes, and write them in a MCC3 file. 
*/
void Ifc2mcc::average(int startfile,int endfile)
{
    int nmodes = 2400;

    double **mcc3_avg;
    mem->allocate(mcc3_avg,nmodes,nmodes);
    
    for (int i=0; i<nmodes; i++){
        for (int j=0; j<nmodes; j++){
            mcc3_avg[i][j]=0.0;
        }
    }

    /* Read ILIST. First number is the number of modes i. */

    int ni;
    
    ifstream fh("ILIST");
    string line;

    int counter = 0;

    getline(fh,line);
    stringstream ss(line);
    ss >> ni;
    if (rank==0) printf(" Taking average of %d MCC3 tensors.\n", ni);
    int *ilist;
    mem->allocate(ilist,ni);
    while (getline(fh, line))
    {

        stringstream ss(line);
        ss >> ilist[counter];
        
        counter++;
        
    }

    fh.close();

    /* Read through MCC3 files */
    bool insideList;
    //int startfile = 1;
    //int endfile = 3;
    for (int f=startfile; f<endfile+1; f++){

        printf(" File: %d\n", f);

        /* Read file and look for relevant mode. */

        char filename[64];
        sprintf (filename, "mcc3/MCC3_%d", f);
        //fh_debug = fopen(debug,"w");

        ifstream fh(filename);
        string line;

        int i,j,k;
        double mcc;
        int counter = 0;

        getline(fh,line); // Ignore first line.
        //stringstream ss2(line);
        //ss2 >> nfc3;
        //if (rank==0) printf(" %d FC3s.\n", nfc3);
        //mem->allocate(fc3,nfc3);
        while (getline(fh, line))
        {

            stringstream ss(line);
            ss >> i >> j >> k >> mcc;
            if (counter % 1000000 == 0){
                printf("  i in file: %d\n", i);
            }

            insideList = false;
            for (int ii=0;ii<ni;ii++){
                if (i==ilist[ii]) insideList=true;
            }
            if (insideList){
                stringstream ss(line);
                ss >> i >> j >> k >> mcc;
                mcc3_avg[j][k]+=mcc/ni;
            }
            
            counter++;
            
        }

        fh.close();

        

    }
    
    //printf(" mcc3_avg[0][0]: %f\n", mcc3_avg[0][0]);
    //printf(" mcc3_avg[0][1]: %f\n", mcc3_avg[0][1]);
    FILE * fh_avg;
    fh_avg=fopen("MCC3_AVG","w");
    for (int i=0; i<nmodes; i++){
        for (int j=0; j<nmodes; j++){
            fprintf(fh_avg,"%d %d %e\n", i,j,mcc3_avg[i][j]);
        }
    }
    fclose(fh_avg);
    
    mem->deallocate(mcc3_avg);
    mem->deallocate(ilist);
}

/*
Extract a few MCC3s associated with a group of modes in a list, ILIST.
*/
void Ifc2mcc::extractFew(int startfile, int endfile)
{

    // Extract list of modes i.
    int ni;
    
    ifstream fh("ILIST");
    string line;

    int counter = 0;

    getline(fh,line);
    stringstream ss(line);
    ss >> ni;
    if (rank==0) printf(" Extracting %d MCC3s.\n", ni*ni*ni);
    int *ilist;
    mem->allocate(ilist,ni);
    while (getline(fh, line))
    {

        stringstream ss(line);
        ss >> ilist[counter];
        
        counter++;
        
    }

    fh.close();


    FILE * fh_mcc3;
    fh_mcc3 = fopen("MCC3","w");
    int nfiles = 4; // Number of files to scan through, ID starting at zero.
    int nmodes = 2400-3; // Number of modes.
    //int nmodes_i = highmode-lowmode+1; // Number of modes "i" we are concerned with in extract MCC3s for. 
    fprintf(fh_mcc3, "%d\n", ni*ni*ni);
    for (int f=startfile; f<endfile; f++){

        printf(" File: %d\n", f);

        /* Read file and look for relevant mode. */

        char filename[64];
        sprintf (filename, "mcc3/MCC3_%d", f);
        //fh_debug = fopen(debug,"w");

        ifstream fh(filename);
        bool isopen = fh.is_open();
        if (!isopen){ 
            printf(" Breaking becaues reached file limit!\n");
            break;
        }
        string line;

        int i,j,k;
        double mcc;
        int counter = 0;

        getline(fh,line); // Ignore first line.
        //stringstream ss2(line);
        //ss2 >> nfc3;
        //if (rank==0) printf(" %d FC3s.\n", nfc3);
        //mem->allocate(fc3,nfc3);
        while (getline(fh, line))
        {

            stringstream ss(line);
            ss >> i >> j >> k >> mcc;
            if (counter % 10000 == 0){
                printf("  i: %d\n", i);
            }

            // Check if i,j, AND k are in ilist.
            bool i_in = false;
            bool j_in = false;
            bool k_in = false;
            for (int ii=0; ii<ni; ii++){
                if (i==ilist[ii]) i_in=true;
                if (j==ilist[ii]) j_in=true;
                if (k==ilist[ii]) k_in=true;
            }
            if (i_in && j_in && k_in){
                if (abs(mcc)>1e45){
                    fprintf(fh_mcc3, "%d %d %d %.15e\n", i,j,k,mcc);
                }
            }
            
            counter++;
            
        }

        fh.close();

        

    }
    fclose(fh_mcc3);

}

/*
Calculate MCC3s associated with a few modes in ILIST.
*/
void Ifc2mcc::calcFew()
{
    printf(" Calculating a few MCCs based on ILIST\n");
    // Extract list of modes i.
    int ni;
    
    ifstream fh("ILIST");
    string line;

    int counter = 0;

    getline(fh,line);
    stringstream ss(line);
    ss >> ni;
    if (rank==0) printf(" Extracting %d MCC3s.\n", ni*ni*ni);
    int *ilist;
    mem->allocate(ilist,ni);
    while (getline(fh, line))
    {

        stringstream ss(line);
        ss >> ilist[counter];
        
        counter++;
        
    }

    fh.close();

    // Proceed with calculating MCC3

    order = mc->order;

    if (rank==0) printf(" Converting IFCs to MCCs, order %d.\n", order);

    natoms = lmp->atom->natoms; // and natoms = 8
    if (rank==0) printf(" natoms: %d\n", natoms);

    /* Read FCs */
    readFcs();

    /* Read eigenvectors */
    readEmat();

    /* Loop through all modes, atoms, directions, and calculate MCC. */

    //printf(" %e\n", lmp->atom->x[1][0]);
    int *type = lmp->atom->type;
    //printf(" %d\n", type[0]);
    double *mass = lmp->atom->mass;

    unsigned long long int n = natoms*3 - 3; // dont count translational modes
    //int nmcc2 = nChoosek(n+2-1,n-1);
    double mass_factor;

    int ii,jj,kk,ll; // dynamical matrix indices e.g. ii = 3*i+a

    int i,j;
    int a,b;

    if (order==3){
        FILE * fh_mcc3;

        char fh_mcc3_filename[64];
        //sprintf (fh_mcc3_filename, "mcc3/MCC3_%d", rank);
        fh_mcc3 = fopen("MCC3_calcFew","w");
        //double mass1 = 4.6637066e-26; // Assume every atom has silicon mass for now.
        //natoms = 8; // and natoms = 8

        n = natoms*3 - 3; // dont count translational modes
        //n = 2000;
        printf(" n: %d\n", n);
        unsigned long long int nmcc3_smaller = nChoosek(n+3-1,n-1);
        //printf(" %d MCC3s with nchoosek.\n", nmcc3_smaller);
        unsigned long long int nmcc3 = ni*ni*ni;
        if (rank==0){ 
            //printf(" %d MCC3s.\n", nmcc3);
            std::cout << "nmcc3: " << nmcc3 << std::endl;
            //std::cout << "nmcc3_smaller: " << nmcc3_smaller << std::endl;
            // Calculate number of bytes required
            std::cout << 8*nmcc3/1e9 << "GB of doubles." << std::endl;
            std::cout << nmcc3*(4+20)/1e9 << "GB of file storage." << std::endl;
            printf(" --- Using nChooseK will reduce by a factor of 6. ---\n");
        }

        fprintf(fh_mcc3, "%d\n", nmcc3); // Number of MCC3s on this proc

        mass_factor = 1.0/(6.0221409e+23*1e3*6.0221409e+23*1e3*6.0221409e+23*1e3); // convert 3 masses to kg in denominator
        //mass_factor=1.0;
        int k;
        int c;
        double fc;
        double mcc3_value;
        double sqrtmasses;
        int counter=0;
        //for (int n1=3; n1<3*natoms; n1++){
        for (int n1=3; n1<3*natoms; n1++){
        //for (int n1=40; n1<41; n1++){
            if (rank==0) printf(" n1: %d\n", n1);
            for (int n2=3; n2<3*natoms; n2++){

                for (int n3=3; n3<3*natoms;n3++){
                    
                    // Check if n1,n2,AND n3 are in ilist.
                    bool n1_in = false;
                    bool n2_in = false;
                    bool n3_in = false;
                    for (int ii=0; ii<ni; ii++){
                        if (n1==ilist[ii]) n1_in=true;
                        if (n2==ilist[ii]) n2_in=true;
                        if (n3==ilist[ii]) n3_in=true;
                    }
                    if (n1_in && n2_in && n3_in){
                    
                        mcc3_value = 0.0;

                        for (int f=0; f<nfc3; f++){ // Loop through nonzero FC3s

                            i = fc3[f].i;
                            j = fc3[f].j;
                            k = fc3[f].k;
                            a = fc3[f].a;
                            b = fc3[f].b;
                            c = fc3[f].c;

                            fc = fc3[f].val*(2.179874099E-18)*1.89e+10*1.89e+10*1.89e+10;
                            //printf(" %d %d %d %d %d %d %f\n", i,j,k,a,b,c,fc3_value);
                                        
                            ii=3*i+a;
                            jj=3*j+b;
                            kk=3*k+c;
                            
        
                            sqrtmasses = sqrt(mass[type[i]]*mass[type[j]]*mass[type[k]]*mass_factor);
                            //printf("%d %d %d\n",ii,jj,kk);
                            //fprintf(fh_debug, "%d %d %d %d %d %d\n",i,a,j,b,k,c);

                            mcc3_value += (fc*emat[ii][n1]*emat[jj][n2]*emat[kk][n3])/sqrtmasses;
                            //fprintf(fh_mcc3, "%d %d %d %d %d %d %e\n", i,a,j,b,k,c,fc3_value);

                         
                            

                        }
                        if (abs(mcc3_value)>1e48){
                            counter++;
                            //fprintf(fh_mcc3, "%d %d %d %.20e\n", n1,n2,n3,mcc3_value);
                        }
                        fprintf(fh_mcc3, "%d %d %d %.20e\n", n1,n2,n3,mcc3_value);

                    }

                }
            }
        }

        //printf(" Found %d MCC3s > 1e48.\n", counter);


        fclose(fh_mcc3);
    }


}

/*
Read MCC file and store values
*/

void Ifc2mcc::go_spatial()
{
    order = mc->order;

    if (rank==0) printf(" Converting IFCs to spatial MCCs, order %d.\n", order);

    natoms = lmp->atom->natoms; // and natoms = 8
    if (rank==0) printf(" natoms: %d\n", natoms);

    /* Read FCs */
    readFcs();

    /* Read eigenvectors */
    readEmat();


    /*
    # Calculate each MCC
    fh = open('MCC3_NEW', 'w')

    for n1 in range(3,3*natoms):
        print(" n1: %d" % (n1))
        for n2 in range(n1,3*natoms):
            for n3 in range(n2,3*natoms):
                value = 0.
                for i in range(0,natoms):
                    for j in range(0,natoms):
                        for k in range(0,natoms):
                            for a in range(0,3):
                                for b in range(0,3):
                                    for c in range(0,3):
                                        value = value + (ifc3[i][j][k][a][b][c]*emat[a+i*3][n1]*emat[b+j*3][n2]*emat[c+k*3][n3])/np.sqrt(m*m*m)
                fh.write("%d %d %d %e\n" % (n1,n2,n3,value))

    fh.close()
    */

    /* Loop through all modes, atoms, directions, and calculate MCC. */

    //printf(" %e\n", lmp->atom->x[1][0]);
    int *type = lmp->atom->type;
    int *tag = lmp->atom->tag;
    //printf(" %d\n", type[0]);
    double *mass = lmp->atom->mass;

    /*
    mem->allocate(mass,lmp->atom->ntypes);
    // Convert mass to SI units
    for (int t=0;t<lmp->atom->ntypes; t++){
        printf(" %d\n", t);
        mass[t] = lmp->atom->mass[t]/(6.0221409e+23*1e3);
    }
    */

    //printf(" %d\n", lmp->atom->ntypes);
    //printf(" %e\n", mass[type[0]]);

    unsigned long long int n = natoms*3 - 3; // dont count translational modes
    //int nmcc2 = nChoosek(n+2-1,n-1);
    double mass_factor;

    int ii,jj,kk,ll; // dynamical matrix indices e.g. ii = 3*i+a

    int i,j;
    int a,b;
        
    //printf("rank: %d\n", rank);
    if (order==2){
        //if (rank==0){ // Do MCC2 calculations on one proc, it's easier to have in a single file.

            FILE * fh_mcc2;
            //fh_mcc2 = fopen("MCC2","w");
            //double mass = 4.6637066e-26; // Assume every atom has silicon mass for now.
            if (rank==0) printf(" natoms: %d\n", natoms);
            //int nmcc2 = nChoosek(n+2-1,n-1);
            int nmcc2 = (3*natoms)*(3*natoms);
            //printf(" %d MCC2s.\n", nmcc2);

            char fh_mcc2_filename[64];
            sprintf (fh_mcc2_filename, "smcc2/SMCC2_%d", rank);
            fh_mcc2 = fopen(fh_mcc2_filename,"w");
            //double mass1 = 4.6637066e-26; // Assume every atom has silicon mass for now.
            //natoms = 8; // and natoms = 8

            int n = natoms*3; // dont count translational modes

            // Split modes across procs

            int *nepp;
            mem->allocate(nepp, mc->nprocs); // number elements (MCC3s) per proc
            for (int p=0; p<mc->nprocs; p++){
                nepp[p] = n/mc->nprocs;
            }

            // divide up the remainder
            for (int p=0; p<(n % mc->nprocs); p++){
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

            if (rank==0){
                printf(" Splitting modes on procs like:\n");
                for (int p=0; p<mc->nprocs; p++){
                    printf("  %d Modes on proc %d.\n", nepp[p],p);
                }
            }
            //printf("  Proc %d: %d-%d.\n", rank,start_indx,end_indx);


            fprintf(fh_mcc2, "%d\n", nepp[rank]*n); // Number of MCC2s on this proc

            //fprintf(fh_mcc2, "%d\n", nmcc2);
            mass_factor = 1.0/(6.0221409e+23*1e3*6.0221409e+23*1e3); // convert 2 masses to kg in denominator

            double fc2_value;
            double mcc2_value;
            double K1, K2; // K_n1n2_AB and K_n2n1_BA
            double mi,mj;
            double phi;
            double term1,term2;
            int itype,jtype;
            for (int n1=start_indx; n1<end_indx; n1++){
            //for (int n1=10360; n1<10361; n1++){
                if (rank==0){
                    if ( (n1 % 1)== 0) printf(" n1: %d\n", n1);
                }
                for (int n2=0; n2<3*natoms; n2++){
                        //if (n1==n2){
                        mcc2_value = 0.0;
                        K1 = 0.0;
                        K2 = 0.0;
                        for (int f=0; f<nfc2; f++){ // Loop through nonzero FC3s
                            i = fc2[f].i;
                            j = fc2[f].j;
                            a = fc2[f].a;
                            b = fc2[f].b;
                            phi = fc2[f].val*(2.179874099E-18)*1.89e+10*1.89e+10; // Convert from Ryd/Bohr^2 to J/m^2
                            //phi = fc2[f].val*13.605698066*(1./0.529177249)*(1./0.529177249); // Convert from Ryd/Bohr^2 to eV/A^2, THIS IS FOR COMPARING TO PAIR_TEP.CPP
                            //printf(" %d %d %d %d %e\n", i,a,j,b,fc2_value);
                            //printf(" %f\n", emat[a+i*3][n1]);
                            //printf(" %f\n", emat[b+j*3][n2]);
                            //printf(" %f\n", mass[type[i]]);
                            //printf(" %f\n", mass[type[j]]);

                            itype = type[i];
                            jtype = type[j];
                            mi = mass[type[i]]*(1.0/(6.02214076e23*1e3)); // atom mass in kg
                            mj = mass[type[j]]*(1.0/(6.02214076e23*1e3)); // atom mass in kg

                            term1 = (phi*emat[3*j+b][n2]*emat[3*i+a][n1])/sqrt(mi*mj); // (J/(m^2 * kg)
                            term2 = (phi*emat[3*i+b][n1]*emat[3*j+a][n2])/sqrt(mi*mj);

                            // Calculate coupling constants.
                            //if (itype==1 && jtype==2){
                            if (itype==1 && jtype==2){
                              mcc2_value += 1.0*term1;
                              /*
                              if (n1==0 && n2==0){
                                //if (i != j){
                                  if (abs(phi) > 0.0){
                                    fprintf(fh_debug, "%d %d %d %d %e\n", i,a,j,b,term1) ;
                                    fprintf(fh_debug, "  %e %e %e %e\n", phi,emat[3*j+b][n2],emat[3*i+a][n1],sqrt(mi*mj) );
                                    fprintf(fh_debug, "  %d %d\n", type[i], type[j]);
                                  }
                                //}
                              }
                              */
                            }
                            //else if (itype==2 && jtype==1){
                            else if (itype==2 && jtype==1){
                              mcc2_value -= 1.0*term2; // THIS WORKS FOR SOME REASON -> Don't need K_n2n1????
                              //k_n2n1 += 0.5*term3; // THIS ALSO WORKS.
                              /*
                              if (n1==0 && n2==0){
                                //if (i != j){
                                  if (abs(phi) > 0.0){
                                    fprintf(fh_debug, "%d %d %d %d %e\n", i,a,j,b,term2) ;
                                    fprintf(fh_debug, "  %e %e %e %e\n", phi,emat[3*j+b][n2],emat[3*i+a][n1],sqrt(mi*mj) );
                                    fprintf(fh_debug, "  %d %d\n", type[i], type[j]);
                                  }
                                //}
                              }
                              */

                            }
                            
                            //else if (type[i]==2 && type[j]
                            //printf("ASDF\n");

                            ii=3*i+a;
                            jj=3*j+b;

                            //fprintf(fh_debug, "n1, n2: %d %d\n", n1,n2);
                            //fprintf(fh_debug, "  i j a b: %d %d %d %d\n", i,j,a,b);
                            //fprintf(fh_debug, "  fc2: %e\n", fc2_value);
                            //fprintf(fh_debug, "  term: %e\n", (fc2_value*emat[a+i*3][n1]*emat[b+j*3][n2])/(sqrt(mass*mass)));
                        }

                        //if (n1==0 && n2==0){
                          //if (i != j){
                            //fprintf(fh_debug, "%d %d %e\n", n1,n2,mcc2_value);
                          //}
                        //}

                        //I used to divide by 2, but we can do this in the heat transfer expression instead.
                        //if (n1 != n2) mcc2_value = 0.5*mcc2_value;

                        //if (abs(mcc2_value) > 1e23){
                        if (abs(mcc2_value)>0){
                            fprintf(fh_mcc2, "%d %d %.20e\n", n1,n2,mcc2_value);
                        }
                    //}
                }
            }

            fclose(fh_mcc2);

        //} // if rank == 0
    } // if order == 2

    if (order==3){
        FILE * fh_mcc3;

        char fh_mcc3_filename[64];
        sprintf (fh_mcc3_filename, "mcc3/MCC3_%d", rank);
        fh_mcc3 = fopen(fh_mcc3_filename,"w");
        //double mass1 = 4.6637066e-26; // Assume every atom has silicon mass for now.
        //natoms = 8; // and natoms = 8

        n = natoms*3; // - 3; // dont count translational modes
        //n = 2000;
        //printf(" n: %d\n", n);
        unsigned long long int nmcc3_smaller = nChoosek(n+3-1,n-1);
        //printf(" %d MCC3s with nchoosek.\n", nmcc3_smaller);
        unsigned long long int nmcc3 = n*n*n;
        if (rank==0){ 
            //printf(" %d MCC3s.\n", nmcc3);
            std::cout << "nmcc3: " << nmcc3 << std::endl;
            //std::cout << "nmcc3_smaller: " << nmcc3_smaller << std::endl;
            // Calculate number of bytes required
            std::cout << 8*nmcc3/1e9 << "GB of doubles." << std::endl;
            std::cout << nmcc3*(4+20)/1e9 << "GB of file storage." << std::endl;
            printf(" --- Using nChooseK will reduce by a factor of 6. ---\n");
        }

        // Split modes across procs

        int *nepp;
        mem->allocate(nepp, mc->nprocs); // number elements (MCC3s) per proc
        for (int p=0; p<mc->nprocs; p++){
            nepp[p] = n/mc->nprocs;
        }

        // divide up the remainder
        for (int p=0; p<(n % mc->nprocs); p++){
            nepp[p] += 1;
        }

        int start_indx = 0; //3;
        for (int p=0; p<rank; p++){
            start_indx += nepp[p];
        }
        int end_indx = 0; //napp[0]-1;
        for (int p=0; p<rank+1; p++){
            end_indx += nepp[p];
        }
        end_indx=end_indx; //+4-1;

        if (rank==0){
            printf(" Splitting modes on procs like:\n");
            for (int p=0; p<mc->nprocs; p++){
                printf("  %d Modes on proc %d.\n", nepp[p],p);
            }
        }
        //printf("  Proc %d: %d-%d.\n", rank,start_indx,end_indx);


        //fprintf(fh_mcc3, "%d\n", nepp[rank]*n*n); // Number of MCC3s on this proc

        mass_factor = 1.0/(6.0221409e+23*1e3*6.0221409e+23*1e3*6.0221409e+23*1e3); // convert 3 masses to kg in denominator
        //mass_factor=1.0;
        int k;
        int c;
        double fc;
        double mi,mj,mk;
        double mcc3_value;
        double sqrtmasses;
        int counter=0;
        double term1;
        //for (int n1=3; n1<3*natoms; n1++){
        for (int n1=start_indx; n1<end_indx; n1++){
        //for (int n1=40; n1<41; n1++){
            if (rank==0) printf(" n1: %d\n", n1);
            for (int n2=0; n2<3*natoms; n2++){

                for (int n3=0; n3<3*natoms;n3++){
                    mcc3_value = 0.0;

                    for (int f=0; f<nfc3; f++){ // Loop through nonzero FC3s
                        //printf("%d %d %d %d\n", n1,n2,n3,f);
                        i = fc3[f].i;
                        j = fc3[f].j;
                        k = fc3[f].k;
                        a = fc3[f].a;
                        b = fc3[f].b;
                        c = fc3[f].c;

                        fc = fc3[f].val*(2.179874099E-18)*1.89e+10*1.89e+10*1.89e+10;
                        //printf(" FC: %d %d %d %d %d %d %f\n", i,a,j,b,k,c,fc);
                                    
                        ii=3*i+a;
                        jj=3*j+b;
                        kk=3*k+c;
                        

                        //printf("%d %d %d %e %e %e\n", type[i],type[j],type[k],mass[type[i]],mass[type[j]],mass[type[k]]);
                        mi = mass[type[i]]*(1.0/(6.02214076e23*1e3)); // atom mass in kg
                        mj = mass[type[j]]*(1.0/(6.02214076e23*1e3)); // atom mass in kg
                        mk = mass[type[k]]*(1.0/(6.02214076e23*1e3)); // atom mass in kg
    
                        sqrtmasses = sqrt(mass[type[i]]*mass[type[j]]*mass[type[k]]*mass_factor);
                        //printf("%d %d %d\n",ii,jj,kk);
                        //fprintf(fh_debug, "%d %d %d %d %d %d\n",i,a,j,b,k,c);
                        //printf("About to calculate termu %d %d %d.\n",ii,jj,kk);

                        term1 = (fc*emat[3*i+a][n1]*emat[3*j+b][n2]*emat[3*k+c][n3])/sqrt(mi*mj*mk); // (J/(m^3 * kg^(3/2))

                        if (type[i]==1 && type[j]==2){
                          mcc3_value += 1.0*term1;
                        }
                        else if (type[i]==2 && type[j]==1){
                          mcc3_value -= 1.0*term1;
                        }              

                        //mcc3_value += (fc*emat[ii][n1]*emat[jj][n2]*emat[kk][n3])/sqrtmasses;
                        //fprintf(fh_mcc3, "%d %d %d %d %d %d %e\n", i,a,j,b,k,c,fc3_value);

                     
                        

                    }
                    if (abs(mcc3_value)>1e48){
                        counter++;
                        //fprintf(fh_mcc3, "%d %d %d %.20e\n", n1,n2,n3,mcc3_value);
                    }
                    fprintf(fh_mcc3, "%d %d %d %.20e\n", n1,n2,n3,mcc3_value);

                }
            }
        }

        //printf(" Found %d MCC3s > 1e48.\n", counter);


        fclose(fh_mcc3);
        mem->deallocate(nepp);
    }

    if (order ==4){

        FILE * fh_mcc4;
        char fh_mcc4_filename[64];
        sprintf (fh_mcc4_filename, "mcc4/MCC4_%d", rank);
        fh_mcc4 = fopen(fh_mcc4_filename,"w");

        n = natoms*3 - 3; // dont count translational modes
        //int nmcc3 = nChoosek(n+3-1,n-1);
        int nmcc4 = n*n*n*n;

        if (rank==0) printf(" %d MCC4s.\n", nmcc4);

        int *nepp;
        mem->allocate(nepp, mc->nprocs); // number elements (MCC3s) per proc
        for (int p=0; p<mc->nprocs; p++){
            nepp[p] = n/mc->nprocs;
        }

        // divide up the remainder
        for (int p=0; p<(n % mc->nprocs); p++){
            nepp[p] += 1;
        }

        int start_indx = 3;
        for (int p=0; p<rank; p++){
            start_indx += nepp[p];
        }
        int end_indx = 0; //napp[0]-1;
        for (int p=0; p<rank+1; p++){
            end_indx += nepp[p];
        }
        end_indx=end_indx+4-1;

        if (rank==0){
            printf(" Splitting modes on procs like:\n");
            for (int p=0; p<mc->nprocs; p++){
                printf("  %d Modes on proc %d.\n", nepp[p],p);
            }
        }
        //printf("  Proc %d: %d-%d.\n", rank,start_indx,end_indx);


        fprintf(fh_mcc4, "%d\n", nepp[rank]*n*n*n); // Number of MCC4s on this proc

        mass_factor = 1.0/(6.0221409e+23*1e3*6.0221409e+23*1e3*6.0221409e+23*1e3*6.0221409e+23*1e3); // convert 4 masses to kg in denominator
        int i,j,k,l;
        int a,b,c,d;
        double fc;
        double val;
        double sqrtmasses;
        int nmodes = 3*natoms;
        //nmodes = 4;
        //for (int n1=3; n1<nmodes; n1++){
        for (int n1=start_indx; n1<end_indx; n1++){
            if (rank==0) printf(" n1: %d\n", n1);
            for (int n2=3; n2<nmodes; n2++){
                for (int n3=3; n3<nmodes;n3++){
                    for (int n4=3; n4<nmodes; n4++){
                        val = 0.0;
                        for (int f=0; f<nfc4; f++){ // Loop through nonzero FC3s
                            ii = fc4[f].i;
                            jj = fc4[f].j;
                            kk = fc4[f].k;
                            ll = fc4[f].l;
                            a = fc4[f].a;
                            b = fc4[f].b;
                            c = fc4[f].c;
                            d = fc4[f].d;
                            fc = fc4[f].val*(2.179874099E-18)*1.89e+10*1.89e+10*1.89e+10*1.89e+10;
                            //printf(" %d %d %d %d %d %d %d %d %f\n", i,j,k,l,a,b,c,d,fc4_value);

                            i=3*ii+a;
                            j=3*jj+b;
                            k=3*kk+c;
                            l=3*ll+d;


                            sqrtmasses = (sqrt(mass[type[ii]]*mass[type[jj]]*mass[type[kk]]*mass[type[ll]]*mass_factor));
                            val += (fc*emat[i][n1]*emat[j][n2]*emat[k][n3]*emat[l][n4])/sqrtmasses;
                            

                        }
                        fprintf(fh_mcc4, "%d %d %d %d %.20e\n", n1,n2,n3,n4,val);
                    }
                }
            }
        }

        fclose(fh_mcc4);
    }
    

}

/*
Read FCs. 
*/

void Ifc2mcc::readFcs()
{

    if (order >= 2){

        /* Read and store FC2s */
        
        ifstream fh("FC2_ASR");
        string line;

        int i,j;
        int a,b;
        double val;
        int counter = 0;

        getline(fh,line);
        stringstream ss2(line);
        ss2 >> nfc2;
        if (rank==0) printf(" %d FC2s.\n", nfc2);
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

        /* Read and store FC3s */

        ifstream fh("FC3");
        string line;

        int i,j,k;
        int a,b,c;
        double val;
        int counter = 0;

        getline(fh,line);
        stringstream ss2(line);
        ss2 >> nfc3;
        if (rank==0) printf(" %d FC3s.\n", nfc3);
        mem->allocate(fc3,nfc3);
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

        /* Read and store FC4s */

        ifstream fh("FC4");
        string line;

        int i,j,k,l;
        int a,b,c,d;
        double val;
        int counter = 0;

        getline(fh,line);
        stringstream ss2(line);
        ss2 >> nfc4;
        if (rank==0) printf(" %d FC4s.\n", nfc4);
        mem->allocate(fc4,nfc4);
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

void Ifc2mcc::go_n1(int n1)
{

  order = 3;

  natoms = lmp->atom->natoms;

  /* Read FCs */
  readFcs();

  /* Read eigenvectors */
  readEmat();

  printf("Done reading EMAT\n");

  //printf(" %e\n", lmp->atom->x[1][0]);
  int *type = lmp->atom->type;
  //printf(" %d\n", type[0]);
  double *mass = lmp->atom->mass;

  
  FILE * fh_mcc3;

  char fh_mcc3_filename[64];
  sprintf (fh_mcc3_filename, "mcc3/MCC3_%d", rank);
  fh_mcc3 = fopen(fh_mcc3_filename,"w");
  //double mass1 = 4.6637066e-26; // Assume every atom has silicon mass for now.
  //natoms = 8; // and natoms = 8

  int n = natoms*3 - 3; // dont count translational modes
  //n = 2000;
  //printf(" n: %d\n", n);
  unsigned long long int nmcc3_smaller = nChoosek(n+3-1,n-1);
  //printf(" %d MCC3s with nchoosek.\n", nmcc3_smaller);
  unsigned long long int nmcc3 = n*n*n;
  if (rank==0){ 
      printf(" --- Calculating MCC3s for n1=%d. ---\n", n1);
  }

  // Split modes across procs

  int *nepp;
  mem->allocate(nepp, mc->nprocs); // number elements (MCC3s) per proc
  for (int p=0; p<mc->nprocs; p++){
      nepp[p] = n/mc->nprocs;
  }

  // divide up the remainder
  for (int p=0; p<(n % mc->nprocs); p++){
      nepp[p] += 1;
  }

  int start_indx = 3;
  for (int p=0; p<rank; p++){
      start_indx += nepp[p];
  }
  int end_indx = 0; //napp[0]-1;
  for (int p=0; p<rank+1; p++){
      end_indx += nepp[p];
  }
  end_indx=end_indx+4-1;

  if (rank==0){
      printf(" Splitting modes on procs like:\n");
      for (int p=0; p<mc->nprocs; p++){
          printf("  %d Modes on proc %d.\n", nepp[p],p);
      }
  }
  //printf("  Proc %d: %d-%d.\n", rank,start_indx,end_indx);


  fprintf(fh_mcc3, "%d\n", nepp[rank]*n*n); // Number of MCC3s on this proc

  double mass_factor = 1.0/(6.0221409e+23*1e3*6.0221409e+23*1e3*6.0221409e+23*1e3); // convert 3 masses to kg in denominator
  //mass_factor=1.0;

  double fc;
  double mcc3_value;
  double sqrtmasses;
  int counter=0;
  int i,j,k,a,b,c;
  int ii,jj,kk;

  for (int n2=start_indx; n2<end_indx; n2++){
      if (rank==0) printf(" n2: %d\n", n2);
      for (int n3=3; n3<3*natoms;n3++){
          mcc3_value = 0.0;

          for (int f=0; f<nfc3; f++){ // Loop through nonzero FC3s

              i = fc3[f].i;
              j = fc3[f].j;
              k = fc3[f].k;
              a = fc3[f].a;
              b = fc3[f].b;
              c = fc3[f].c;

              fc = fc3[f].val*(2.179874099E-18)*1.89e+10*1.89e+10*1.89e+10;
              //printf(" %d %d %d %d %d %d %f\n", i,j,k,a,b,c,fc3_value);
                          
              ii=3*i+a;
              jj=3*j+b;
              kk=3*k+c;
              

              sqrtmasses = sqrt(mass[type[i]]*mass[type[j]]*mass[type[k]]*mass_factor);
              //printf("%d %d %d\n",ii,jj,kk);
              //fprintf(fh_debug, "%d %d %d %d %d %d\n",i,a,j,b,k,c);

              mcc3_value += (fc*emat[ii][n1]*emat[jj][n2]*emat[kk][n3])/sqrtmasses;
              //fprintf(fh_mcc3, "%d %d %d %d %d %d %e\n", i,a,j,b,k,c,fc3_value);

           
              

          }
          if (abs(mcc3_value)>1e48){
              counter++;
              //fprintf(fh_mcc3, "%d %d %d %.20e\n", n1,n2,n3,mcc3_value);
          }
          if (abs(mcc3_value)>1e44){
              fprintf(fh_mcc3, "%d %d %d %.20e\n", n1,n2,n3,mcc3_value);
          }
      }
  }

  //printf(" Found %d MCC3s > 1e48.\n", counter);


  fclose(fh_mcc3);
  mem->deallocate(nepp);

}

/*
Extract SMCC2s and write them in a SMCC2 file. 
*/
void Ifc2mcc::extract_smcc2(int nfiles)
{

    FILE * fh_mcc2;
    fh_mcc2 = fopen("SMCC2","w");
    natoms = lmp->atom->natoms; 
    int nmodes = 3*natoms; //-3; // Number of modes.
    //int nmodes_i = highmode-lowmode+1; // Number of modes "i" we are concerned with in extract MCC3s for. 
    //fprintf(fh_mcc3, "%d\n", nmodes_i*nmodes*nmodes);
    
    for (int f=0; f<nfiles; f++){

        printf(" File: %d\n", f);

        /* Read file and look for relevant mode. */

        char filename[64];
        sprintf (filename, "smcc2/SMCC2_%d", f);
        //fh_debug = fopen(debug,"w");

        ifstream fh(filename);
        string line;

        int i,j;
        double mcc;
        int counter = 0;

        getline(fh,line); // Ignore first line.
        //stringstream ss2(line);
        //ss2 >> nfc3;
        //if (rank==0) printf(" %d FC3s.\n", nfc3);
        //mem->allocate(fc3,nfc3);
        while (getline(fh, line))
        {

            stringstream ss(line);
            ss >> i >> j >> mcc;
            if (counter % 10000 == 0){
                printf("  i: %d\n", i);
            }

            //if (lowmode <= i && i <= highmode){
                //if (abs(mcc)>1e24){
                if (abs(mcc)>0){
                    fprintf(fh_mcc2, "%d %d %.15e\n", i,j,mcc);
                }
            //}
            
            counter++;
            
        }

        fh.close();

        

    }
    fclose(fh_mcc2);

}



void Ifc2mcc::readEmat()
{

    ifstream readfile;

    if (rank==0) printf(" Reading Eigenvectors.\n");

    //int natoms = 8;

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

unsigned long long int Ifc2mcc::nChoosek( unsigned long long int n, unsigned long long int k )
{
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    unsigned long long int result = n;
    for( unsigned long long int i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}

/*
Calculate generalized velocities, similar to SMCC2s.
*/

void Ifc2mcc::go_gv(int alpha)
{
    order = 2;

    if (rank==0) printf(" Converting IFCs to GVs.\n");

    natoms = lmp->atom->natoms; // and natoms = 8
    if (rank==0) printf(" natoms: %d\n", natoms);

    /* Read FCs */
    readFcs();

    /* Read eigenvectors */
    readEmat();

    /* Read frequencies */
    if (rank==0) printf(" Reading frequencies.\n");
    double *freq;
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



    //printf(" %e\n", lmp->atom->x[1][0]);
    int *type = lmp->atom->type;
    int *tag = lmp->atom->tag;
    //printf(" %d\n", type[0]);
    double *mass = lmp->atom->mass;
    double **x = lmp->atom->x;

    /*
    mem->allocate(mass,lmp->atom->ntypes);
    // Convert mass to SI units
    for (int t=0;t<lmp->atom->ntypes; t++){
        printf(" %d\n", t);
        mass[t] = lmp->atom->mass[t]/(6.0221409e+23*1e3);
    }
    */

    //printf(" %d\n", lmp->atom->ntypes);
    //printf(" %e\n", mass[type[0]]);

    unsigned long long int n = natoms*3 - 3; // dont count translational modes
    //int nmcc2 = nChoosek(n+2-1,n-1);
    double mass_factor;

    int ii,jj,kk,ll; // dynamical matrix indices e.g. ii = 3*i+a

    int i,j;
    int a,b;
    double disp; // displacement between two atoms, used in the GV calculation. 
    //printf("ASDF\n");
    //printf("rank: %d\n", rank);
    if (order==2){
        //if (rank==0){ // Do MCC2 calculations on one proc, it's easier to have in a single file.

            FILE * fh_mcc2;
            //fh_mcc2 = fopen("MCC2","w");
            //double mass = 4.6637066e-26; // Assume every atom has silicon mass for now.
            if (rank==0) printf(" natoms: %d\n", natoms);
            //int nmcc2 = nChoosek(n+2-1,n-1);
            int nmcc2 = (3*natoms)*(3*natoms);
            //printf(" %d MCC2s.\n", nmcc2);

            char fh_mcc2_filename[64];
            sprintf (fh_mcc2_filename, "gv/GV_%d", rank);
            fh_mcc2 = fopen(fh_mcc2_filename,"w");
            //double mass1 = 4.6637066e-26; // Assume every atom has silicon mass for now.
            //natoms = 8; // and natoms = 8

            int n = natoms*3; // dont count translational modes

            // Split modes across procs

            int *nepp;
            mem->allocate(nepp, mc->nprocs); // number elements (MCC3s) per proc
            for (int p=0; p<mc->nprocs; p++){
                nepp[p] = n/mc->nprocs;
            }

            // divide up the remainder
            for (int p=0; p<(n % mc->nprocs); p++){
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

            if (rank==0){
                printf(" Splitting modes on procs like:\n");
                for (int p=0; p<mc->nprocs; p++){
                    printf("  %d Modes on proc %d.\n", nepp[p],p);
                }
            }
            //printf("  Proc %d: %d-%d.\n", rank,start_indx,end_indx);


            fprintf(fh_mcc2, "%d\n", nepp[rank]*n); // Number of MCC2s on this proc
            //fprintf(fh_mcc2, "%d\n", nmcc2);
            mass_factor = 1.0/(6.0221409e+23*1e3*6.0221409e+23*1e3); // convert 2 masses to kg in denominator

            double fc2_value;
            double mcc2_value;
            double K1, K2; // K_n1n2_AB and K_n2n1_BA
            double mi,mj;
            double phi;
            double term1,term2;
            int itype,jtype;

            for (int n1=start_indx; n1<end_indx; n1++){
            //for (int n1=10360; n1<10361; n1++){
                if (rank==0){
                    if ( (n1 % 1)== 0) printf(" n1: %d\n", n1);
                }
                for (int n2=0; n2<3*natoms; n2++){
                        //if (n1==n2){
                        mcc2_value = 0.0;
                        K1 = 0.0;
                        K2 = 0.0;
                        for (int f=0; f<nfc2; f++){ // Loop through nonzero FC3s
                            //printf("%d %d %d %d %f\n", i,j,a,b,phi);
                            i = fc2[f].i;
                            j = fc2[f].j;
                            a = fc2[f].a;
                            b = fc2[f].b;
                            phi = fc2[f].val*(2.179874099E-18)*1.89e+10*1.89e+10; // Convert from Ryd/Bohr^2 to J/m^2
                            //phi = fc2[f].val*13.605698066*(1./0.529177249)*(1./0.529177249); // Convert from Ryd/Bohr^2 to eV/A^2, THIS IS FOR COMPARING TO PAIR_TEP.CPP
                            //printf(" %d %d %d %d %e\n", i,a,j,b,fc2_value);
                            //printf(" %f\n", emat[a+i*3][n1]);
                            //printf(" %f\n", emat[b+j*3][n2]);
                            //printf(" %f\n", mass[type[i]]);
                            //printf(" %f\n", mass[type[j]]);
                            //printf("Got IFC.\n");
                            itype = type[i];
                            jtype = type[j];
                            mi = mass[type[i]]*(1.0/(6.02214076e23*1e3)); // atom mass in kg
                            mj = mass[type[j]]*(1.0/(6.02214076e23*1e3)); // atom mass in kg

                            //printf("Got mass.\n");
                            disp = x[i][alpha]-x[j][alpha];
                            if (std::abs(disp) > lmp->domain->boxhi[alpha]/2.0 && disp > 0) disp -= lmp->domain->boxhi[alpha];
                            else if (std::abs(disp) > lmp->domain->boxhi[alpha]/2.0 && disp < 0) disp += lmp->domain->boxhi[alpha];

                            // printf("Got disp.\n");
                            mcc2_value += disp*1e-10*( (phi*emat[3*i+a][n1]*emat[3*j+b][n2])/sqrt(mi*mj) ); // m/s

                            //printf("Added to MCC2.\n");

                            
                            //else if (type[i]==2 && type[j]
                            //printf("ASDF\n");

                            ii=3*i+a;
                            jj=3*j+b;

                            //fprintf(fh_debug, "n1, n2: %d %d\n", n1,n2);
                            //fprintf(fh_debug, "  i j a b: %d %d %d %d\n", i,j,a,b);
                            //fprintf(fh_debug, "  fc2: %e\n", fc2_value);
                            //fprintf(fh_debug, "  term: %e\n", (fc2_value*emat[a+i*3][n1]*emat[b+j*3][n2])/(sqrt(mass*mass)));
                        }

                        //if (n1==0 && n2==0){
                          //if (i != j){
                            //fprintf(fh_debug, "%d %d %e\n", n1,n2,mcc2_value);
                          //}
                        //}

                        //I used to divide by 2, but we can do this in the heat transfer expression instead.
                        //if (n1 != n2) mcc2_value = 0.5*mcc2_value;

                        mcc2_value = mcc2_value*0.5*(1.0/sqrt(2.0*3.14159265359*freq[n1]*1e12*2.0*3.14159265359*freq[n2]*1e12));
                        if (abs(mcc2_value) > 0.0){
                            fprintf(fh_mcc2, "%d %d %.20e\n", n1,n2,mcc2_value);
                        }
                    //}
                }
            }

            fclose(fh_mcc2);

        //} // if rank == 0
    } // if order == 2

    mem->deallocate(freq);

    if (rank==0){

        FILE * fh_mcc2;
        fh_mcc2 = fopen("GV","w");
        natoms = lmp->atom->natoms; 
        int nmodes = 3*natoms; //-3; // Number of modes.
        //int nmodes_i = highmode-lowmode+1; // Number of modes "i" we are concerned with in extract MCC3s for. 
        //fprintf(fh_mcc3, "%d\n", nmodes_i*nmodes*nmodes);
        
        for (int f=0; f<mc->nprocs; f++){

            printf(" File: %d\n", f);

            /* Read file and look for relevant mode. */

            char filename[64];
            sprintf (filename, "gv/GV_%d", f);
            //fh_debug = fopen(debug,"w");

            ifstream fh(filename);
            string line;

            int i,j;
            double mcc;
            int counter = 0;

            getline(fh,line); // Ignore first line.
            //stringstream ss2(line);
            //ss2 >> nfc3;
            //if (rank==0) printf(" %d FC3s.\n", nfc3);
            //mem->allocate(fc3,nfc3);
            while (getline(fh, line))
            {

                stringstream ss(line);
                ss >> i >> j >> mcc;
                if (counter % 10000 == 0){
                    printf("  i: %d\n", i);
                }

                //if (lowmode <= i && i <= highmode){
                    //if (abs(mcc)>1e5){
                        fprintf(fh_mcc2, "%d %d %.15e\n", i,j,mcc);
                    //}
                //}
                
                counter++;
                
            }

            fh.close();

            

        }
        fclose(fh_mcc2);
    }

}
