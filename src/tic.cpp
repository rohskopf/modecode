/*
This class computes thermal interface conductance (TIC), from SMCC2s and MCC3s.
*/

#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <random>
#include "mpi.h"

#include "tic.h"
#include "mem.h"
#include "in.h"
//#include "config.h"

// LAMMPS include files
#include "lammps.h"
#include "input.h"
#include "atom.h"
#include "library.h"
#include "memory.h"
#include "domain.h"

#include <ctime>

using namespace std;

using namespace MC_NS;

Tic::Tic(MC *mc) : Ptrs(mc) {
    //fh_debug = fopen("DEBUG_ASR","w");
    //fh_fc2 = fopen("FC2_ASR","w");

    rank = mc->rank;

    kb = 1.38064852e-23;
    hbar = 1.054571817e-34;
    pi = 3.14159265359;

}

Tic::~Tic() 
{

    mem->deallocate(linewidths);
    mem->deallocate(linewidths_p);
    mem->deallocate(freq);
    //mem->deallocate(nepp);

};

/*
Main function for computing TIC. 
*/

void Tic::go()
{

    natoms = lmp->atom->natoms; // and natoms = 8
    nmodes = 3*natoms;
    type = lmp->atom->type;
    tag = lmp->atom->tag;
    //printf(" %d\n", type[0]);
    mass = lmp->atom->mass;
    double boxhi_x,boxhi_y,boxhi_z;
    boxhi_x = lmp->domain->boxhi[0];
    boxhi_y = lmp->domain->boxhi[1];
    boxhi_z = lmp->domain->boxhi[2];
    // Cross-sectional area.
    double ac = boxhi_x*boxhi_y;
    ac = ac*1e-10*1e-10;
    if (rank==0) printf(" Cross-section area: %e A\n", ac*1e10*1e10);


    if (rank==0) printf(" Computing TIC for %d atom system.\n", natoms); 

    mem->allocate(freq,nmodes);
    mem->allocate(linewidths, nmodes);
    mem->allocate(linewidths_p, nmodes);
    for (int n=0; n<nmodes; n++){
        linewidths_p[n] = 0.0; //1e-15; // Default value gives relaxation time of 500 ps. 
    }

    // Read EMAT, FC3s, and frequencies since these are inputs.
    readEmat();
    readFcs();
    readFrequencies();
    readSMCC2();

    // Compute linewidths for all modes.
    if (rank==0) printf("Calculating linewidths.\n");
    calcLineWidths();
    MPI_Allreduce(linewidths_p,linewidths,nmodes,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    // Loop over all SMCC2s to calculate TIC.
    
    if (rank==0){

        // Print linewidths (units of THz).
        fh_linewidths = fopen("linewidths.dat", "w");
        for (int n=0; n<nmodes; n++){
            fprintf(fh_linewidths, "%e %e\n", freq[n],linewidths[n]);
        }
        fclose(fh_linewidths);

        double conductance;
        int n1,n2;
        double lw1,lw2,freq1,freq2; // linewidths and frequencies for a pair.
        double tau;
        double summ = 0.0;
        double coeff;
        for (int s=0; s<nsmcc2; s++){
            //printf("s: %d\n", s);
            n1 = smcc2[s].n1;
            n2 = smcc2[s].n2;
            // Convert linewidths from THz to rad/s, just like frequency.
            lw1 = 2.0*pi*linewidths[n1]*1e12;
            lw2 = 2.0*pi*linewidths[n2]*1e12;
            freq1 = 2.0*pi*freq[n1]*1e12;  
            freq2 = 2.0*pi*freq[n2]*1e12;
            tau = (lw1+lw2)/( ((lw1+lw2)*(lw1+lw2)) + ( (freq1-freq2)*(freq1-freq2) ) );
            if (freq1 > 0.0 && freq2 > 0.0){
                tau = (lw1+lw2)/( ((lw1+lw2)*(lw1+lw2)) + ( (freq1-freq2)*(freq1-freq2) ) );
                coeff = (smcc2[s].val*smcc2[s].val)/(freq1*freq2);
                //printf("   %e %e %e\n", coeff,tau,coeff*tau);
                summ = summ + (coeff*tau);
            }
        }
        conductance = (kb/(8.0*ac))*summ;
        //printf("summ: %e\n", summ);
        printf("Conductance: %e\n", conductance);
    }
    


}

void Tic::calcLineWidths()
{

    int *nepp;
    mem->allocate(nepp, mc->nprocs); // number elements (MCC3s) per proc
    for (int p=0; p<mc->nprocs; p++){
        nepp[p] = nmodes/mc->nprocs;
    }

    // divide up the remainder
    for (int p=0; p<(nmodes % mc->nprocs); p++){
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
    printf("  Proc %d: %d-%d.\n", rank,start_indx,end_indx);


    //fprintf(fh_mcc3, "%d\n", nepp[rank]*n*n); // Number of MCC3s on this proc

    double summ;
    double delta1,delta2;
    double mcc3;
    double dist2,dist3; // distribution functions of modes n2 and n3.
    double term1,term2;
    //if (rank==0) printf(" nepp[end-1] start end: %d %d %d\n", nepp[mc->nprocs-1],start_indx, end_indx);
    for (int n1=start_indx; n1<end_indx; n1++){
    //for (int n1=1000; n1<1001; n1++){
        if (rank==0) printf(" n1: %d\n", n1);
        summ = 0.0;
        for (int n2=3; n2<nmodes; n2++){
            if (rank==0 && (n2%5==0) ) printf("    %d\n", n2);
            for (int n3=3; n3<nmodes; n3++){
                delta1 = 0.0;
                delta2 = 0.0;
                if ( (abs(freq[n1]-freq[n2]-freq[n3]) < 1e-1) )
                {
                    delta1 = 1.0/(1e-1*1e12); // Need to divide by 2pi?
                } 
                
                if ( (abs(freq[n1]+freq[n2]-freq[n3]) < 1e-1) )
                {
                    delta2 = 1.0/(1e-1*1e12);
                }
                if ( (delta1>0.0) || (delta2>0.0) ){
                    // Calculate MCC3 for n1,n2,n3
                    // if (rank==0) printf("  Calculate MCC3 for %d-%d-%d: %f %f %f\n",n1,n2,n3,freq[n1],freq[n2],freq[n3]);
                    mcc3 = calcMCC3(n1,n2,n3);
                    dist2 = calcDistribution(n2);
                    dist3 = calcDistribution(n3);
                    //printf(" %e %e %e\n", mcc3, dist2, dist3);
                    //if (abs(mcc3) > 1e47) printf(" Yup!\n");
                    term1 = 0.5*(1+dist2+dist3)*delta1;
                    term2 = (dist2-dist3)*delta2;
                    //if (rank==0) printf("   %e\n", mcc3);
                    //printf("  %e\n", (term1+term2) );
                    //printf("  %e\n", ( (mcc3*mcc3)/(2.0*pi*freq[n2]*1e12*2.0*pi*freq[n3]*1e12) )*(term1+term2) );
                    summ += ( (mcc3*mcc3)/(2.0*pi*freq[n2]*1e12*2.0*pi*freq[n3]*1e12) )*(term1+term2);
                }

            }
        }
        linewidths_p[n1] = ( (pi*hbar*hbar)/(8.0*2.0*pi*freq[n1]*1e12))*summ; // Units of J.

        // Linewidths must have units of rad/s (same as freq).
        // Convert from J to eV to THz.
        linewidths_p[n1] = linewidths_p[n1]*6.242e+18*241.8; // THz
        //printf("   %e %e %e\n", (pi*hbar*hbar)/(8.0*2.0*pi*freq[n1]*1e12), linewidths[n1], summ);
    }
    
    mem->deallocate(nepp);
}

/*
Calculate 3rd order mode coupling constant (MCC3) for a particular triplet of modes.
*/
double Tic::calcMCC3(int n1, int n2, int n3)
{

    double mcc3 = 0.0;
    int i,j,k,a,b,c;
    double fc;
    int ii,jj,kk;
    double mi,mj,mk;
    double mass_factor = 1.0/(6.0221409e+23*1e3*6.0221409e+23*1e3*6.0221409e+23*1e3); // convert 3 masses to kg in denominator
    //printf(" Looping over %d FC3s.\n", nfc3);

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

        //sqrtmasses = sqrt(mass[type[i]]*mass[type[j]]*mass[type[k]]*mass_factor);
        //printf("%d %d %d\n",ii,jj,kk);
        //fprintf(fh_debug, "%d %d %d %d %d %d\n",i,a,j,b,k,c);
        //printf("About to calculate termu %d %d %d.\n",ii,jj,kk);

        mcc3 += (fc*emat[3*i+a][n1]*emat[3*j+b][n2]*emat[3*k+c][n3])/sqrt(mi*mj*mk); // (J/(m^3 * kg^(3/2))

    }    
    
    return mcc3;
}

/*
Calculate distribution for a particular mode n.
*/
double Tic::calcDistribution(int n)
{
    double temperature = 5; // K
    double dist = exp(-1.0*(1.0/(kb*temperature))*(hbar*2.0*freq[n]*1e12) );
    return dist;
}

/*
Read eigenvector matrix. 
To index "a" Cartesian component of eigenvector on atom "i" in mode "n1", do:
    emat[a+i*3][n1]
*/

void Tic::readEmat()
{
    //printf(" Reading eigenvector matrix for %d atoms.\n",natoms);

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


void Tic::readFcs()
{

    /* Read number of FC3s */
    ifstream fh("FC3");
    string line;
    nfc3 = 0;
    while (getline(fh, line))
    {
        nfc3++;
    }
    nfc3 = nfc3-1; // subtract 1 for first line
    //printf(" Found %d FC3s.\n", nfc3);
    fh.close();

    mem->allocate(fc3,nfc3);

    /* Read and store FC3s */
    int i,j,k;
    int a,b,c;
    double val;
    int counter = 0;
    fh.open("FC3");
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

/*
Read FREQUENCIES file.
*/
void Tic::readFrequencies()
{
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
Read SMCC2s.
*/
void Tic::readSMCC2()
{

    /* Read number of FC3s */
    ifstream fh("SMCC2");
    string line;
    nsmcc2 = 0;
    while (getline(fh, line))
    {
        nsmcc2++;
    }
    nsmcc2 = nsmcc2-1; // subtract 1 for first line
    //printf(" Found %d FC3s.\n", nfc3);
    fh.close();

    mem->allocate(smcc2,nsmcc2);

    /* Read and store FC3s */
    int n1,n2;
    double val;
    int counter = 0;
    fh.open("SMCC2");
    getline(fh,line);
    while (getline(fh, line))
    {

        stringstream ss(line);
        ss >> n1 >> n2 >> val;
        
        smcc2[counter].n1=n1;
        smcc2[counter].n2=n2;
        smcc2[counter].val=val;
        //if (counter==0) printf(" %f\n", val);
        
        //printf("%f\n", fc3[counter].val);
        //printf(" %d\n", counter);
        counter++;
        
    }

    fh.close();

}
