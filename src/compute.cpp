/*
This class computes quantities based on IFCs, eigenvectors, MCCs, etc.
Current functions include:
calcPR: Calculate participation ratio based on eigenvectors.
*/

#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <random>
#include "mpi.h"

#include "compute.h"
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

Compute::Compute(MC *mc) : Ptrs(mc) {
    //fh_debug = fopen("DEBUG_ASR","w");
    //fh_fc2 = fopen("FC2_ASR","w");
    pr_call=false;
    esp_call=false;
    
    pi = 3.1415926535;
    hbar = 1.0545718e-34;
    kb = 1.38064852e-23;
}

Compute::~Compute() 
{

    //fclose(fh_debug);
    //fclose(fh_fc2);

    //if (order >= 3) mem->deallocate(fc3);
    //if (order >= 4) mem->deallocate(fc4);

    if (pr_call){
        mem->deallocate(emat);
        mem->deallocate(prs);
    }

    if (esp_call){
        mem->deallocate(emat);
        mem->deallocate(esp);
    }
    

};

/*
Compute mode participiation ratios.
*/

void Compute::participationRatio()
{

    printf(" Computing participation ratios for %d atom system.\n", natoms); 

    mem->allocate(prs,natoms*3);
    FILE * fh_pr;
    fh_pr = fopen("pr.dat","w");

    // First check unit magnitude of eigenvectors.
    double sum_col;
    for (int n=0; n<3*natoms; n++){
        sum_col = 0.0;
        for (int i=0; i<natoms; i++){
            for (int a=0; a<3; a++){
                sum_col += emat[3*i+a][n]*emat[3*i+a][n];
            }
        }
        sum_col = round( sum_col * 10000.0 ) / 10000.0;
        if (sum_col != 1.0) printf(" Mode %d doesn't satisfy column unit magnitude.\n");
    }
    double sum_row;
    for (int i=0; i<natoms; i++){
        for (int a=0; a<3; a++){
            sum_row = 0.0;
            for (int n=0; n<3*natoms; n++){
                sum_row += emat[3*i+a][n]*emat[3*i+a][n];
            }
            sum_row = round( sum_row * 10000.0 ) / 10000.0;
            //printf(" sum_row: %f\n", sum_row);
            if (sum_row != 1.0) printf(" Mode %d doesn't satisfy row unit magnitude.\n");
        }
    }

    // Loop over all modes to calculate PR for each.
    double eia; // "a" component of eigenvector on atom "i", for a particular mode. 
    double edotin; // Dot product of eigenvector on atom i, in mode n, with itself.
    double edotin_sq;
    double sum_edotin;
    double sum_edotin_sq; // Sum of squared edotin.
    double numerator;
    double denominator;
    for (int n=0; n<3*natoms; n++){

        sum_edotin = 0.0;
        sum_edotin_sq = 0.0;

        for (int i=0; i<natoms; i++){
            edotin = 0.0;
            for (int a=0; a<3; a++){
                eia = emat[3*i+a][n];
                edotin+=(eia*eia);
            }
            edotin_sq = edotin*edotin;
            sum_edotin += edotin;
            sum_edotin_sq += edotin_sq;
        }

        prs[n] = (sum_edotin*sum_edotin)/(natoms*sum_edotin_sq);

        fprintf(fh_pr,"%d %f\n", n+1, prs[n]);
    }

    fclose(fh_pr);

}

/*
Compute eigenvector spatial parameters.
*/

void Compute::eigenvectorSpatialParameter()
{
  int *type = lmp->atom->type;
  natoms = lmp->atom->natoms;
  readEmat();
  mem->allocate(esp,natoms*3);


  FILE * fh;
  fh = fopen("esp.dat","w");


   
  double total;
  for (int n=0; n<3*natoms; n++){
    printf(" %d\n", n);
    esp[n] = 0.0;
    total = 0.0;
    
    for (int i=0; i<natoms; i++){
      total += sqrt(emat[3*i+0][n]*emat[3*i+0][n] + emat[3*i+1][n]*emat[3*i+1][n] + emat[3*i+2][n]*emat[3*i+2][n]);
      if (type[i]==2){
        esp[n] += sqrt(emat[3*i+0][n]*emat[3*i+0][n] + emat[3*i+1][n]*emat[3*i+1][n] + emat[3*i+2][n]*emat[3*i+2][n]);
      }
    }
    
    printf(" Writing.\n");
 
    fprintf(fh, "%d %.6f\n", n,esp[n]/total);
  }
  

  fclose(fh);
  
}

/*
Read eigenvector matrix. 
To index "a" Cartesian component of eigenvector on atom "i" in mode "n1", do:
    emat[a+i*3][n1]
*/

void Compute::readEmat()
{
    printf(" Reading eigenvector matrix for %d atoms.\n",natoms);

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


void Compute::readFcs()
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

/*
Sum directional mode coupling constants.
*/
void Compute::sumDMCC()
{

    printf(" Reading GVs.\n");
    natoms = lmp->atom->natoms;
    mem->allocate(kn,natoms*3);
    
    for (int n=0; n<3*natoms; n++){
        kn[n] = 0.0;
    }
    

    // Read number of GVs
    int n1,n2;
    double val;
    ifstream fh("../GV");
    string line;
    //ngv = 0;
    while (getline(fh, line))
    {
        stringstream ss(line);
        ss >> n1 >> n2 >> val;
        kn[n1] += abs(val);
        //if (n1==n_indx) ngv++;
    }
    //ngv = ngv-1; // subtract 1 for first line
    //printf(" Found %d FC3s.\n", nfc3);
    fh.close();

    FILE * fh_kn;
    fh_kn = fopen("kn.dat","w");
    
    for (int n=3; n<3*natoms; n++){
        fprintf(fh_kn, "%d %e\n", n, kn[n]);
    }
    
    fclose(fh_kn);

}

void Compute::linewidths(int n1, double temperature)
{
    printf(" n=%d\n", n1);
    printf(" Reading MCC3.\n");
    
    // declare arrays
    struct mcc3_struct{
        int i,j,k;
        double val;
    };
    mcc3_struct *mcc3; // array of MCC3 values
    double *freq;
    struct lw_struct{
        int i,j,k;
        double val;
    };
    lw_struct *lw; // array of linewidths values
    
    /* MCC3s */
    
    ifstream fh_mcc3("MCC3");
    string line;
    int nmcc3 = 0;
    while (getline(fh_mcc3, line))
    {
        nmcc3=nmcc3+1;
    }

    printf("  Found %d MCC3s.\n", nmcc3);

    mem->allocate(mcc3,nmcc3);

    fh_mcc3.close();

    ifstream readfile6;
    readfile6.open("MCC3");

    for (int n=0; n<nmcc3; n++){
      readfile6 >> mcc3[n].i >> mcc3[n].j >> mcc3[n].k >> mcc3[n].val;
    }

    readfile6.close();
    
    printf(" Reading FREQUENCIES\n");
    
    /* FREQUENCIES */
    
    ifstream fh_freq("FREQUENCIES");
    int nmodes = 0;
    while (getline(fh_freq, line))
    {
        nmodes=nmodes+1;
    }

    printf("  Found %d frequencies.\n", nmodes);

    mem->allocate(freq,nmodes);

    ifstream readfile3;

    for (int i = 0; i < nmodes; i++) {
      freq[i]=0.0;
    }

    readfile3.open("FREQUENCIES");
    //readfile2.open("ev_real.txt");

    if (!readfile3.is_open()) {
        printf("Unable to open FREQUENCIES.\n");
        exit(1);
    }

    //printf("natoms: %d\n",  natoms);
    for (int i=0;i<nmodes;i++){
      readfile3>>freq[i];
    }
    
    /* Start linewidth calculation */
    mem->allocate(lw, nmcc3);
    
    double summ;
    double delta1,delta2;
    //double mcc3;
    double dist2,dist3; // distribution functions of modes n2 and n3.
    double term1,term2;
    
    // Loop over MCC3s
    int n2,n3;
    for (int s=0; s<nmcc3; s++){
        n2 = mcc3[s].j;
        n3 = mcc3[s].k;
        if ( (abs(freq[n1]-freq[n2]-freq[n3]) < 1e-1) )
        {
            delta1 = 1.0/(1e-1*1e12); // Need to divide by 2pi?
            printf("%d %d %d\n", s, n2,n3);
        } 
        else{
            delta1 = 0.0;
        }
        
        if ( (abs(freq[n1]+freq[n2]-freq[n3]) < 1e-1) )
        {
            delta2 = 1.0/(1e-1*1e12);
            printf("%d %d %d\n", s,n2,n3);
        }
        else{
            delta2 = 0.0;
        }

        // Calculate MCC3 for n1,n2,n3
        // if (rank==0) printf("  Calculate MCC3 for %d-%d-%d: %f %f %f\n",n1,n2,n3,freq[n1],freq[n2],freq[n3]);
        //mcc3 = calcMCC3(n1,n2,n3);
        dist2 = calcDistribution(n2,temperature,freq[n2]);
        dist3 = calcDistribution(n3,temperature,freq[n3]);
        //printf(" %e %e %e\n", mcc3, dist2, dist3);
        //if (abs(mcc3) > 1e47) printf(" Yup!\n");
        term1 = 0.5*(1+dist2+dist3)*delta1;
        term2 = (dist2-dist3)*delta2;
        //if (rank==0) printf("   %e\n", mcc3);
        //printf("%e\n", (term1+term2) );
        //printf("  %e\n", mcc3);
        //printf("  %e %e\n", dist2, dist3);
        //printf("  %e\n", ( (mcc3*mcc3)/(2.0*pi*freq[n2]*1e12*2.0*pi*freq[n3]*1e12) )*(term1+term2) );
        lw[s].j = n2;
        lw[s].k = n3;
        lw[s].val =  ( (mcc3[s].val*mcc3[s].val)/(2.0*pi*freq[n2]*1e12*2.0*pi*freq[n3]*1e12) )*(term1+term2);
        lw[s].val = ( (pi*hbar*hbar)/(8.0*2.0*pi*freq[n1]*1e12))*lw[s].val; // Units of J.
        // Linewidths must have units of rad/s (same as freq).
        // Convert from J to eV to THz.
        //linewidths_p[n1] = linewidths_p[n1]*6.242e+18*241.8; // THz
        // Convert from J to eV.
        lw[s].val = lw[s].val*6.242e+18; // eV
        // Convert from eV to rad/s, using E=hbar*omega -> omega = E/hbar.
        //lw[s].val = (lw[s].val/(6.582119569e-16)); // rad/s
        //printf("   %e %e\n", linewidths[n1], summ);
    
        // Use rad/s or Hz units for frequencies?
        //summ += ( (mcc3*mcc3)/(1.0*freq[n2]*1e12*1.0*freq[n3]*1e12) )*(term1+term2);
        
    }
    
    // Print linewidths
    FILE * fh_lw;
    fh_lw = fopen("lw.dat","w");
    for (int s=0; s<nmcc3; s++){
        fprintf(fh_lw, "%.6e %.6e %.6e\n", freq[lw[s].j], freq[lw[s].k], lw[s].val);
    
    }
    
    fclose(fh_lw);
    
    // deallocate
    mem->deallocate(mcc3);
    mem->deallocate(freq);
    mem->deallocate(lw);
}

/*
Calculate BE distribution for a particular mode n.
*/
double Compute::calcDistribution(int n, double temperature, double freq_n)
{
    double dist;
    int quantum_setting = 0;
    if (quantum_setting==0){ // Classical MB distribution.
        dist = exp(-1.0*(1.0/(kb*temperature))*(hbar*2.0*freq_n*1e12) );
    }
    if (quantum_setting==1){ // Quantum BE distribution
        double exponent = exp(1.0*(1.0/(kb*temperature))*(hbar*2.0*freq_n*1e12) );
        dist = 1.0/(exponent-1);
    }

    return dist;
}
