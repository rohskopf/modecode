#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <random>
#include "mpi.h"

#include "verify.h"
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

#include <ctime>

using namespace std;

using namespace LAMMPS_NS;

using namespace MC_NS;

Verify::Verify(MC *mc) : Ptrs(mc) {
    fh_debug = fopen("DEBUG_VERIFY","w");
}

Verify::~Verify() 
{

    fclose(fh_debug);
    mem->deallocate(mcc3);

};

/*
Read MCC file and store values
*/
void Verify::readMCC(const char* mcc3_filename)
{

    int natoms = lmp->atom->natoms;
    mem->allocate(mcc3,3*natoms,3*natoms,3*natoms);
    
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            for (int k=0; k<3; k++){
                mcc3[i][j][k] = 0.0; // translational modes have no interaction
            }
        }
    }  
    
    //char filename[64];
    //sprintf (filename, "%s", mcc3_filename);
    
    ifstream fh(mcc3_filename);     
    string line;
    int i,j,k;
    double value;
    int counter = 0;
    while (getline(fh, line))
    {
        stringstream ss(line);
        ss >> i >> j >> k >> value;
        if (i==j && j==k){
            mcc3[i][j][k]=value;
            counter++;
        }
        else if (i==j && i!=k){ //iij case
            mcc3[i][j][k]=value;
            mcc3[i][k][j]=value;
            mcc3[k][i][j]=value;
            counter += 3;
        }
        else if (i!=j && j==k){ //ijj case
            mcc3[i][j][k]=value;
            mcc3[j][i][k]=value;
            mcc3[j][k][i]=value;
            counter += 3;
        }
        else if (i!=j && j!=k){ //ijk case
            mcc3[i][j][k]=value;
            mcc3[i][k][j]=value;
            mcc3[j][i][k]=value;
            mcc3[j][k][i]=value;
            mcc3[k][i][j]=value;
            mcc3[k][j][i]=value;
            counter += 6;
        }

    }
    fh.close();
    //printf(" %d\n", counter);
}

/*
Compare MTEP and LAMMPS potential against random mode displacements.
Could compare total and anharmonic energy.
*/
void Verify::compareRandom(int nconfigs, double mag, int seed)
{

    printf(" Verifying %d configs with magnitude=%f.\n", nconfigs, mag);

    double **u;
    double *q; // mode amplitudes
    //printf(" %d atoms.\n", lmp->atom->natoms);
    mem->allocate(u,lmp->atom->natoms,3);
    nmodes = lmp->atom->natoms*3;
    mem->allocate(q,nmodes);
    double *mass;
    mem->allocate(mass,lmp->atom->ntypes);
    for (int t=1; t<=lmp->atom->ntypes; t++){
        mass[t] = lmp->atom->mass[t];
        //printf("mass[%d]: %f\n", t,lmp->atom->mass[t]);
    }
    int *type;
    mem->allocate(type,lmp->atom->natoms);
    for (int i=0; i<lmp->atom->natoms; i++){
        type[i]=lmp->atom->type[i];
    }

    std::default_random_engine engine(seed); // engine for generating random numbers
    std::uniform_real_distribution<double> dist(-1.0*mag,1.0*mag); // for choosing random [-mag,mag]
    //double rando = dist(engine);
    //printf(" rando: %f\n", rando);

    double mass_convert;

    double pe_lmp, pe_mtep;

    // First calculate cohesive energy for LAMMPS potential
    double e_pot_min = calcCohesiveEnergy();
    //printf(" mass[%d] = %f\n", type[0],mass[type[0]]);
    //printf("type[%d]: %d\n", 0,type[0]);
    double pi = 3.14159265359;

    // Loop through all congfigurations, generate random atom displacements from which mode displacements are calculated
    for (int m=0; m<nconfigs; m++){

        // Generate random displacements for this config
        for (int i=0; i<lmp->atom->natoms; i++){
            u[i][0] = dist(engine);
            u[i][1] = dist(engine);
            u[i][2] = dist(engine);
            //printf(" %f %f %f\n", u[i][0],u[i][1],u[i][2]);
            //printf(" mass[type[%d]] = %f\n", n,mass[type[n]]);
        }

        // Convert to mode displacmements
        for (int n=0; n<nmodes; n++){
            q[n]=0.0;
            //printf(" %d\n", n);
            for (int i=0; i<lmp->atom->natoms; i++){
                //printf(" %d\n", i);
                mass_convert = (mass[type[i]]/(6.02214076e23))*1e-3;
                //printf(" %d\n", i);
                //mass_convert = mass[type[i]];
                q[n] += sqrt(mass_convert)*mc->in->emat[0+(i*3)][n]*u[i][0]*1e-10;
                q[n] += sqrt(mass_convert)*mc->in->emat[1+(i*3)][n]*u[i][1]*1e-10;
                q[n] += sqrt(mass_convert)*mc->in->emat[2+(i*3)][n]*u[i][2]*1e-10;
            }
            //printf(" q[%d]: %e\n", n,q[n]);
        }

        // Convert back to atomic displacements to verify conversion
        /*
        for (int i=0; i<lmp->atom->natoms; i++){
            mass_convert = (mass[type[i]]/(6.02214076e23))*1e-3;
            u[i][0] = 0.0;
            u[i][1] = 0.0;
            u[i][2] = 0.0;
            for (int n=0; n<nmodes; n++){
                u[i][0] += (1.0/sqrt(mass_convert))*mc->in->emat[0+(i*3)][n]*q[n];
                u[i][1] += (1.0/sqrt(mass_convert))*mc->in->emat[1+(i*3)][n]*q[n];
                u[i][2] += (1.0/sqrt(mass_convert))*mc->in->emat[2+(i*3)][n]*q[n];
                //printf(" mass[type[%d]] = %f\n", n,mass[type[n]]);
            }
            printf(" %e %e %e\n", u[i][0]*1e10,u[i][1]*1e10,u[i][2]*1e10);
        }
        */

        // Calculate potential using LAMMPS
        pe_lmp = calcEnergyLmp(u);

        // Calculate potential using MTEP
        pe_mtep = calcEnergyMtep(q);

        printf(" %e %e\n", pe_lmp-e_pot_min,pe_mtep);
    }

    mem->deallocate(u);
    mem->deallocate(q);
    mem->deallocate(mass);
    mem->deallocate(type);
    

}

/*
Calcualte cohesive energy in LAMMPS, which is just energy of equilibrium config
*/
double Verify::calcCohesiveEnergy()
{

    double *x;
    double *e_ptr;
    double e;
    mem->allocate(x,lmp->atom->natoms*3);
    for (int i=0;i<lmp->atom->natoms;i++){
        //printf(" %f %f %f\n", u[i][0],u[i][1],u[i][2]); 
        x[0+i*3]=mc->in->x0[i][0];
        x[1+i*3]=mc->in->x0[i][1];
        x[2+i*3]=mc->in->x0[i][2];
    }

    lammps_scatter_atoms(lmp,"x",1,3,x); // Simply change positions if nbonds!=0.
    lmp->input->one("run 0");
    e_ptr = (double *) lammps_extract_compute(lmp, "P", 0, 0); // style = 0 for global data, type = 0 for scalar quantity
    e= e_ptr[0];
    printf(" Cohesive energy of LAMMPS potential: %f\n", e);

    mem->deallocate(x);
    return e;

}

/*
Calcualte potential energy using LAMMPS
Input: Nx3 array of atomic displacements
*/
double Verify::calcEnergyLmp(double **u)
{
    double *x;
    double *e_ptr;
    double e;
    mem->allocate(x,lmp->atom->natoms*3);
    for (int i=0;i<lmp->atom->natoms;i++){
        //printf(" %f %f %f\n", u[i][0],u[i][1],u[i][2]); 
        x[0+i*3]=mc->in->x0[i][0] + u[i][0];
        x[1+i*3]=mc->in->x0[i][1] + u[i][1];
        x[2+i*3]=mc->in->x0[i][2] + u[i][2];
    }

    lammps_scatter_atoms(lmp,"x",1,3,x); // Simply change positions if nbonds!=0.
    lmp->input->one("run 0");
    e_ptr = (double *) lammps_extract_compute(lmp, "P", 0, 0); // style = 0 for global data, type = 0 for scalar quantity
    e= e_ptr[0];
    //printf(" %f\n", e);

    mem->deallocate(x);
    return e;
}

/*
Calcualte potential energy using MTEP
Input: 3N array of mode amplitudes
*/
double Verify::calcEnergyMtep(double *q)
{

    double e =0.0;
    double e_anh=0.0;

    double pi = 3.14159265359;

    for (int i=3; i<nmodes; i++){
        //printf(" q[%d]: %e\n",i,q[i]);
        for (int j=3; j<nmodes; j++){
            if (i==j){
                printf(" q[%d]: %e\n", i,q[i]);
                e += 0.5*4*pi*pi*mc->in->freqs[i]*mc->in->freqs[i]*q[i]*q[i]*6.242e+18; // Need to convert harmonic part to eV
            }
            for (int k=3; k<nmodes; k++){
                e += (1.0/6.0)*mcc3[i][j][k]*q[i]*q[j]*q[k]*1e30; // Need to convert mode coordinates back to Angstrom since MCC3s are in eV/A^3
                //printf("%e\n", (1.0/6.0)*mcc3[i][j][k]*q[i]*q[j]*q[k]*6.242e+18);
                e_anh += (1.0/6.0)*mcc3[i][j][k]*q[i]*q[j]*q[k]*1e30;
            }
        }
    }
    //e = 0.5*e;

    // Convert to eV
    //e = e*6.242e+18;
    //e_anh = e_anh*6.242e+18;
    printf(" e_harm: %e\n", e-e_anh);
    printf(" e_anha: %e\n", e_anh);
    return e;
}
