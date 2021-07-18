/*
 potential.cpp

 Copyright (c) 2018 Andrew Rohskopf

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include "mpi.h"
#include <math.h>       /* sqrt */
#include <random>

#include "potential.h"
#include "memory.h"
#include "input.h"
#include "compute.h"
#include "mode.h"

using namespace std;

using namespace EM3_NS;

Potential::Potential(EM3 *em3) : Pointers(em3) {


    fh_debug = fopen("DEBUG_potential", "w");
}

Potential::~Potential() 
{

    fclose(fh_debug);

};

void Potential::calculate()
{

    
    //printf(" Calculating potential.\n");

    double **xa = mode->xa; // If you change these pointers, it changes the original mode->variable!
    double **xa0 = input->xa0;
    double **va = mode->va;
    double **fa = mode->fa;
    double **fa2 = mode->fa2;
    double **fa3 = mode->fa3;
    double **fa4 = mode->fa4;

    double *xm = mode->xm; // If you change these pointers, it changes the original mode->variable!
    double *fm = mode->fm;
    double *fm2 = mode->fm2;
    double *fm3 = mode->fm3;
    double *fm4 = mode->fm4;

    int natoms = input->natoms;
    int nmodes = 3*natoms;

    // Zero potential and force intitially

    pea=0.0;
    for (int i=0; i<natoms; i++){
        for (int a=0;a<3;a++){
            fa[i][a]=0.0;
            fa2[i][a]=0.0;
            fa3[i][a]=0.0;
            fa4[i][a]=0.0;
        }
    }

    // Calculate potential
    int i,j,k,l;
    int a,b,c,d;
    double fc;

    /*
    for (int w=0; w<input->nfc2; w++){
        i = input->fc2[w].i;
        a = input->fc2[w].a;
        //fprintf(fh_debug, "%d\n", a);
        //a=2;
        f[i][a] = 1.0;
    }
    */
    
    double uijb;
    int ii,jj;
    if (input->space==0 || input->space==2){
        if (input->order >=2){
            for (int w=0; w<input->nfc2; w++){
                i=input->fc2[w].i;
                j=input->fc2[w].j;
                a=input->fc2[w].a;
                b=input->fc2[w].b;
                fc=input->fc2[w].val*(2.179874099E-18)*1.89e+10*1.89e+10;
                uijb = (xa[i][b]-xa0[i][b]) - (xa[j][b]-xa0[j][b]);
                pea += fc*(xa[i][a]-xa0[i][a])*(xa[j][b]-xa0[j][b]);
                fa[i][a] -= fc*(xa[j][b]-xa0[j][b]);
                fa2[i][a] -= fc*(xa[j][b]-xa0[j][b]);
                //f[i][a] += 0.5*fc*uijb;

                ii=3*i+a;
                jj=3*j+b;  
                
                
            }
        }

        if (input->order>=3){
            for (int w=0; w<input->nfc3; w++){

                i=input->fc3[w].i;
                j=input->fc3[w].j;
                k=input->fc3[w].k;

                a=input->fc3[w].a;
                b=input->fc3[w].b;
                c=input->fc3[w].c;

                fc=input->fc3[w].val*(2.179874099E-18)*1.89e+10*1.89e+10*1.89e+10;
                //pe += fc*(x[i][a]-x0[i][a])*(x[j][b]-x0[j][b])*(x[k][c]-x0[k][c]);
                fa[i][a] -= (1.0/2.0)*fc*(xa[j][b]-xa0[j][b])*(xa[k][c]-xa0[k][c]);
                fa3[i][a] -= (1.0/2.0)*fc*(xa[j][b]-xa0[j][b])*(xa[k][c]-xa0[k][c]);
            }
        }
        
        if (input->order>=4){
            for (int w=0; w<input->nfc4; w++){

                i=input->fc4[w].i;
                j=input->fc4[w].j;
                k=input->fc4[w].k;
                l=input->fc4[w].l;

                a=input->fc4[w].a;
                b=input->fc4[w].b;
                c=input->fc4[w].c;
                d=input->fc4[w].d;

                fc=input->fc4[w].val*(2.179874099E-18)*1.89e+10*1.89e+10*1.89e+10*1.89e+10;
                //pe += fc*(x[i][a]-x0[i][a])*(x[j][b]-x0[j][b])*(x[k][c]-x0[k][c])*(x[l][d]-x0[l][d]);
                fa[i][a] -= (1.0/6.0)*fc*(xa[j][b]-xa0[j][b])*(xa[k][c]-xa0[k][c])*(xa[l][d]-xa0[l][d]);
                fa4[i][a] -= (1.0/6.0)*fc*(xa[j][b]-xa0[j][b])*(xa[k][c]-xa0[k][c])*(xa[l][d]-xa0[l][d]);
            }
        }
    }
    

    /* Calculate mode forces */

    if (input->space==1 || input->space==2){
        //mode->atom2Mode(); // Convert atom coordinates to mode coordinates
        
        pem = 0.0;
        for (int i=0; i<nmodes; i++){
            fm[i] = 0.0;
            fm2[i] = 0.0;
            fm3[i] = 0.0;
            fm4[i] = 0.0;
        }

        
        double mcc;
        if (input->order>=2){
            for (int w=0; w<input->nmcc2; w++){
                i = input->mcc2[w].i;
                mcc = input->mcc2[w].val;
                pem += mcc*xm[i]*xm[i];
                fm[i] -= mcc*xm[i];
                fm2[i] -= mcc*xm[i];
                //printf(" %d %e %e\n", i, q[i], mcc);
            }
        }
        
        if (input->order>=3){
            for (int w=0; w<input->nmcc3; w++){
                i=input->mcc3[w].i;
                j=input->mcc3[w].j;
                k=input->mcc3[w].k;
                mcc = input->mcc3[w].val;
                fm[i] -= (1.0/2.0)*mcc*xm[j]*xm[k];
                fm3[i] -= (1.0/2.0)*mcc*xm[j]*xm[k];
                
                /*
                if (i==j && j==k){ // iii
                    f[i] -= input->mcc3[w].val*q[i]*q[i];
                    f3[i] -= input->mcc3[w].val*q[i]*q[i];
                }
                
                if (i==j && i!=k){ //iij
                    f[i] -= input->mcc3[w].val*q[j]*q[k];
                    f[k] -= input->mcc3[w].val*q[i]*q[j];
                    f3[i] -= input->mcc3[w].val*q[j]*q[k];
                    f3[k] -= input->mcc3[w].val*q[i]*q[j];
                }
                
                if (i!=j && j==k){ // ijj
                    f[i] -= input->mcc3[w].val*q[j]*q[k];
                    f[j] -= input->mcc3[w].val*q[i]*q[k];
                    f3[i] -= input->mcc3[w].val*q[j]*q[k];
                    f3[j] -= input->mcc3[w].val*q[i]*q[k];
                }    
                
                if (i!=j && j!=k){ // ijk
                    f[i] -= input->mcc3[w].val*q[j]*q[k];
                    f[j] -= input->mcc3[w].val*q[i]*q[k];
                    f[k] -= input->mcc3[w].val*q[i]*q[j];
                    f3[i] -= input->mcc3[w].val*q[j]*q[k];
                    f3[j] -= input->mcc3[w].val*q[i]*q[k];
                    f3[k] -= input->mcc3[w].val*q[i]*q[j];
                }
                */
                //f[i] -= input->mcc3[f].val*q[
            }
        }

        if (input->order>=4){
            for (int w=0; w<input->nmcc4; w++){
                i=input->mcc4[w].i;
                j=input->mcc4[w].j;
                k=input->mcc4[w].k;
                l=input->mcc4[w].l;
                mcc = input->mcc4[w].val;
                fm[i] -= (1.0/6.0)*mcc*xm[j]*xm[k]*xm[l];
                fm4[i] -= (1.0/6.0)*mcc*xm[j]*xm[k]*xm[l];
                //printf(" %e\n", mcc);
            }
        }

        //mode->mode2Atom(); // convert mode quantities to atom quantities
    }

}
