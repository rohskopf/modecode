/*
 input.cpp

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

#include "input.h"
#include "memory.h"
#include "compute.h"

using namespace std;

using namespace EM3_NS;

Input::Input(EM3 *em3) : Pointers(em3) {
    fh_debug = fopen("DEBUG_input", "w");
}

Input::~Input() 
{
    //em3->memory->deallocate(xa);
    //printf("ASDF\n");
    //em3->memory->deallocate(va);
    //printf("ASDF\n");
    em3->memory->deallocate(xa0);
    //printf("ASDF\n");

    memory->deallocate(vm);
    memory->deallocate(xm);
    memory->deallocate(xa);
    memory->deallocate(va);
    memory->deallocate(mass);
    memory->deallocate(type);

    //memory->deallocate(freqs);
    if (space==1 || space==2){
        memory->deallocate(emat);
        if (order >= 2) memory->deallocate(mcc2);
        if (order >= 3) memory->deallocate(mcc3);
        if (order >= 4) memory->deallocate(mcc4);
    }

    if (space==0 || space==2){
        if (order >= 2) memory->deallocate(fc2);
        if (order >= 3) memory->deallocate(fc3);
        if (order >= 4) memory->deallocate(fc4);
    }

    fclose(fh_debug);

};

void Input::readinput()
{
    /* Read INPUT file */

    string line;

    // Declare scalar inputs
    double value;

    // Open INPUT file
    ifstream INPUT("INPUT");
    // Ignore the first line
    getline(INPUT, line); 
    string characters;
    // Get input variables
    int i=1;
    //for (int i=1; i<=12; i++)
    //{
        //getline(INPUT, line);
    while (getline(INPUT,line)){
        switch (i)
        {
            case 1:{ 
                stringstream ss(line);
                ss >> characters >> nsteps;
            }
            case 2:{ 
                stringstream ss(line);
                ss >> characters >> space;
            }
            case 3:{ 
                stringstream ss(line);
                ss >> characters >> nout;
            }
            case 4:{ 
                stringstream ss(line);
                ss >> characters >> order;
            }
            case 5:{ 
                stringstream ss(line);
                ss >> characters >> temp;
            }
            case 6:{ 
                stringstream ss(line);
                ss >> characters >> dt;
            }
            case 7:{
                stringstream ss(line);
                ss >> characters >> epsilon;
            }
            case 8:{
                stringstream ss(line);
                ss >> characters >> sigma;
            }

  
        } // switch (i)
        i++;
    } // for (int i=1..)

    INPUT.close();

}

void Input::readconfig()
{
    /* Read CONFIG file */

    ifstream config("CONFIG");
    string line;

    getline(config, line);
    stringstream ss(line);
    ss >> natoms >> ntypes;
    printf(" %d atom types.\n", ntypes);
    memory->allocate(mass,ntypes);

    getline(config, line);
    stringstream ss4(line);
    for (int t=0; t<ntypes; t++){
        ss4 >> mass[t];
        //printf("%e\n", mass[t]);
    }
    // Convert mass to SI units
    for (int t=0; t<ntypes;t++){
        mass[t] = mass[t]*(1.0/6.0221409e+23)*1e-3;
    }

    getline(config, line);
    stringstream ss2(line);
    ss2 >> box[0] >> box[1] >> box[2];
    //fprintf(fh_debug, "Box: %f %f %f\n", box[0], box[1], box[2]);

    printf(" Reading %d atoms in box (%.2f, %.2f, %.2f)\n", natoms, box[0],box[1],box[2]);

    em3->memory->allocate(xa, natoms, 3); 
    em3->memory->allocate(xa0, natoms, 3); 
    memory->allocate(type,natoms);
    for (int i=0; i<natoms; i++){
        for (int j=0; j<3; j++){
            xa[i][j] = 0.0;
        }
    }

    for (int i=0; i<natoms; i++){
        for (int j=0; j<3; j++){
            xa0[i][j] = 0.0;
        }
    }

    for (int i=0; i<natoms; i++){

        getline(config, line);
        //std::cout << line << std::endl;
        stringstream ss3(line);
        ss3 >> type[i] >> xa[i][0] >> xa[i][1] >> xa[i][2];
        
    }
    config.close();


    // Store original positions
    /*
    for (int i=0; i<natoms; i++){
        for (int j=0; j<3; j++){
            x0[i][j] = x[i][j];
        }
        //printf("%f %f %f\n", x[i][0],x[i][1],x[i][2]);
    }
    */

    /* Read EQUIL */
    
    ifstream readfile1;

    printf(" Reading EQUIL.\n");

    readfile1.open("EQUIL");
    //readfile.open("ev_real.txt");
  
    if (!readfile1.is_open()) {
        cout<<"Unable to open the file!"<<endl;
        exit(1);
    }

    for (int i=0; i<natoms; i++){

        getline(readfile1, line);
        //std::cout << line << std::endl;
        stringstream ss3(line);
        ss3 >> xa0[i][0] >> xa0[i][1] >> xa0[i][2];
        
    }
    
    readfile1.close();
    

    
    //ifstream readfile1;

    //printf(" Reading Amplitudes.\n");

    memory->allocate(xm,natoms*3);

    for (int i = 0; i < natoms*3; i++) {
        xm[i]=0.0;
    }

    //readfile1.open("AMPS");
    //readfile.open("ev_real.txt");
    /*
    if (!readfile1.is_open()) {
        cout<<"Unable to open the file!"<<endl;
        exit(1);
    }

    for (int i=0;i<3*natoms;i++){
        readfile1>>q[i];
        //printf(" %e\n", q[i]);
    }
    */
    //readfile1.close();

    // Convert everything to SI units
    for (int i=0; i<natoms; i++){
        for (int a=0; a<3; a++){
            xa[i][a] = xa[i][a]*1e-10;
            xa0[i][a] = xa0[i][a]*1e-10;

        }
    }
    

}

void Input::initialize()
{

    double k = 1.38064852e-23; // Boltzmann constant [J/K]

    // Convert inputs to SI units

    dt = dt*1e-15;

    //tau = sigma*sqrt(m/epsilon);
    //velocity = sqrt(epsilon/m);
    //force = epsilon/sigma;
    //pressure = epsilon/(sigma*sigma*sigma);
    //temperature = epsilon/k;
    /*
    printf("Units measured in:\n");
    printf("  length: %e m\n", sigma);
    printf("  energy: %e J\n", epsilon);
    printf("  mass: %e kg\n", m);
    printf("  time: %e s\n", tau );
    printf("  velocity: %e m/s\n", velocity );
    printf("  force: %e N\n", force );
    printf("  pressure: %e N/m^2\n", pressure );
    printf("  temperature: %f K\n", temperature );
    */

    volume = box[0]*box[1]*box[2];

    //temp = temp/temperature;

    //dt = dt/tau;
    
    //printf("Tau: %e\n", tau);
    //printf("Cutoff [sigma]: %f\n", rc);
    printf(" Temperature [K]: %f\n", temp);
    printf(" Box [A]: %f %f %f\n", box[0],box[1],box[2]);
    printf(" Timestep [s]: %e\n", dt);

    // Initialize MB velocities

    em3->memory->allocate(va, natoms, 3); 
    for (int i=0; i<natoms; i++){
        for (int j=0; j<3; j++){
            va[i][j] = 0.0;
        }
    }

    // random device class instance, source of 'true' randomness for initializing random seed
    std::random_device rd; 

    // Mersenne twister PRNG, initialized with seed from previous random device instance
    //std::mt19937 gen(rd()); 
    std::mt19937 gen(1); 
    
    // instance of class std::normal_distribution with specific mean and stddev
    std::normal_distribution<float> d(0.0, 1.0); 

    double normal;
    for (int i=0; i<natoms; i++){
        for (int j=0; j<3; j++){

            // get random number with normal distribution using gen as random source
            normal = d(gen); 
            //v[i][j] = sqrt(temp)*normal;
            va[i][j] = normal;
        }
    }

    // Scale each velocity component
    //printf(" %e %e %e\n", mass[0],mass[1],temp);
    for (int i=0; i<natoms; i++){
        for (int j=0; j<3; j++){
            va[i][j] = sqrt(k*temp/mass[type[i]-1])*va[i][j];
            //fprintf(fh_debug, "%e\n", va[i][j]);
        }
    }

    // Calculate velocity center of mass

    double vx_com = 0.0;
    double vy_com = 0.0;
    double vz_com = 0.0;
    double totalMass = 0.0;
    for (int i=0; i<natoms; i++){
        vx_com += mass[type[i]-1]*va[i][0];
        vy_com += mass[type[i]-1]*va[i][1];
        vz_com += mass[type[i]-1]*va[i][2];
        totalMass += mass[type[i]-1];
    }
    vx_com = vx_com/(totalMass);
    vy_com = vy_com/(totalMass);
    vz_com = vz_com/(totalMass);
    //fprintf(fh_debug, "%f %f %f\n", vx_com, vy_com, vz_com);

    // Subtract velocity COM from all components
    for (int i=0; i<natoms; i++){

        va[i][0] -= vx_com;
        va[i][1] -= vy_com;
        va[i][2] -= vz_com;

    }

    
    // Calculate velocity center of mass for sanity check
    vx_com = 0.0;
    vy_com = 0.0;
    vz_com = 0.0;
    for (int i=0; i<natoms; i++){
        vx_com += mass[type[i]-1]*va[i][0];
        vy_com += mass[type[i]-1]*va[i][1];
        vz_com += mass[type[i]-1]*va[i][2];
    }
    vx_com = vx_com/(totalMass);
    vy_com = vy_com/(totalMass);
    vz_com = vz_com/(totalMass);
    fprintf(fh_debug, "%f %f %f\n", vx_com, vy_com, vz_com);
    

    // Scale each velocity component
    //printf(" %e %e %e\n", mass[0],mass[1],temp);
    /*
    for (int i=0; i<natoms; i++){
        for (int j=0; j<3; j++){
            va[i][j] = sqrt(k*temp/mass[type[i]-1])*va[i][j];
            //fprintf(fh_debug, "%e\n", va[i][j]);
        }
    }
    */

    // Calculate velocity center of mass for sanity check
    /*
    for (int i=0; i<natoms; i++){
        vx_com += va[i][0];
        vy_com += va[i][1];
        vz_com += va[i][2];
    }
    vx_com = vx_com/(natoms);
    vy_com = vy_com/(natoms);
    vz_com = vz_com/(natoms);
    fprintf(fh_debug, "%f %f %f\n", vx_com, vy_com, vz_com);
    */
    

    
    // Kinetic energy is half times the sum of squared speeds times mass

    double ke=0.0;
    for (int i=0; i<natoms; i++){

        ke += mass[type[i]-1]*(va[i][0]*va[i][0] + va[i][1]*va[i][1] + va[i][2]*va[i][2])/2.0;

    }

    printf(" ke: %e eV.\n", ke*6.242e+18);

    double temp_actual = ke*(2.0/3.0)/(natoms*k);
    //printf(" temp/temp_actual: %f\n", temp/temp_actual);

    // Recale the velocities
    for (int i=0; i<natoms; i++){
        va[i][0] = va[i][0]*sqrt(temp/temp_actual);
        va[i][1] = va[i][1]*sqrt(temp/temp_actual);
        va[i][2] = va[i][2]*sqrt(temp/temp_actual);
    }

    ke=0.0;
    for (int i=0; i<natoms; i++){

        ke += mass[type[i]-1]*(va[i][0]*va[i][0] + va[i][1]*va[i][1] + va[i][2]*va[i][2])/2.0;

    }

    printf(" ke: %e eV.\n", ke*6.242e+18);

    // Convert to mode velocities

    memory->allocate(vm,natoms*3);

    if (space==1 || space==2){
        for (int n=0; n<natoms*3; n++){
            //printf(" n: %d\n", n);
            vm[n]=0.0;
            for (int j=0; j<natoms; j++){
                //printf("  j: %d\n", j);
                vm[n] += sqrt(mass[type[j]-1])*emat[0+3*j][n]*va[j][0];
                vm[n] += sqrt(mass[type[j]-1])*emat[1+3*j][n]*va[j][1];
                vm[n] += sqrt(mass[type[j]-1])*emat[2+3*j][n]*va[j][2];
            }
        }
    }

    /*
    for (int n=0; n<natoms*3; n++){
        printf(" v[%d]: %e\n", n,v[n]);
    }
    */

    /*
    for (int n=0; n<natoms; n++){
        printf(" va: %f %f %f\n", va[n][0],va[n][1],va[n][2]);
    }
    */

    // Initially excite one mode as a test
    /*
    for (int n=0; n<3*natoms; n++){
        if (n==3){
            v[n]=1e-10;
        }
        else v[n]=0.0;
    }
    */

  // Excite a single mode

  

}

/*
Read frequencies and MCCs.
*/
void Input::readParams()
{

    if (space == 1 || space == 2){

        ifstream readfile1;
        ifstream readfile2;
        
        printf(" Reading Eigenvectors.\n");

        memory->allocate(emat,natoms*3,natoms*3);

        for (int i = 0; i < natoms*3; i++) {
            for (int j = 0; j < natoms*3; j++) {
                emat[i][j] = 0.0;
            }
        }

        readfile2.open("EMAT");
        //readfile2.open("ev_real.txt");
      
        if (!readfile2.is_open()) {
            cout<<"Unable to open the file!"<<endl;
            exit(1);
        }

        //printf("natoms: %d\n",  natoms);
        for (int i=0;i<3*natoms;i++){
            for (int j=0;j<3*natoms;j++){
                readfile2>>emat[i][j];
            }
        }

        /*
        // Add random numbers to eigenvectors
        std::srand(std::time(nullptr)); // use current time as seed for random generator
        for (int i=0;i<3*natoms;i++){
            for (int j=0;j<3*natoms;j++){
                emat[i][j] += ((double) rand() / (RAND_MAX)) + 2.0;
                //printf("%f\n", ((double) rand() / (RAND_MAX)));
            }
        }
        */
        
        readfile2.close();

        int i,j,k,l;
        int a,b,c,d;

        if (order >=2){
            /* Read MCC2s */
            ifstream fh("MCC2");   
          
            string line;
            getline(fh,line);
            stringstream ss(line);
            ss >> nmcc2;
            printf(" Reading %d MCC2s.\n", nmcc2);

            memory->allocate(mcc2,nmcc2);
            double value;
            int counter = 0;
            
            while (getline(fh, line))
            {
                stringstream ss(line);
                ss >> i >> j >> value;
                mcc2[counter].i=i;
                mcc2[counter].val=value;
                counter++;

            }
            fh.close();
        }
        
        if (order >=3){
            /* Read MCC3s */
            ifstream fh("MCC3");   
            string line;
            getline(fh,line);
            std::cout<<line<<std::endl;
            stringstream ss3(line);
            ss3 >> nmcc3;
            printf(" Reading %d MCC3s.\n", nmcc3);

            memory->allocate(mcc3,nmcc3);
            int counter = 0;
            double value;
            while (getline(fh, line))
            {
                stringstream ss(line);
                ss >> i >> j >> k >> value;
                mcc3[counter].i=i;
                mcc3[counter].j=j;
                mcc3[counter].k=k;
                mcc3[counter].val=value;
                counter++;

            }
            fh.close();
        }

        if (order >= 4){
            /* Read MCC4s */
            ifstream fh("MCC4");   
          
            string line;
            getline(fh,line);
            stringstream ss4(line);
            ss4 >> nmcc4;
            printf(" Reading %d MCC4s.\n", nmcc4);

            memory->allocate(mcc4,nmcc4);
            int counter = 0;
            double value;
            while (getline(fh, line))
            {
                stringstream ss(line);
                ss >> i >> j >> k >> l >> value;
                mcc4[counter].i=i;
                mcc4[counter].j=j;
                mcc4[counter].k=k;
                mcc4[counter].l=l;
                mcc4[counter].val=value;
                counter++;

            }
            fh.close();
        }

    }

    if (space==0 || space==2){

        int i,j,k,l;
        int a,b,c,d; 

        if (order >=2){
            /* Read IFC2s */
            ifstream fh("FC2_FD");
            string line;
            nfc2 = -1; // -1 to subtract first line
            while (getline(fh, line))
            {
                nfc2++;

            }
            fh.close();
            printf(" Found %d FC2s.\n", nfc2);
            memory->allocate(fc2,nfc2);
            fh.open("FC2_FD");
            getline(fh,line);
            int counter=0;
            double value;
            while (getline(fh, line))
            {
                stringstream ss(line);
                ss >> i >> a >> j >> b >> value;
                fc2[counter].i=i-1;
                fc2[counter].j=j-1;
                fc2[counter].a=a-1;
                fc2[counter].b=b-1;
                fc2[counter].val=value;
                counter++;

            }
            fh.close();
        }

        if (order >=3){
            /* Read IFC3s */
            ifstream fh("FC3_FD");
            string line;
            nfc3 = -1; // -1 to subtract first line
            while (getline(fh, line))
            {
                nfc3++;

            }
            fh.close();
            printf(" Found %d FC3s.\n", nfc3);
            memory->allocate(fc3,nfc3);
            fh.open("FC3_FD");
            getline(fh,line);
            int counter=0;
            double value;
            while (getline(fh, line))
            {
                stringstream ss(line);
                ss >> i >> a >> j >> b >> k >> c >> value;
                fc3[counter].i=i-1;
                fc3[counter].j=j-1;
                fc3[counter].k=k-1;
                fc3[counter].a=a-1;
                fc3[counter].b=b-1;
                fc3[counter].c=c-1;
                fc3[counter].val=value;
                counter++;

            }
            fh.close();
        }

        if (order >=4){
            /* Read IFC4s */
            ifstream fh("FC4_FD");
            string line;
            nfc4 = -1; // -1 to subtract first line
            while (getline(fh, line))
            {
                nfc4++;

            }
            fh.close();
            printf(" Found %d FC4s.\n", nfc4);
            memory->allocate(fc4,nfc4);
            fh.open("FC4_FD");
            getline(fh,line);
            int counter=0;
            double value;
            while (getline(fh, line))
            {
                stringstream ss(line);
                ss >> i >> a >> j >> b >> k >> c >> l >> d >> value;
                fc4[counter].i=i-1;
                fc4[counter].j=j-1;
                fc4[counter].k=k-1;
                fc4[counter].l=l-1;
                fc4[counter].a=a-1;
                fc4[counter].b=b-1;
                fc4[counter].c=c-1;
                fc4[counter].d=d-1;
                fc4[counter].val=value;
                counter++;

            }
            fh.close();
        }

    }
    
    

}


