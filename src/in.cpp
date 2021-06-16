#include <string>
#include <sstream>
#include <vector>
#include <algorithm>    // std::find
#include <fstream>
#include "mpi.h"
#include <cmath>

#include "in.h"
#include "mem.h"
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

#include <cstddef> // for MPI offsets

using namespace std;

using namespace LAMMPS_NS;

using namespace MC_NS;

In::In(MC *mc) : Ptrs(mc) {

    rank = mc->rank;

    //char debug[64];
    //sprintf (debug, "debug_in/DEBUG_IN%d", rank);

    //fh_debug = fopen(debug,"w");
}

In::~In() 
{

    mem->deallocate(x0);

    //fclose(fh_debug);
};


void In::readInput()
{

  //std::cout << "Reading INPUT on proc " << pops->rank << std::endl;

  string line;

  // Declare scalar inputs
  double value;

  // Open INPUT file
  ifstream INPUT("INPUT");
  // Ignore the first line
  getline(INPUT, line); 
  string characters;

  // Now loop through the LAMMPS setup commands, and store the commands in "commands", delimited by a ","
  getline(INPUT, line);
  while (line != "----------------------------------------LAMMPS PAIR STYLE COMMANDS")
  {
    commands.append(line);
    commands.append(",");
    getline(INPUT,line);
    //debug << commands << endl;
  }
  // Continue looping through to get the pair style commands
  getline(INPUT, line);
  while (line != "----------------------------------------END OF FILE")
  {
    nontab_commands.append(line);
    nontab_commands.append(",");
    getline(INPUT,line);
    //debug << commands << endl;
  }

  INPUT.close();

  // Loop through LAMMPS commands and input them 
  size_t n = commands.find(',');
  string substring;
  string rest = commands;
  //writefile << commands << endl;
  while (n != std::string::npos)
  {
    substring = rest.substr(0,n);
    // Convert strings to const char * types so that lammps can read them
    const char* substring_cc = substring.c_str();
    //debug << "SUBSTRING: " << substring_cc << endl;
    // Input the command
    lmp->input->one(substring_cc);
    //debug << "input succeeded" << endl;
    // Find location of next command and repeat
    rest = rest.substr(n+1);
    //debug << "REST: " << rest << endl;
    n = rest.find(',');
    //debug << "n: " << n << endl;
  }
  //writefile << "LAMMPS SETUP COMMANDS SUCCESSFULLY INPUTTED" << endl;

  // Loop through LAMMPS commands and input them 
  n = nontab_commands.find(',');
  rest = nontab_commands;
  //writefile << commands << endl;
  while (n != std::string::npos)
  {
    substring = rest.substr(0,n);
    // Convert strings to const char * types so that lammps can read them
    const char* substring_cc = substring.c_str();
    //debug << "SUBSTRING: " << substring_cc << endl;
    // Input the command
    lmp->input->one(substring_cc);
    //debug << "input succeeded" << endl;
    // Find location of next command and repeat
    rest = rest.substr(n+1);
    //debug << "REST: " << rest << endl;
    n = rest.find(',');
    //debug << "n: " << n << endl;
  }
  
    /* Store original positions as equilibrium positions*/
    
    int natoms = lmp->atom->natoms;
    mem->allocate(x0,natoms,3);
    for (int i = 0; i < natoms; i++) {
        for (int j = 0; j < 3; j++) {
            x0[i][j] = lmp->atom->x[i][j];
        }
    }
    

}

void In::readConfig()
{
    /* Read CONFIG file */

    /*
    ifstream config("CONFIG");
    string line;

    getline(config, line);
    stringstream ss(line);
    ss >> natoms >> ntypes;
    printf(" %d atom types.\n", ntypes);
    mem->allocate(mass,ntypes);

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

    for (int i=0; i<natoms; i++){

        getline(config, line);
        //std::cout << line << std::endl;
        stringstream ss3(line);
        ss3 >> type[i];// >> xa[i][0] >> xa[i][1] >> xa[i][2];
        
    }
    config.close();

    */
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
    /*
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
    */
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
    /*
    // Convert everything to SI units
    for (int i=0; i<natoms; i++){
        for (int a=0; a<3; a++){
            xa[i][a] = xa[i][a]*1e-10;
            xa0[i][a] = xa0[i][a]*1e-10;

        }
    }
    */

}


void In::evaluate()
{

    int mode = 200;

    int natoms = lmp->atom->natoms;
    double *pe_p;
    double pe;

    double *x; 

    mem->allocate(x,natoms*3);

    

    FILE * fh_xyz;
    fh_xyz = fopen("dump.xyz","w");

    /*
    int eindx=0;
    for (int n=0; n<natoms; n++){
        fprintf(fh_debug, "%f %f %f\n", emat[0+eindx][mode], emat[1+eindx][mode], emat[2+eindx][mode]);
        eindx += 3;
    }
    */

    for (int d=-1.0*ndisps+1; d<ndisps; d++){

        fprintf(fh_xyz, "%d\n", natoms);
        fprintf(fh_xyz, "Atoms. Timestep: %d\n", d-1);

        int eindx = 0;
        for (int n=0; n<natoms; n++){
            x[eindx+0] = x0[n][0]+d*delta*emat[0+eindx][mode];
            x[eindx+1] = x0[n][1]+d*delta*emat[1+eindx][mode];
            x[eindx+2] = x0[n][2]+d*delta*emat[2+eindx][mode];
            fprintf(fh_xyz, "%d %f %f %f\n", lmp->atom->type[n], x[eindx+0],x[eindx+1],x[eindx+2]);

            eindx+=3;
        }

        // Scatter the atoms with the new positions

        lammps_scatter_atoms(lmp,"x",1,3,x);

        lmp->input->one("run 0");

        // Extract the energy

        pe_p = (double *) lammps_extract_compute(lmp, "P", 0, 0); // style = 0 for global data, type = 0 for scalar quantity
        pe = pe_p[0];
        printf("pe: %f\n", pe);

    }


    mem->deallocate(x);

    fclose(fh_xyz);


}


void In::calcFC2()
{


    char fc2name[64];
    char fc3name[64];
    char fc4name[64];
    sprintf (fc2name, "FC2_PROC%d", rank);
    sprintf (fc3name, "FC3_PROC%d", rank);
    sprintf (fc4name, "FC4_PROC%d", rank);
    //std::cout << debug << std::endl;
    //fh_debug = fopen(debug, "w");

    //printf(" Rank: %d\n", rank);

    delta = mc->delta;
    rc = mc->cutoff;
    order = mc->order;
    // Tolerances are hardcoded for now... May change later.
    tol2 = 1e-10;
    tol3 = 1e-10;
    tol4 = 1e-10;
    //printf("delta: %f\n",delta);
    int natoms = lmp->atom->natoms;
    double *pe_p;
    double pe;
    double pe1,pe2,pe3,pe4;
    //double *x; 
    double **hessian;
    double phi;

    double ux,uy,uz;
    double sum_sqd;
    double msd;
    double sum_msd = 0;

    //mem->allocate(x,natoms*3);

    //FILE * fh_msd;
    //fh_msd = fopen("mode_sqrtmsd.txt", "w");

    //mem->allocate(hessian, natoms*3, natoms*3);
    int i,j,k,l;
    int a,b,c,d;
    
    double rij,rik,ril;

    
	if (order ==2){

        if (rank==0) fh_fc2 = fopen("FC2","w");

        if (rank==0) printf(" Calculating FC2.\n");
        if (rank==0) printf(" Tolerance = %e.\n",tol2);

        int nfc2=0;
		i=0;
		j=0;
		for (int ni=0; ni<natoms; ni++){
			//if (rank==0) printf(" ni: %d\n", ni);
			for (int a=0; a<3; a++){
                i = 3*ni+a;
				for (int nj=0; nj<natoms;nj++){
					rij = calcRij(x0[ni][0],x0[ni][1],x0[ni][2],x0[nj][0],x0[nj][1],x0[nj][2]);
					if (rij<rc){
					    for (int b=0; b<3;b++){
						    j=3*nj+b;
						    if (j!=i){ // Don't include i=j here, since we'll do ASR corrections after
								nfc2++;
							}
						}
					}
				}
			}
		}
        if (rank==0) printf(" %d FC2s.\n", nfc2);

        int *nepp;
        mem->allocate(nepp, mc->nprocs); // number elements (FCs) per proc
        for (int p=0; p<mc->nprocs; p++){
            nepp[p] = nfc2/mc->nprocs;
        }

        // divide up the remainder
        for (int p=0; p<(nfc2 % mc->nprocs); p++){
            nepp[p] += 1;
        }

        if (rank==0){
            printf(" Splitting FC2s on procs like:\n");
            for (int p=0; p<mc->nprocs; p++){
                printf("  %d FC2s on proc %d.\n", nepp[p],p);
            }
        }

        int start_indx = 0;
        for (int p=0; p<rank; p++){
            start_indx += nepp[p];
        }
        int end_indx = 0; //napp[0]-1;
        for (int p=0; p<rank+1; p++){
            end_indx += nepp[p];
        }
        end_indx=end_indx-1;
        printf("  Proc %d: %d-%d.\n", rank,start_indx,end_indx);


        if (rank==0) printf(" Creating FC2 list on each proc.\n");
        mem->allocate(fc2,nepp[rank]);

        int fc2_count=0;
        int fc2_indx=0; // indx of a FC2 on this proc
		i=0; // i and j represent an atom-direction pair in this loop
		j=0;
		for (int ni=0; ni<natoms; ni++){
			//if (rank==0) printf(" ni: %d\n", ni);
			for (int a=0; a<3; a++){
                i = 3*ni+a;
				for (int nj=0; nj<natoms;nj++){
					rij = calcRij(x0[ni][0],x0[ni][1],x0[ni][2],x0[nj][0],x0[nj][1],x0[nj][2]);
					if (rij<rc){
					    for (int b=0; b<3;b++){
						    j=3*nj+b;
						    if (j!=i){ // Don't include i=j here, since we'll do ASR corrections after
								//fprintf(fh_debug, "%d %d\n", i,j);
                                //if (rank==0) printf(" %d\n", fc2_indx);
                                if (start_indx<=fc2_count && fc2_count <= end_indx){
                                    fc2[fc2_indx].i=ni;
                                    fc2[fc2_indx].j=nj;
                                    fc2[fc2_indx].a=a;
                                    fc2[fc2_indx].b=b;
                                    //fprintf(fh_debug, " %d %d %d %d %d %d\n", fc2_indx, fc2_count,ni,a,nj,b);
                                    fc2_indx++;
                                }
                                fc2_count++;
							}
						}
					}
				}
			}
		}

        // Loop over all FC2s on this proc, calculate and store FC
        // Also count the number of FCs above tolerance value on each proc
        int nfc2_tol_proc = 0;
        printf(" Progress from 0-%d:\n", nepp[rank]);
        for (int w=0; w<nepp[rank]; w++){
            if (rank==0 && w%10==0) printf(" %d\n", w);
            i=fc2[w].i;
            j=fc2[w].j;
            a=fc2[w].a;
            b=fc2[w].b;
            fc2[w].fc = calcD1(i,j,a,b);
            if (abs(fc2[w].fc)>tol2){
                nfc2_tol_proc++;
                //printf("%d %d %d %d\n", i+1,a+1,j+1,b+1);
            }
            
            
        }

        // Gatherv all fc2 arrays onto the root proc
        MPI_Aint disp[5] = {offsetof(fc2_struct, i), offsetof(fc2_struct, j), offsetof(fc2_struct,a), offsetof(fc2_struct,b), offsetof(fc2_struct,fc)};
        int blocklen[5] = {1,1,1,1,1};
        MPI_Datatype type[5] = {MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_DOUBLE};
        MPI_Datatype fc2_dt; // fc2 datatype
        MPI_Type_create_struct(5, blocklen, disp, type, &fc2_dt);
        MPI_Type_commit(&fc2_dt);
        
        // Gathering the data
        fc2_struct *fc2_all = nullptr;
        
        int nprocs = mc->nprocs;
        
        if (rank == 0){
            fc2_all = new fc2_struct[nfc2];
        }
        //MPI_Gather(fc2, nepp[rank], fc2_dt, fc2_all, nepp[rank], fc2_dt, 0, MPI_COMM_WORLD); // only works if every proc has same number of fc2

        // Calculate receive displacements
        int *receive_displacements;
        mem->allocate(receive_displacements,nprocs);
        //int receive_displacements[4];//  = { 0, 0, 1, 3 };
        receive_displacements[0]=0;
        int location = 0;
        for (int w=1; w<nprocs; w++){
            location += nepp[w-1];
            //printf(" %d\n", location);
            receive_displacements[w]=location;
        }

        MPI_Gatherv(fc2,nepp[rank],fc2_dt,fc2_all,nepp,receive_displacements,fc2_dt,0,MPI_COMM_WORLD);

        /* Count the number of FCs above tolerance value on each proc */
        /*
        int nfc2_tol_proc = 0;
        for (int w=0; w<nepp[rank]; w++){
            //printf(" %f\n", fc2[w].fc);
            if (abs(fc2[w].fc) > tol2){
                //printf(" YUP!\n");
                nfc2_tol_proc++;
            }
        }
        */
        //printf(" Proc %d, %d\n", rank,nfc2_tol_proc);

        // Sum all nfc2_tol_proc to the root proc
        int nfc2_tol_total;
        MPI_Reduce(&nfc2_tol_proc, &nfc2_tol_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank==0) printf(" Total FC2 above tol: %d\n", nfc2_tol_total);

        // Print the FCs that are above the tol
        if (rank==0){
            fprintf(fh_fc2, "%d\n",nfc2_tol_total);
            for (int w=0; w<nfc2; w++){
                if (abs(fc2_all[w].fc)>tol2){
                    fprintf(fh_fc2, "%d    %d    %d    %d    %.10e\n", fc2_all[w].i+1,fc2_all[w].a+1,fc2_all[w].j+1,fc2_all[w].b+1,fc2_all[w].fc);
                }
            }
        }
        
        mem->deallocate(fc2);
        mem->deallocate(nepp);
        mem->deallocate(receive_displacements);
        if (rank==0) fclose(fh_fc2);

	} // if order==2


	/* Need a separate 3rd order loop, since i,j,k must all be within cutoff */
	if (order ==3){

        if (rank==0) fh_fc3 = fopen("FC3","w");
        //fprintf(fh_fc3, "atom1 xyz1 atom2 xyz2 atom3 xyz3 ifc3\n");

        if (rank==0) printf(" Calculating FC3.\n");
        if (rank==0) printf(" Tolerance = %e.\n",tol3);

        int nfc3=0;
		int ii,jj,kk;
        double rjk;
		for (int ni=0; ni<natoms; ni++){
			//if (rank==0) printf(" ni: %d\n", ni);
			for (int a=0; a<3; a++){
                ii = 3*ni+a;
				for (int nj=0; nj<natoms; nj++){
					rij = calcRij(x0[ni][0],x0[ni][1],x0[ni][2],x0[nj][0],x0[nj][1],x0[nj][2]);
					if (rij<rc){
					    for (int b=0; b<3; b++){
						    jj = 3*nj+b;
						    for (int nk=0; nk<natoms; nk++){

						        rik = calcRij(x0[ni][0],x0[ni][1],x0[ni][2],x0[nk][0],x0[nk][1],x0[nk][2]);
                                rjk = calcRij(x0[nj][0],x0[nj][1],x0[nj][2],x0[nk][0],x0[nk][1],x0[nk][2]);
						        if (rik<rc && rjk<rc){
							        for (int c=0;c<3;c++){
								        kk=3*nk+c;
                                        //fprintf(fh_debug, "%d %d %d\n", i,j,k);
								        //if (jj>=ii && kk>=jj){

                                            /*
                                            if (ii==jj && jj==kk){
                                                fprintf(fh_debug, "%d %d %d ----- iii\n", ii,jj,kk);
                                            }
                                            else if (ii==jj && jj!=kk){
                                                fprintf(fh_debug, "%d %d %d ----- iik\n", ii,jj,kk);
                                                fprintf(fh_debug, "%d %d %d\n", ii,kk,ii);
                                                fprintf(fh_debug, "%d %d %d\n", kk,ii,ii);
                                            }
                                            else if (ii!=jj && jj==kk){
                                                fprintf(fh_debug, "%d %d %d ----- ijj\n", ii,jj,kk);
                                                fprintf(fh_debug, "%d %d %d\n", jj,ii,jj);
                                                fprintf(fh_debug, "%d %d %d\n", jj,jj,ii);
                                            }
                                            else if (ii!=jj && jj!=kk){
                                                fprintf(fh_debug, "%d %d %d ----- ijk\n", ii,jj,kk);
                                                fprintf(fh_debug, "%d %d %d\n", ii,kk,jj);
                                                fprintf(fh_debug, "%d %d %d\n", jj,ii,kk);
                                                fprintf(fh_debug, "%d %d %d\n", jj,kk,ii);
                                                fprintf(fh_debug, "%d %d %d\n", kk,ii,jj);
                                                fprintf(fh_debug, "%d %d %d\n", kk,jj,ii);
                                            }
                                            */

											nfc3++;

										//}
									}
									
								}
							}
						}
					}
				}
			}
		}
        if (rank==0) printf(" %d FC3s.\n",nfc3);
        
        int *nepp;
        mem->allocate(nepp, mc->nprocs); // number elements (FCs) per proc
        for (int p=0; p<mc->nprocs; p++){
            nepp[p] = nfc3/mc->nprocs;
        }

        // divide up the remainder
        for (int p=0; p<(nfc3 % mc->nprocs); p++){
            nepp[p] += 1;
        }

        if (rank==0){
            printf(" Splitting FC3s on procs like:\n");
            for (int p=0; p<mc->nprocs; p++){
                printf("  %d FC3s on proc %d.\n", nepp[p],p);
            }
        }

        int start_indx = 0;
        for (int p=0; p<rank; p++){
            start_indx += nepp[p];
        }
        int end_indx = 0; //napp[0]-1;
        for (int p=0; p<rank+1; p++){
            end_indx += nepp[p];
        }
        end_indx=end_indx-1;
        //printf("  Proc %d: %d-%d.\n", rank,start_indx,end_indx);


        if (rank==0) printf(" Creating FC3 list on each proc.\n");
        mem->allocate(fc3,nepp[rank]);

        int fc3_count=0;
        int fc3_indx=0; // indx of a FC3 on this proc
		i=0;
		j=0;
		k=0;
		for (int ni=0; ni<natoms; ni++){
			//if (rank==0) printf(" ni: %d\n", ni);
			for (int a=0; a<3; a++){
                ii = 3*ni+a;
				for (int nj=0; nj<natoms; nj++){
					rij = calcRij(x0[ni][0],x0[ni][1],x0[ni][2],x0[nj][0],x0[nj][1],x0[nj][2]);
					if (rij<rc){
					    for (int b=0; b<3; b++){
						    jj = 3*nj+b;
						    for (int nk=0; nk<natoms; nk++){

						        rik = calcRij(x0[ni][0],x0[ni][1],x0[ni][2],x0[nk][0],x0[nk][1],x0[nk][2]);
                                rjk = calcRij(x0[nj][0],x0[nj][1],x0[nj][2],x0[nk][0],x0[nk][1],x0[nk][2]);
						        if (rik<rc && rjk<rc){
							        for (int c=0;c<3;c++){
								        kk=3*nk+c;
                                        //fprintf(fh_debug, "%d %d %d\n", i,j,k);
								        //if (jj>=ii && kk>=jj){
									

											//fprintf(fh_debug, "%d %d %d %f %f\n", i,j,k, rij, rik);
                                            if (start_indx<=fc3_count && fc3_count <= end_indx){
                                                //if (rank==0) printf(" %d\n", fc3_indx);
                                                fc3[fc3_indx].i=ni;
                                                fc3[fc3_indx].j=nj;
                                                fc3[fc3_indx].k=nk;
                                                fc3[fc3_indx].a=a;
                                                fc3[fc3_indx].b=b;
                                                fc3[fc3_indx].c=c;
                                                //fprintf(fh_debug, " %d %d %d %d %d %d\n", fc2_indx, fc2_count,ni,a,nj,b);
                                                fc3_indx++;
                                            }
                                            fc3_count++;

										//}
									}
									
								}
							}
						}
					}
				}
			}
		}

        // Loop over all FC3s on this proc, calculate and store FC
        int nfc3_tol_proc = 0;
        for (int w=0; w<nepp[rank]; w++){
            if (rank==0 && w%10==0) printf(" %d\n", w);
            i=fc3[w].i;
            j=fc3[w].j;
            k=fc3[w].k;
            a=fc3[w].a;
            b=fc3[w].b;
            c=fc3[w].c;
            fc3[w].fc = calcD2(i,j,k,a,b,c);
            if (abs(fc3[w].fc) > tol3) nfc3_tol_proc++;
        }

        // Gatherv all fc3 arrays onto the root proc
        MPI_Aint disp[7] = {offsetof(fc3_struct, i), offsetof(fc3_struct, j), offsetof(fc3_struct, k), \
                            offsetof(fc3_struct, a), offsetof(fc3_struct, b), offsetof(fc3_struct, c), \
                            offsetof(fc3_struct, fc) };
        int blocklen[7] = {1,1,1,1,1,1,1};
        MPI_Datatype type[7] = {MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_DOUBLE};
        MPI_Datatype fc3_dt; // fc3 datatype
        MPI_Type_create_struct(7, blocklen, disp, type, &fc3_dt);
        MPI_Type_commit(&fc3_dt);
        
        // Gathering the data
        fc3_struct *fc3_all = nullptr;
        
        int nprocs = mc->nprocs;
        
        if (rank == 0){
            fc3_all = new fc3_struct[nfc3];
        }

        // Calculate receive displacements
        int *receive_displacements;
        mem->allocate(receive_displacements,nprocs);
        //int receive_displacements[4];//  = { 0, 0, 1, 3 };
        receive_displacements[0]=0;
        int location = 0;
        for (int w=1; w<nprocs; w++){
            location += nepp[w-1];
            //printf(" %d\n", location);
            receive_displacements[w]=location;
        }

        MPI_Gatherv(fc3,nepp[rank],fc3_dt,fc3_all,nepp,receive_displacements,fc3_dt,0,MPI_COMM_WORLD);

        //printf(" Proc %d, %d\n", rank,nfc2_tol_proc);

        // Sum all nfc3_tol_proc to the root proc
        int nfc3_tol_total;
        MPI_Reduce(&nfc3_tol_proc, &nfc3_tol_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank==0) printf(" Total FC3 above tol: %d\n", nfc3_tol_total);

        // Print the FCs that are above the tol
        //int ii,jj,kk;
        if (rank==0){
            fprintf(fh_fc3, "%d\n",nfc3_tol_total);
            for (int w=0; w<nfc3; w++){
                if (abs(fc3_all[w].fc)>tol3){
                    ii=3*fc3_all[w].i+fc3_all[w].a;
                    jj=3*fc3_all[w].j+fc3_all[w].b;
                    kk=3*fc3_all[w].k+fc3_all[w].c;
                    fprintf(fh_fc3, "%d %d %d %d %d %d %.10e\n", fc3_all[w].i+1,fc3_all[w].a+1,fc3_all[w].j+1,fc3_all[w].b+1, \
                                                           fc3_all[w].k+1,fc3_all[w].c+1,fc3_all[w].fc);
                    //fprintf(fh_fc3, "%d %d %d %.10e\n", ii,jj,kk,fc3_all[w].fc);
                }
            }
        }
        
        mem->deallocate(fc3);
        mem->deallocate(nepp);
        mem->deallocate(receive_displacements);
        if (rank==0) fclose(fh_fc3);
        
	}

	if (order ==4){

        if (rank==0) fh_fc4 = fopen("FC4","w");
        //fprintf(fh_fc4, "atom1 xyz1 atom2 xyz2 atom3 xyz3 atom4 xyz4 ifc4\n");

        if (rank==0) printf(" Calculating FC4.\n");
        if (rank==0) printf(" Tolerance = %e.\n",tol4);

        int nfc4=0;
		i=0;
		j=0;
		k=0;
        l=0;
        double rjk,rjl,rkl;
		for (int ni=0; ni<natoms; ni++){
			//if (rank==0) printf(" ni: %d\n", ni);
			for (int a=0; a<3; a++){
                i = 3*ni+a;
				for (int nj=0; nj<natoms; nj++){

			        rij = calcRij(x0[ni][0],x0[ni][1],x0[ni][2],x0[nj][0],x0[nj][1],x0[nj][2]);
			        if (rij<rc){

					    for (int b=0; b<3; b++){
						    j = 3*nj+b;
						    for (int nk=0; nk<natoms; nk++){

					            rik = calcRij(x0[ni][0],x0[ni][1],x0[ni][2],x0[nk][0],x0[nk][1],x0[nk][2]);
                                rjk = calcRij(x0[nj][0],x0[nj][1],x0[nj][2],x0[nk][0],x0[nk][1],x0[nk][2]);
					            if (rik<rc && rjk<rc){

							        for (int c=0;c<3;c++){
								        k=3*nk+c;
                                        //fprintf(fh_debug, "%d %d %d\n", i,j,k);
                                        for (int nl=0; nl<natoms; nl++){

                                            ril = calcRij(x0[ni][0],x0[ni][1],x0[ni][2],x0[nl][0],x0[nl][1],x0[nl][2]);
                                            rjl = calcRij(x0[nj][0],x0[nj][1],x0[nj][2],x0[nl][0],x0[nl][1],x0[nl][2]);
                                            rkl = calcRij(x0[nk][0],x0[nk][1],x0[nk][2],x0[nl][0],x0[nl][1],x0[nl][2]);
							                if (ril<rc && rjl<rc && rkl<rc){

                                                for (int d=0;d<3;d++){
                                                    l=3*nl+d;
								                    //if (j>=i && k>=j && l>=k){
									
                                                                    nfc4++;

                                                    //}
								                }
									        }
									
								        }
                                    }
                                }
				            }
					    }
					}
				}
			}
		}
        /*
		for (int ni=0; ni<natoms; ni++){
			if (rank==0) printf(" ni: %d\n", ni);
			for (int a=0; a<3; a++){
                i = 3*ni+a;
				for (int nj=0; nj<natoms; nj++){

                    for (int b=0; b<3; b++){
                        j = 3*nj+b;
                        for (int nk=0; nk<natoms; nk++){

                            for (int c=0; c<3;c++){

                                k=3*nk+c;

                                for (int nl=0; nl<natoms; nl++){

                                    for (int d=0; d<3; d++){
            
                                        l=3*nl+d;

                                        if (j>=i && k>=j && l>=k){

                                            rij = calcRij(x0[ni][0],x0[ni][1],x0[ni][2],x0[nj][0],x0[nj][1],x0[nj][2]);
			                                if (rij<rc3){
                                                rik = calcRij(x0[ni][0],x0[ni][1],x0[ni][2],x0[nk][0],x0[nk][1],x0[nk][2]);
                                                if (rik<rc3){
                                                    ril = calcRij(x0[ni][0],x0[ni][1],x0[ni][2],x0[nl][0],x0[nl][1],x0[nl][2]);
							                        if (ril<rc3){
									
                                                        nfc4++;

                                                    }
								                }
									        }
									
								        }
                                    }
                                }
				            }
					    }
					}
				}
			}
		}
        */
        if (rank==0) printf(" %d FC4s.\n",nfc4);

        int *nepp;
        mem->allocate(nepp, mc->nprocs); // number elements (FCs) per proc
        for (int p=0; p<mc->nprocs; p++){
            nepp[p] = nfc4/mc->nprocs;
        }

        // divide up the remainder
        for (int p=0; p<(nfc4 % mc->nprocs); p++){
            nepp[p] += 1;
        }

        if (rank==0){
            printf(" Splitting FC4s on procs like:\n");
            for (int p=0; p<mc->nprocs; p++){
                printf("  %d FC4s on proc %d.\n", nepp[p],p);
            }
        }

        int start_indx = 0;
        for (int p=0; p<rank; p++){
            start_indx += nepp[p];
        }
        int end_indx = 0; //napp[0]-1;
        for (int p=0; p<rank+1; p++){
            end_indx += nepp[p];
        }
        end_indx=end_indx-1;
        //printf("  Proc %d: %d-%d.\n", rank,start_indx,end_indx);


        if (rank==0) printf(" Creating FC4 list on each proc.\n");
        mem->allocate(fc4,nepp[rank]);

        int fc4_count=0;
        int fc4_indx=0; // indx of a FC4 on this proc
		i=0;
		j=0;
		k=0;
        l=0;
        
		for (int ni=0; ni<natoms; ni++){
			//if (rank==0) printf(" ni: %d\n", ni);
			for (int a=0; a<3; a++){
                i = 3*ni+a;
				for (int nj=0; nj<natoms; nj++){

			        rij = calcRij(x0[ni][0],x0[ni][1],x0[ni][2],x0[nj][0],x0[nj][1],x0[nj][2]);
			        if (rij<rc){

					    for (int b=0; b<3; b++){
						    j = 3*nj+b;
						    for (int nk=0; nk<natoms; nk++){

					            rik = calcRij(x0[ni][0],x0[ni][1],x0[ni][2],x0[nk][0],x0[nk][1],x0[nk][2]);
                                rjk = calcRij(x0[nj][0],x0[nj][1],x0[nj][2],x0[nk][0],x0[nk][1],x0[nk][2]);
					            if (rik<rc && rjk<rc){

							        for (int c=0;c<3;c++){
								        k=3*nk+c;
                                        //fprintf(fh_debug, "%d %d %d\n", i,j,k);
                                        for (int nl=0; nl<natoms; nl++){

                                            ril = calcRij(x0[ni][0],x0[ni][1],x0[ni][2],x0[nl][0],x0[nl][1],x0[nl][2]);
                                            rjl = calcRij(x0[nj][0],x0[nj][1],x0[nj][2],x0[nl][0],x0[nl][1],x0[nl][2]);
                                            rkl = calcRij(x0[nk][0],x0[nk][1],x0[nk][2],x0[nl][0],x0[nl][1],x0[nl][2]);
							                if (ril<rc && rjl<rc && rkl<rc){

                                                for (int d=0;d<3;d++){
                                                    l=3*nl+d;
								                    //if (j>=i && k>=j && l>=k){
									
                                                        if (start_indx<=fc4_count && fc4_count <= end_indx){
                                                            //if (rank==0) printf(" %d\n", fc4_indx);
                                                            fc4[fc4_indx].i=ni;
                                                            fc4[fc4_indx].j=nj;
                                                            fc4[fc4_indx].k=nk;
                                                            fc4[fc4_indx].l=nl;
                                                            fc4[fc4_indx].a=a;
                                                            fc4[fc4_indx].b=b;
                                                            fc4[fc4_indx].c=c;
                                                            fc4[fc4_indx].d=d;
                                                            //fprintf(fh_debug, " %d %d %d %d %d %d %d %d\n", ni,nj,nk,nl,a,b,c,d);
                                                            fc4_indx++;
                                                        }
                                                        fc4_count++;

                                                    //}
								                }
									        }
									
								        }
                                    }
                                }
				            }
					    }
					}
				}
			}
		}
        
        // Loop over all FC4s on this proc, calculate and store FC
        int nfc4_tol_proc=0;
        for (int w=0; w<nepp[rank]; w++){
            if (rank==0 && w%10==0) printf(" %d\n", w);
            i=fc4[w].i;
            j=fc4[w].j;
            k=fc4[w].k;
            l=fc4[w].l;
            a=fc4[w].a;
            b=fc4[w].b;
            c=fc4[w].c;
            d=fc4[w].d;
            fc4[w].fc = calcD3(i,j,k,l,a,b,c,d);
            //fprintf(fh_debug, " %d %d %d %d %d %d %d %d | %e\n", i,a,j,b,k,c,l,d,fc4[w].fc);
            if (abs(fc4[w].fc) > tol4) nfc4_tol_proc++;
        }

        // Gatherv all fc4 arrays onto the root proc
        MPI_Aint disp[9] = {offsetof(fc4_struct, i), offsetof(fc4_struct, j), offsetof(fc4_struct, k), offsetof(fc4_struct, l), \
                            offsetof(fc4_struct, a), offsetof(fc4_struct, b), offsetof(fc4_struct, c), offsetof(fc4_struct, d), \
                            offsetof(fc4_struct, fc) };
        int blocklen[9] = {1,1,1,1,1,1,1,1,1};
        MPI_Datatype type[9] = {MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_DOUBLE};
        MPI_Datatype fc4_dt; // fc4 datatype
        MPI_Type_create_struct(9, blocklen, disp, type, &fc4_dt);
        MPI_Type_commit(&fc4_dt);
        
        // Gathering the data
        fc4_struct *fc4_all = nullptr;
        
        int nprocs = mc->nprocs;
        
        if (rank == 0){
            fc4_all = new fc4_struct[nfc4];
        }

        // Calculate receive displacements
        int *receive_displacements;
        mem->allocate(receive_displacements,nprocs);
        //int receive_displacements[4];//  = { 0, 0, 1, 3 };
        receive_displacements[0]=0;
        int location = 0;
        for (int w=1; w<nprocs; w++){
            location += nepp[w-1];
            //printf(" %d\n", location);
            receive_displacements[w]=location;
        }

        MPI_Gatherv(fc4,nepp[rank],fc4_dt,fc4_all,nepp,receive_displacements,fc4_dt,0,MPI_COMM_WORLD);

        /* Count the number of FCs above tolerance value on each proc */
        /*
        int nfc4_tol_proc = 0;
        for (int w=0; w<nepp[rank]; w++){
            if (abs(fc4[w].fc) > tol4){
                //printf(" YUP!\n");
                nfc4_tol_proc++;
            }
        }
        */
        //printf(" Proc %d, %d\n", rank,nfc2_tol_proc);

        // Sum all nfc4_tol_proc to the root proc
        int nfc4_tol_total;
        MPI_Reduce(&nfc4_tol_proc, &nfc4_tol_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank==0) printf(" Total FC4 above tol: %d\n", nfc4_tol_total);

        // Print the FCs that are above the tol
        if (rank==0){
            fprintf(fh_fc4, "%d\n",nfc4_tol_total);
            for (int w=0; w<nfc4; w++){
                if (abs(fc4_all[w].fc)>tol4){
                    fprintf(fh_fc4, "%d %d %d %d %d %d %d %d %.20e\n", fc4_all[w].i+1,fc4_all[w].a+1,fc4_all[w].j+1,fc4_all[w].b+1,\
                                                                       fc4_all[w].k+1,fc4_all[w].c+1,fc4_all[w].l+1,fc4_all[w].d+1,fc4_all[w].fc);
                }
            }
        }
        
        mem->deallocate(fc4);
        mem->deallocate(nepp);
        mem->deallocate(receive_displacements);
        if (rank==0) fclose(fh_fc4);

	}

    // Calculate self interactions via ASR
    /*
    double fc2;
    for (i=0; i<natoms; i++){
    }
    */

    //calcD1(0,0,0,0);
    //calcD2(0,0,0,0,0,0);
    //calcD2(0,0,4,0,0,1);

    /*
    int nmodes = 3*natoms;
    for (i=3; i<nmodes; i++){
        printf("i: %d\n",i);
        for (j=i; j<nmodes; j++){
            printf("  j: %d\n", j);
            //printf("  j: %d\n", j);
            for (k=j; k<nmodes; k++){
                //printf("    k: %d\n",k);
                //printf("%d %d %d\n", i,j,k);
                //fprintf(fh_debug, "ijk: %d,%d,%d\n", i,j,k);
                calcD3(i,j,k);
            }
        }
    }
    */

    //mem->deallocate(x);

	if (rank==0) printf(" Finished.\n");

}

/*
Calculate the negative first derivative of force
Inputs: i,j: atom indices (starting from 0)
        a,b: Cartesian direction indices
Outputs: Derivative dFia/dujb
*/

double In::calcD1(int i, int j, int a, int b)
{

    int atomIndices[2] = {i,j};
    int cartIndices[2] = {a,b};

    double disps[1];
    double f1,f2;
    double d1;

    disps[0] = 1.0*delta;
    f1 = calcF(atomIndices,cartIndices,disps,1);
    //printf(" f1: %f\n", f1);
    disps[0] = -1.0*delta;
    f2 = calcF(atomIndices,cartIndices,disps,1);
    //printf(" f2: %f\n", f2);

    d1 = -1.0*(f1-f2)/(2*delta);

    d1 = d1*(0.073499)*(0.529177*0.529177); // Convert to Ryd/bohr^2
    
    //printf(" a,b: %d,%d\n", a,b);
    //printf(" d1: %f\n", d1);
    //if (abs(d1)>tol2) fprintf(fh_fc2, "%d    %d    %d    %d    %.10e\n", i+1,a+1,j+1,b+1,d1);
    //fprintf(fh_fc2, "%d    %d    %d    %d    %.10e\n", i+1,a+1,j+1,b+1,d1);

    return d1;

}

/*
Calculate the negative second derivative of force
Inputs: i,j: atom indices (starting from 0)
        a,b: Cartesian direction indices
Outputs: Derivative dFia/dujb
*/

double In::calcD2(int i, int j, int k, int a, int b, int c)
{

    int atomIndices[3] = {i,j,k};
    int cartIndices[3] = {a,b,c};

    double d2=0;


    double h = delta;
    if (i==j && j==k && a==b && b==c){
        double disps[1];
        double f1,f2,f3;
        disps[0] = 1.0*delta;
        f1 = calcF(atomIndices,cartIndices,disps,1);
        disps[0] = 0.0*delta;
        f2 = calcF(atomIndices,cartIndices,disps,1);
        disps[0] = -1.0*delta;
        f3 = calcF(atomIndices,cartIndices,disps,1);
        //fprintf(fh_fc3, " %f %f %f\n", f1,f2,f3);
        d2 = (-1.0/(delta*delta))*(f1-2*f2+f3);

        //fprintf(fh_fc3, "%d    %d    %d    %d    %d    %d    %.10e\n", i+1,a+1,j+1,b+1,k+1,c+1,d2);
    }
    else{
        double disps[2];
        double f1,f2,f3,f4;
        disps[0] = 0.5*delta;
        disps[1] = 0.5*delta;
        f1 = calcF(atomIndices,cartIndices,disps,2);
        disps[0] = 0.5*delta;
        disps[1] = -0.5*delta;
        f2 = calcF(atomIndices,cartIndices,disps,2);
        disps[0] = -0.5*delta;
        disps[1] = 0.5*delta;
        f3 = calcF(atomIndices,cartIndices,disps,2);
        disps[0] = -0.5*delta;
        disps[1] = -0.5*delta;
        f4 = calcF(atomIndices,cartIndices,disps,2);
        d2 = (-1.0/(h*h))*(f1-f2-f3+f4);        
    }

    d2 = d2*(0.073499)*(0.529177*0.529177*0.529177); // Convert to Ryd/bohr^3   
    //d2 = 2.0*d2; // Need to multiply by 2, for some weird reason. Think I'll need to multiply ifc4 by 3.
    // Only need to multiply by 2 when comparing directly to TEP results? Weird...

    //if (abs(d2)>tol3) fprintf(fh_fc3, "%d    %d    %d    %d    %d    %d    %.10e\n", i+1,a+1,j+1,b+1,k+1,c+1,d2);


    return d2;

}

/*
Calculate the negative third derivative of force
Inputs: i,j: atom indices (starting from 0)
        a,b: Cartesian direction indices
Outputs: Derivative dFia/dujb
*/

double In::calcD3(int i, int j, int k, int l, int a, int b, int c, int d)
{

    int atomIndices[4] = {i,j,k,l};
    int cartIndices[4] = {a,b,c,d};

    double d3=0;


    double h = delta;

    double disps[3];
    double f1,f2,f3,f4,f5,f6,f7,f8;

    disps[0] = 1.0*delta;
    disps[1] = 1.0*delta;
    disps[2] = 1.0*delta;
    f1 = calcF(atomIndices,cartIndices,disps,3);
    disps[0] = 1.0*delta;
    disps[1] = 1.0*delta;
    disps[2] = -1.0*delta;
    f2 = calcF(atomIndices,cartIndices,disps,3);
    disps[0] = 1.0*delta;
    disps[1] = -1.0*delta;
    disps[2] = 1.0*delta;
    f3 = calcF(atomIndices,cartIndices,disps,3);
    disps[0] = 1.0*delta;
    disps[1] = -1.0*delta;
    disps[2] = -1.0*delta;
    f4 = calcF(atomIndices,cartIndices,disps,3);
    disps[0] = -1.0*delta;
    disps[1] = 1.0*delta;
    disps[2] = 1.0*delta;
    f5 = calcF(atomIndices,cartIndices,disps,3);
    disps[0] = -1.0*delta;
    disps[1] = 1.0*delta;
    disps[2] = -1.0*delta;
    f6 = calcF(atomIndices,cartIndices,disps,3);
    disps[0] = -1.0*delta;
    disps[1] = -1.0*delta;
    disps[2] = 1.0*delta;
    f7 = calcF(atomIndices,cartIndices,disps,3);
    disps[0] = -1.0*delta;
    disps[1] = -1.0*delta;
    disps[2] = -1.0*delta;
    f8 = calcF(atomIndices,cartIndices,disps,3);

    d3 = (-1.0/(8*h*h*h))*(f1-f2-f3+f4-f5+f6+f7-f8);        

    d3 = d3*(0.073499)*(0.529177*0.529177*0.529177*0.529177); // Convert to Ryd/bohr^4
    //d3 = 3.0*d3; // Need to multiply by 2, for some weird reason. Think I'll need to multiply ifc4 by 3.

    //if (abs(d3)>tol4) fprintf(fh_fc4, "%d    %d    %d    %d    %d    %d    %d    %d    %.10e\n", i+1,a+1,j+1,b+1,k+1,c+1,l+1,d+1,d3);

    return d3;
}



/*
Calculate the force on atom i due to a bunch of displacements.
Atom i indices are the first indices in atomIndices and cartIndices.
i.e. find atom i via atomIndices[0] and cartIndices[0].
Inputs: atomIndices: indices of atoms to displace
        cartIndices: Cartesian indices to displace (x,y,z)
        ndisps: number of atoms to move when calculating the force
*/

double In::calcF(int atomIndices[], int cartIndices[], double disps[], int ndisps)
{

    //printf("  Mode: %d\n", indices[0]);
    //printf("  Disp: %f\n", disps[0]);
    //printf("num: %d\n", num);
    //printf("modes: %d %d %d\n", indices[0],indices[1],indices[2]);
    //printf("disps: %f %f %f\n", disps[0],disps[1],disps[2]);

    int natoms = lmp->atom->natoms;
    tag = lmp->atom->tag;
    type = lmp->atom->type;
    double *mass = lmp->atom->mass;
    double *pe_p; // Pointer for LAMMPS extraction
    double pe;
    double **x; 
    mem->allocate(x,natoms,3);
    for (int n=0; n<natoms; n++){
        x[n][0]=x0[n][0];
        x[n][1]=x0[n][1];
        x[n][2]=x0[n][2];
    }
    char set_x[50]; // used for sprintf, setting atom x
    char set_y[50]; // used for sprintf, setting atom y
    char set_z[50]; // used for sprintf, setting atom z
    int nchars_x; // number of characters in buffer
    int nchars_y; // number of characters in buffer
    int nchars_z; // number of characters in buffer

    // Displace all necessary atoms 
    int atomIndx,cartIndx;
    for (int d=0; d<ndisps; d++){
        atomIndx = atomIndices[d+1];
        cartIndx = cartIndices[d+1];
        x[atomIndx][cartIndx] += disps[d];
        //fprintf(fh_fc3," x[%d][%d] += %f\n", atomIndx,cartIndx,disps[d]);
    }

    // Set atom positions in LAMMPS
    for (int n=0; n<natoms; n++){
        nchars_x = sprintf(set_x, "set atom %d x %f", n+1, x[n][0]);
        nchars_y = sprintf(set_y, "set atom %d y %f", n+1, x[n][1]);
        nchars_z = sprintf(set_z, "set atom %d z %f", n+1, x[n][2]);

        //std::cout << set_x << std::endl;
        //std::cout << set_y << std::endl;
        //std::cout << set_z << std::endl;
        
        lmp->input->one(set_x);
        lmp->input->one(set_y);
        lmp->input->one(set_z);    
    }

    lmp->input->one("run 0");

    // Get required force on atom i and component a
    
    double f = lmp->atom->f[atomIndices[0]][cartIndices[0]];

    //fprintf(fh_fc3, " f: %f\n", f);

    mem->deallocate(x);

    return f;

}


/*
Calculate 3rd derivative of mode interactions.
Inputs: i,j,k - mode indices (starting from 0).
*/

double In::calcD3_MTEP(int i, int j, int k)
{

    // Check to see what case this is (e.g. iii, iij, ijk, etc.) to apply appropriate finite diff.

    double d3;

    if (i==j && j==k){

        int modeIndices[1] = {i};
        //printf("  iii\n");
        double si=delta/sqrt(freqs[i]);
        //printf("%e\n", delta);
        //printf(" freq: %e\n", freqs[i]);
        //double modeDisps[3] = {si,si,si};
        //double modeDisps[1] = {2.0*si};
        double modeDisps[1];
        modeDisps[0] = 2.0*si;
        double e1 = calcE_MTEP(modeIndices,modeDisps,1); //2
        modeDisps[0] = -1.0*si;
        double e2 = calcE_MTEP(modeIndices,modeDisps,1); //-1
        modeDisps[0] = 1.0*si;
        double e3 = calcE_MTEP(modeIndices,modeDisps,1); //1
        modeDisps[0] = -2.0*si;
        double e4 = calcE_MTEP(modeIndices,modeDisps,1); //-2
        //printf("%f\n", si);
        d3 = (1.0/(2.0*si*si*si))*(e1+(2.0*e2)-(2.0*e3)-e4);
        //printf("-----      d3: %.20f\n", d3);
    }
    else if ( (i==j && i!=k) || (i==k && i!=j) || (j==k && j!=i) ){

        //printf("  iij\n");

        //printf("%d, %d, %d\n", i,j,k);

        int modeIndices[2];
        double modeDisps[2];
    
        double si,sj;

        if (i==j && i!=k){ //iij
            //printf("iij\n");
            modeIndices[0] = i;
            modeIndices[1] = k;
            si=delta/sqrt(freqs[i]);
            sj=delta/sqrt(freqs[k]);
        }
        if (i==k && i!=j){ //iji
            printf("I don't think this case happens!\n");
            double si=delta/sqrt(freqs[i]);
            double sj=delta/sqrt(freqs[j]);
        }
        if (j==k && j!=i){ //ijj
            modeIndices[0] = k;
            modeIndices[1] = i;
            //printf("ijj\n");
            si=delta/sqrt(freqs[k]);
            sj=delta/sqrt(freqs[i]);
        }
        //calcE(modeIndices,3);

        modeDisps[0] = 1.0*si;
        modeDisps[1] = 1.0*sj;
        double e1 = calcE_MTEP(modeIndices,modeDisps,2); //1,1
        modeDisps[0] = 1.0*si;
        modeDisps[1] = -1.0*sj;
        double e2 = calcE_MTEP(modeIndices,modeDisps,2); //1,-1
        modeDisps[0] = 0.0*si;
        modeDisps[1] = 1.0*sj;
        double e3 = calcE_MTEP(modeIndices,modeDisps,2); //0,1
        modeDisps[0] = 0.0*si;
        modeDisps[1] = -1.0*sj;
        double e4 = calcE_MTEP(modeIndices,modeDisps,2); //0,-1
        modeDisps[0] = -1.0*si;
        modeDisps[1] = 1.0*sj;
        double e5 = calcE_MTEP(modeIndices,modeDisps,2); //-1,1
        modeDisps[0] = -1.0*si;
        modeDisps[1] = -1.0*sj;
        double e6 = calcE_MTEP(modeIndices,modeDisps,2); //-1,-1

        d3 = (1.0/(2.0*si*si*sj))*(e1 - e2 - (2.0*e3) + (2.0*e4) + e5 - e6);
        //printf("-----      d3: %.20f\n", d3);

    }

    else if (i!=j && j!=k){
        int modeIndices[3] = {i,j,k};
        //printf("  ijk\n");
        double si=delta/sqrt(freqs[i]);
        double sj=delta/sqrt(freqs[i]);
        double sk=delta/sqrt(freqs[i]);
        //printf(" freq: %f\n", freqs[i]);
        //double modeDisps[3] = {si,si,si};
        //double modeDisps[1] = {2.0*si};
        double modeDisps[3];
        modeDisps[0] = 1.0*si;
        modeDisps[1] = 1.0*sj;
        modeDisps[2] = 1.0*sk;
        double e1 = calcE_MTEP(modeIndices,modeDisps,3); //1,1,1
        modeDisps[0] = 1.0*si;
        modeDisps[1] = 1.0*sj;
        modeDisps[2] = -1.0*sk;
        double e2 = calcE_MTEP(modeIndices,modeDisps,3); //1,1,-1
        modeDisps[0] = 1.0*si;
        modeDisps[1] = -1.0*sj;
        modeDisps[2] = 1.0*sk;
        double e3 = calcE_MTEP(modeIndices,modeDisps,3); //1,-1,1
        modeDisps[0] = 1.0*si;
        modeDisps[1] = -1.0*sj;
        modeDisps[2] = -1.0*sk;
        double e4 = calcE_MTEP(modeIndices,modeDisps,3); //1,-1,-1
        modeDisps[0] = -1.0*si;
        modeDisps[1] = 1.0*sj;
        modeDisps[2] = 1.0*sk;
        double e5 = calcE_MTEP(modeIndices,modeDisps,3); //-1,1,1
        modeDisps[0] = -1.0*si;
        modeDisps[1] = 1.0*sj;
        modeDisps[2] = -1.0*sk;
        double e6 = calcE_MTEP(modeIndices,modeDisps,3); //-1,1,-1
        modeDisps[0] = -1.0*si;
        modeDisps[1] = -1.0*sj;
        modeDisps[2] = 1.0*sk;
        double e7 = calcE_MTEP(modeIndices,modeDisps,3); //-1,-1,1
        modeDisps[0] = -1.0*si;
        modeDisps[1] = -1.0*sj;
        modeDisps[2] = -1.0*sk;
        double e8 = calcE_MTEP(modeIndices,modeDisps,3); //-1,-1,-1
        //printf("%f\n", si);
        d3 = (1.0/(8.0*si*sj*sk))*(e1-e2-e3+e4-e5+e6+e7-e8);
        //printf("-----      d3: %.20f\n", d3);
    }

    //if (d3 > 1e-6){
        //fprintf(fh_mcc3, "%d %d %d %.6e\n", i,j,k,d3);
    //}

}

/*
Calculate system energy corresponding to displacement of modes in input array.
Input: indices - array of indices of modes to displace.
       
       num - number of modes to displace.
*/
double In::calcE_MTEP(int indices[], double disps[], int num){

    //printf("  Mode: %d\n", indices[0]);
    //printf("  Disp: %f\n", disps[0]);
    //printf("num: %d\n", num);
    //printf("modes: %d %d %d\n", indices[0],indices[1],indices[2]);
    //printf("disps: %f %f %f\n", disps[0],disps[1],disps[2]);

    int natoms = lmp->atom->natoms;
    tag = lmp->atom->tag;
    type = lmp->atom->type;
    double *mass = lmp->atom->mass;
    double *pe_p; // Pointer for LAMMPS extraction
    double pe;
    double **x; 
    mem->allocate(x,natoms,3);
    for (int n=0; n<natoms; n++){
        x[n][0]=x0[n][0];
        x[n][1]=x0[n][1];
        x[n][2]=x0[n][2];
    }

    char set_x[50]; // used for sprintf, setting atom x
    char set_y[50]; // used for sprintf, setting atom y
    char set_z[50]; // used for sprintf, setting atom z
    int nchars_x; // number of characters in buffer
    int nchars_y; // number of characters in buffer
    int nchars_z; // number of characters in buffer
    
    // Loop over modes 

    double mass_convert;

    for (int i=0; i<num; i++){

        int eindx = 0;
        for (int n=0; n<natoms; n++){
            /*
            x[eindx+0] = x0[n][0]+disps[i]*emat[0+eindx][i];
            x[eindx+1] = x0[n][1]+disps[i]*emat[1+eindx][i];
            x[eindx+2] = x0[n][2]+disps[i]*emat[2+eindx][i];
            */
            //fprintf(fh_xyz, "%d %f %f %f\n", lmp->atom->type[n], x[eindx+0],x[eindx+1],x[eindx+2]);
            //fprintf(fh_debug,"%f %f %f\n",x[eindx+0],x[eindx+1],x[eindx+2]);
            //fprintf(fh_debug,"%d, %d\n",n, tag[n]);

            mass_convert = (mass[type[n]]/(6.02214076e23))*1e-3;
            //printf(" mass_convert: %e\n", mass_convert);

            x[n][0] += (1.0/sqrt(mass_convert))*disps[i]*emat[0+eindx][indices[i]];
            x[n][1] += (1.0/sqrt(mass_convert))*disps[i]*emat[1+eindx][indices[i]];
            x[n][2] += (1.0/sqrt(mass_convert))*disps[i]*emat[2+eindx][indices[i]];

            //fprintf(fh_debug, "disp: %f\n", disps[i]);
            
            /*
            nchars_x = sprintf(set_x, "set atom %d x %f", n+1, x0[n][0]+disps[i]*emat[0+eindx][indices[i]]);
            nchars_y = sprintf(set_y, "set atom %d y %f", n+1, x0[n][1]+disps[i]*emat[1+eindx][indices[i]]);
            nchars_z = sprintf(set_z, "set atom %d z %f", n+1, x0[n][2]+disps[i]*emat[2+eindx][indices[i]]);
            */

            nchars_x = sprintf(set_x, "set atom %d x %f", n+1, x[n][0]);
            nchars_y = sprintf(set_y, "set atom %d y %f", n+1, x[n][1]);
            nchars_z = sprintf(set_z, "set atom %d z %f", n+1, x[n][2]);
            
            lmp->input->one(set_x);
            lmp->input->one(set_y);
            lmp->input->one(set_z);
            

            eindx+=3;
        }

        
    }

    //lmp->input->one("set atom 1 x 10.944");
    lmp->input->one("run 0");

    // Extract the potential energy

    pe_p = (double *) lammps_extract_compute(lmp, "P", 0, 0); // style = 0 for global data, type = 0 for scalar quantity
    pe = pe_p[0];
    //printf("     pe: %.15f\n", pe);
    // Convert to Joules
    //pe = pe*1.60218e-19;

    mem->deallocate(x);

    return pe;
    
}

/*
Calculate radial distance between pairs of atoms
*/
double In::calcRij(double xi, double yi, double zi, double xj, double yj, double zj){

    double boxlo_x, boxlo_y, boxlo_z, boxhi_x, boxhi_y, boxhi_z;
    boxlo_x = lmp->domain->boxlo[0];
    boxlo_y = lmp->domain->boxlo[1];
    boxlo_z = lmp->domain->boxlo[2];
    boxhi_x = lmp->domain->boxhi[0];
    boxhi_y = lmp->domain->boxhi[1];
    boxhi_z = lmp->domain->boxhi[2];

    double boxinvx = 1.0/boxhi_x;
    double boxinvy = 1.0/boxhi_y;
    double boxinvz = 1.0/boxhi_z;

    double xij = xi - xj;
    double yij = yi - yj;
    double zij = zi - zj;

    int nearintx = std::round(xij*boxinvx);
    int nearinty = std::round(yij*boxinvy);
    int nearintz = std::round(zij*boxinvz);

    xij = xij - boxhi_x*nearintx;
    yij = yij - boxhi_y*nearinty;
    zij = zij - boxhi_z*nearintz;

    double rij = sqrt(xij*xij + yij*yij + zij*zij);

    return rij;

}
