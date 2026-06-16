#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h> 
#include <cstdlib> 
#include <cmath>
#include <vector>
#include "matrix.h"
#include "computationalElement.h"
#include "physicalElement.h"
#include "vectors.h"
#include "utilities.h"
#include "inputOutput.h"
#include "cases.h"
#include "global.h"
//
int main(int argc, char** argv)
{
//    int provided;
//    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Init(&argc, &argv);
    int w_size;
    MPI_Comm_size(MPI_COMM_WORLD, &w_size);
    int w_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &w_rank);
// Header
    std::cout << "*** d3d START ***\n";
// input of computational ambient data
    global glb;
    std::string caseName, s; int Nth;
    int nameSize;
    if (w_rank==0)
    {
        s="./ambient.dat"; std::ifstream ambientFile(s); chk(ambientFile,s); // ambient data
        ambientFile >> caseName; //reading the case name
        nameSize=caseName.size();
        ambientFile.close();
        std::cout << "Running case: " << caseName << std::endl;
    }
    MPI_Bcast(&Nth,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&nameSize,1,MPI_INT,0,MPI_COMM_WORLD);
    caseName.resize(nameSize);
    MPI_Bcast(&caseName[0],nameSize,MPI_CHAR,0,MPI_COMM_WORLD);
    // imput of physical parametrs
    if (w_rank==0)
    {
        s="./"+caseName+"/input/physics.dat"; std::ifstream physicsFile(s); chk(physicsFile,s); // physics data
        readPhysics(physicsFile,glb);
        physicsFile.close();
    }
    MPI_Bcast(&glb.phy,6,MPI_DOUBLE,0,MPI_COMM_WORLD);
// input of run data
    if (w_rank==0)
    {
        s="./"+caseName+"/input/run.dat"; std::ifstream runDataFile(s); chk(runDataFile,s); // run setting data input file
        readRun(runDataFile,glb);
        runDataFile.close();
        std::cout << "pol. order N=" << glb.sch[0] << std::endl;
        std::cout << "Num. of Modes=" << nModes(glb.sch[0]) << std::endl;
        std::cout << "src=" << glb.sch[3] << std::endl;
        std::cout << "dt=" << glb.dt << std::endl;
        std::cout << "LES model=" << glb.sch[1] << std::endl;
    }
    MPI_Bcast(&glb.dt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(glb.ctr,5,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(glb.sch,4,MPI_INT,0,MPI_COMM_WORLD);
    computationalElement cc(glb);
//  input of boundary conditions
    int nBC; std::string process=std::to_string(w_rank);
    s="./"+caseName+"/input/BC_pr"+process+".dat"; std::ifstream inputFileBC(s); chk(inputFileBC,s); // boundary conditions definition input file
    inputFileBC >> nBC; nBC++; skipLine(inputFileBC, 1); // reading number of boundary conditions
    boundaryCondition BC[nBC]; for (int i=1; i<nBC; i++) {BC[i].input(inputFileBC);}   
//***************************************
//  input of geometry (points, elements, links and boundary conditions)
    long nNodes, nCells, totCells=0;
    s="./"+caseName+"/input/mesh_pr"+process+".dat"; std::ifstream inputFileMesh(s); chk(inputFileMesh,s); // mesh input files
    s="./"+caseName+"/input/link_pr"+process+".dat"; std::ifstream inputFileLink(s); chk(inputFileLink,s); // link input files
    inputFileMesh >> nNodes >> nCells; // reading number of grid points and number of physical cells
    vector3D* xN=new vector3D[nNodes]; // grid point array    
    physicalElement* e=new physicalElement[nCells]; // physical cells array    
    double volume=readMesh(inputFileMesh,inputFileLink,nNodes,xN,nCells,e,glb,&cc,w_rank);
    inputFileMesh.close(); inputFileLink.close();
    double totVolume=0.;
    MPI_Reduce(&volume, &totVolume, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // initialization    
    MPI_Reduce(&nCells, &totCells, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); // initialization    
    if (w_rank==0) {std::cout << "Volume=" << totVolume << std::endl; std::cout << "Num. of Cells=" << totCells << std::endl;}
// preparation for mpi messaging
    int nJnPr=0; // number of joined processes
    std::vector<int>jnPr; // contains the ranks of connected processes
    std::vector<std::vector<int>> jnFaces; // contains the nJnPr vectors containing the indices of the joined cells and faces
    bool newJ; intMatrix jn;
    for (long iC=0; iC<nCells; iC++)
    {
        jn=e[iC].getJoin();
        for (int iS=0; iS<4; iS++)
        {
            if (jn.get(iS,0)!=w_rank)
            {
                newJ=true;           
                for (int i=0; i<nJnPr; i++)
                {
                    if (jnPr[i]==jn.get(iS,0))
                    {
                        newJ=false;
                        jnFaces[i].push_back(jn.get(iS,1));
                        jnFaces[i].push_back(jn.get(iS,2));
//                        if (w_rank<jnPr[i]) {e[iC].setJoin(iS,1,jnFaces[i].size()-1);}                        
                    }
                }
                if (newJ)
                {
                    nJnPr++; jnPr.push_back(jn.get(iS,0));
                    jnFaces.push_back({jn.get(iS,1)});
                    jnFaces[nJnPr-1].push_back(jn.get(iS,2));
//                    if (w_rank<jn.get(iS,0)) {e[iC].setJoin(iS,1,0);} 
                }
            }
        }
    }
    std::cout << "my rank=" << w_rank << " NjnFaces=" << jnFaces[0].size() << std::endl;
    std::vector <MPI_Request> rqs(nJnPr);
    for (int i=0; i<nJnPr; i++)
    {
        if (w_rank<jnPr[i])
        {
            MPI_Isend(jnFaces[i].data(),jnFaces[i].size(),MPI_INT,jnPr[i],1,MPI_COMM_WORLD,&rqs[i]);
        }
        else
        {
            MPI_Irecv(jnFaces[i].data(),jnFaces[i].size(),MPI_INT,jnPr[i],1,MPI_COMM_WORLD,&rqs[i]);
        }        
    }
    MPI_Waitall(nJnPr,rqs.data(),MPI_STATUSES_IGNORE);   
//
    initialConditions(caseName,nCells,e,glb,w_rank);
    matrix H;
// start time marching simulation       
    std::ofstream outputFileRes,outputFileHist;
    if (glb.ctr[0]==0)
    {
        outputFileRes.open("./"+caseName+"/output/residual.dat",std::ios::out);
        outputFileHist.open("./"+caseName+"/output/history.dat",std::ios::out);
    }
    else
    {
        outputFileRes.open("./"+caseName+"/output/residual.dat",std::ios::app);
        outputFileHist.open("./"+caseName+"/output/history.dat",std::ios::app);
    }
    outputFileHist << std::setprecision(12);
    double t_start; if (w_rank==0) {t_start=omp_get_wtime();}
    bool dmpR, dmpH=false;
    if (glb.ctr[0]==0) {printOut(caseName,e,nCells,0,glb,w_rank);}
// start iterations in time
    for (int i=glb.ctr[0]+1; i<(glb.ctr[1]+1); i++)
    {
        for (int ii=0; ii<5; ii++)
        {
           if ((ii==4)&&(i % glb.ctr[3]==0)){dmpR=true;} else {dmpR=false;}
           if ((ii==0)&&(i % glb.ctr[4]==0)){dmpH=true;} else {dmpH=false;}
// compute (for all elements) conservative and auxiliary (primitive) variables on side quadrature points
#pragma omp parallel for schedule(static)
            for (long iC=0; iC<nCells; iC++)
            {
                e[iC].step_0(BC); // computes conservative and auxiliary (primitive) variables on side quadrature points
            }
// computes (far all elements) the auxiliary variables gradients and physical fluxes
#pragma omp parallel for schedule(static)
            for (long iC=0; iC<nCells; iC++)
            {
                e[iC].step_I(caseName,e,BC,&dmpH); // computes the auxiliary variable gradients and physical fluxes on all quadrature points
            }
// computes numerical fluxes and advances the solution
#pragma omp parallel for schedule(static)
            for (long iC=0; iC<nCells; iC++)
            {
                {e[iC].step_II(glb.dt,ii,e,dmpR);}
            }
        }
        if ((i % glb.ctr[2]==0)||(i==glb.ctr[1]))
        {
            printOut(caseName,e,nCells,i,glb,w_rank);
        }
        if (dmpR)
        {
            double residual=0., totResidual, eCFL, maxCFL=0., maxMaxCFL;
            for (long iC=0; iC<nCells; iC++)
            {
                residual+=std::abs(e[iC].getResidual());
                eCFL=e[iC].CFL(glb.dt);
                if (eCFL>maxCFL) {maxCFL=eCFL;}
            }
            MPI_Reduce(&residual,&totResidual,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
            MPI_Reduce(&maxCFL,&maxMaxCFL,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
            if (w_rank==0)
            {
                totResidual=totResidual/totVolume;
                outputFileRes << i << " " << totResidual << " " << maxCFL << std::endl;
                std::cout << "Res.=" << totResidual << "  CFL=" << maxMaxCFL << std::endl;
            }
            dmpR=false;
        }
        if (i % glb.ctr[4]==0)
        {
            H=e[0].getHist();
            for (long iC=1; iC<nCells; iC++)
            {
               H+=e[iC].getHist();
            }
            double h[H.nC()], glH[H.nC()]={0.};
            for (int j=0; j<H.nC(); j++) {h[j]=H.get(j);}
            MPI_Reduce(h,glH,H.nC(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
            if(w_rank==0)
            {
                outputFileHist << i*glb.dt;
                for (int j=0; j<H.nC(); j++)
                {
                    outputFileHist << " " << glH[j]/totVolume;
                }
                outputFileHist << std::endl;
            }
        }
    }
    delete[] e; delete[] xN;
    if (w_rank==0) {double t_end = omp_get_wtime(); std::cout << "Elapsed time: " << t_end - t_start << " s\n";}
    //
    MPI_Finalize();
    return 0;
}
