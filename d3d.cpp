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
    int commSize;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
// Header
    std::cout << "*** d3d START ***\n";
// input of computational ambient data
    global glb;
    std::string caseName, s; int Nth;
    int nameSize;
    if (myRank==0)
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
    if (myRank==0)
    {
        s="./"+caseName+"/input/physics.dat"; std::ifstream physicsFile(s); chk(physicsFile,s); // physics data
        readPhysics(physicsFile,glb);
        physicsFile.close();
    }
    MPI_Bcast(&glb.phy,6,MPI_DOUBLE,0,MPI_COMM_WORLD);
// input of run data
    if (myRank==0)
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
    int nBC; std::string process=std::to_string(myRank);
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
    double volume=readMesh(inputFileMesh,inputFileLink,nNodes,xN,nCells,e,glb,&cc,myRank);
    inputFileMesh.close(); inputFileLink.close();
    double totVolume=0.;
    MPI_Reduce(&volume, &totVolume, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // initialization    
    MPI_Reduce(&nCells, &totCells, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); // initialization    
    if (myRank==0) {std::cout << "Volume=" << totVolume << std::endl; std::cout << "Num. of Cells=" << totCells << std::endl;}
    MPI_Barrier(MPI_COMM_WORLD);   
// preparation for mpi messaging
    int nJnPr=0; // number of joined processes
    std::vector<int>jnPr; // contains the ranks of connected processes
    std::vector<std::vector<int>> jnFaces; // contains the nJnPr vectors containing the indices of the joined cells and faces
    matrix* qSnd=nullptr; matrix* qRcv=nullptr;
    matrix* qASnd=nullptr; matrix* qARcv=nullptr;
    matrix* fSnd=nullptr; matrix* fRcv=nullptr;
    bool newJ; intMatrix jn;
    int Npq2;
    for (long iC=0; iC<nCells; iC++)
    {
        jn=e[iC].getJoin();
        Npq2=e[iC].nQuadPoints();
        for (int iS=0; iS<4; iS++)
        {
            if (jn.get(iS,0)!=myRank) // if the cell side is connected to another process
            {
                newJ=true;           
                for (int i=0; i<nJnPr; i++)
                {
                    if (jnPr[i]==jn.get(iS,0)) // if the joined process is already in the list of joined processes
                    {
                        newJ=false;
                        jnFaces[i].push_back(jn.get(iS,1)); // adds the index of the joined cell to the linking data vector of the corresponding process
                        jnFaces[i].push_back(jn.get(iS,2)); // adds the index of the joined face to the linking data vector of the corresponding process
                        // if the present process rank is lower than the joined one, the second element in the face join vector (the one that normally indicates
                        // the joined cell) is substituted with the index of the joined process in jnPr vector, while the third element is substituted with the
                        // index of the face in the "interface wall" (the set of faces joined between the two processes)
                        if (myRank<jnPr[i]) {e[iC].setJoin(iS,1,i);}                    
                        if (myRank<jnPr[i]) {e[iC].setJoin(iS,2,jnFaces[i].size()/2-1);}                    
                    }
                }
                if (newJ) // if the joined process was not in the list of joined processes
                {
                    nJnPr++; jnPr.push_back(jn.get(iS,0));
                    jnFaces.push_back({jn.get(iS,1),jn.get(iS,2)}); // adds a new linking indices vector with the indices of the joined cell and face
                    if (myRank<jn.get(iS,0)) {e[iC].setJoin(iS,1,nJnPr-1);} 
                    if (myRank<jn.get(iS,0)) {e[iC].setJoin(iS,2,0);} 
                }
            }
        }
    }
    if (nJnPr>0)
    {
        qSnd=new matrix[nJnPr]; qRcv=new matrix[nJnPr];
        qASnd=new matrix[nJnPr]; qARcv=new matrix[nJnPr];
        fSnd=new matrix[nJnPr]; fRcv=new matrix[nJnPr];
    }  
    for (int i=0; i<nJnPr; i++) // allocate the buffers for data exchange between the processes
    {
        qSnd[i].dim((jnFaces[i].size()/2)*Npq2,5);
        qRcv[i].dim((jnFaces[i].size()/2)*Npq2,5);
        qASnd[i].dim((jnFaces[i].size()/2)*Npq2,4);
        qARcv[i].dim((jnFaces[i].size()/2)*Npq2,4);
        fSnd[i].dim((jnFaces[i].size()/2)*Npq2,5);
        fRcv[i].dim((jnFaces[i].size()/2)*Npq2,5);
    }
    MPI_Request rqs[4*nJnPr];
    for (int i=0; i<nJnPr; i++) // exchange the linking data between the processes
    {
        if (myRank<jnPr[i]) // if the present process rank is lower than the joined one
        {
            MPI_Isend(jnFaces[i].data(),jnFaces[i].size(),MPI_INT,jnPr[i],1,MPI_COMM_WORLD,&rqs[i]); // the linking indices are sent
        }
        else // if the present process rank is larger than the joined one
        {
            MPI_Irecv(jnFaces[i].data(),jnFaces[i].size(),MPI_INT,jnPr[i],1,MPI_COMM_WORLD,&rqs[i]); // the linking indices are received
        }
    }
    MPI_Waitall(nJnPr,rqs,MPI_STATUSES_IGNORE);
    MPI_Barrier(MPI_COMM_WORLD);
//  
    for (int i=0; i<nJnPr; i++)
    {
        if (myRank>jnPr[i]) // if the present process rank is larger than the joined one
        {
            for (int j=0; j<jnFaces[i].size();j+=2)
            {
                e[jnFaces[i][j]].setJoin(jnFaces[i][j+1],1,i); // the second element in the face join vector of the joined cell
                e[jnFaces[i][j]].setJoin(jnFaces[i][j+1],2,j/2); // the second element in the face join vector of the joined cell
                // is substituted with the index of the face in the interface wall
            }
        }
    }
    std::vector<std::vector<int>>().swap(jnFaces);
//
    initialConditions(caseName,nCells,e,glb,myRank);
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
    double t_start; if (myRank==0) {t_start=omp_get_wtime();}
    bool dmpR, dmpH=false;
    if (glb.ctr[0]==0) {printOut(caseName,e,nCells,0,glb,myRank);}
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
                e[iC].step_0(BC,myRank,qSnd,qASnd); // computes conservative and auxiliary (primitive) variables on side quadrature points
            }
// computes (far all elements) the auxiliary variables gradients and physical fluxes
    for (int i=0; i<nJnPr; i++) // exchange the data between the processes
    {
        MPI_Isend(qSnd[i].data(),qSnd[i].size(),MPI_DOUBLE,jnPr[i],1,MPI_COMM_WORLD,&rqs[i*4]); // the linking indices are sent
        MPI_Isend(qASnd[i].data(),qASnd[i].size(),MPI_DOUBLE,jnPr[i],2,MPI_COMM_WORLD,&rqs[i*4+1]); // the linking indices are sent
        MPI_Irecv(qRcv[i].data(),qRcv[i].size(),MPI_DOUBLE,jnPr[i],1,MPI_COMM_WORLD,&rqs[i*4+2]); // the linking indices are received
        MPI_Irecv(qARcv[i].data(),qARcv[i].size(),MPI_DOUBLE,jnPr[i],2,MPI_COMM_WORLD,&rqs[i*4+3]); // the linking indices are received
    }
    MPI_Waitall(4*nJnPr,rqs,MPI_STATUSES_IGNORE);
    MPI_Barrier(MPI_COMM_WORLD);   

#pragma omp parallel for schedule(static)
            for (long iC=0; iC<nCells; iC++)
            {
                e[iC].step_I(caseName,e,BC,&dmpH,myRank,qARcv,fSnd); // computes the auxiliary variable gradients and physical fluxes on all quadrature points
            }
// computes numerical fluxes and advances the solution

    for (int i=0; i<nJnPr; i++) // exchange the data between the processes
    {
        MPI_Isend(fSnd[i].data(),fSnd[i].size(),MPI_DOUBLE,jnPr[i],1,MPI_COMM_WORLD,&rqs[i*2]); // the linking indices are sent
        MPI_Irecv(fRcv[i].data(),fRcv[i].size(),MPI_DOUBLE,jnPr[i],1,MPI_COMM_WORLD,&rqs[i*2+1]); // the linking indices are received
    }
    MPI_Waitall(2*nJnPr,rqs,MPI_STATUSES_IGNORE);
    MPI_Barrier(MPI_COMM_WORLD);

#pragma omp parallel for schedule(static)
            for (long iC=0; iC<nCells; iC++)
            {
                e[iC].step_II(glb.dt,ii,e,dmpR,myRank,qRcv,fRcv);
            }
        }
        if ((i % glb.ctr[2]==0)||(i==glb.ctr[1]))
        {
            printOut(caseName,e,nCells,i,glb,myRank);
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
            if (myRank==0)
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
            if(myRank==0)
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
    delete[] e; e=nullptr; delete[] xN; xN=nullptr;
    delete[] qSnd; qSnd=nullptr; delete[] qRcv; qRcv=nullptr;
    delete[] qASnd; qASnd=nullptr; delete[] qARcv; qARcv=nullptr;
    if (myRank==0) {double t_end = omp_get_wtime(); std::cout << "Elapsed time: " << t_end - t_start << " s\n";}
    //
    MPI_Finalize();
    return 0;
}
