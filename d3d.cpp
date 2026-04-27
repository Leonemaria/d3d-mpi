#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h> 
#include <cstdlib> 
#include <cmath>
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
    MPI_Init(&argc, &argv);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
// Header
    std::cout << "*** d3d START ***\n";
// input of computational ambient data
    global glb;
    std::string caseName, s; int Nth;
    int nameSize;
    if (world_rank==0)
    {
        std::cout << "sono zero\n";
    }
    s="./ambient.dat"; std::ifstream ambientFile(s); chk(ambientFile,s); // ambient data
    ambientFile >> Nth; skipLine(ambientFile, 1); //reading the number of threads
    std::cout << "N threads=" << Nth << std::endl;
    omp_set_num_threads(Nth);
    ambientFile >> caseName; //reading the case name
    nameSize=caseName.size();
    ambientFile.close();
    std::cout << "Running case: " << caseName << std::endl;
//    }
    omp_set_num_threads(2);
    MPI_Bcast(&nameSize,1,MPI_INT,0,MPI_COMM_WORLD);
    caseName.resize(nameSize);
    MPI_Bcast(&caseName[0],nameSize,MPI_CHAR,0,MPI_COMM_WORLD);
    // imput of physical parametrs
    s="./"+caseName+"/input/physics.dat"; std::ifstream physicsFile(s); chk(physicsFile,s); // physics data
    readPhysics(physicsFile,glb);
    physicsFile.close();
// input of run data
    s="./"+caseName+"/input/run.dat"; std::ifstream runDataFile(s); chk(runDataFile,s); // run setting data input file
    readRun(runDataFile,glb);
    runDataFile.close();
    std::cout << "src=" << glb.sch[3] << std::endl;
    computationalElement cc(glb);
    std::cout << "dt=" << glb.dt << std::endl;
    std::cout << "LES model=" << glb.sch[1] << std::endl;
//  input of boundary conditions
    int nBC;
    s="./"+caseName+"/input/BC.dat"; std::ifstream inputFileBC(s); chk(inputFileBC,s); // boundary conditions definition input file
    inputFileBC >> nBC; nBC++; skipLine(inputFileBC, 1); // reading number of boundary conditions
    std::cout << "nBC=" << nBC << std::endl;
    boundaryCondition BC[nBC]; for (int i=1; i<nBC; i++) {std::cout << "BC N." << i << std::endl; BC[i].input(inputFileBC);}   
//***************************************
//  input of geometry (points, elements, links and boundary conditions)
    long nNodes, nCells;
    s="./"+caseName+"/input/mesh.dat"; std::ifstream inputFileMesh(s); chk(inputFileMesh,s); // mesh input files
    s="./"+caseName+"/input/link.dat"; std::ifstream inputFileLink(s); chk(inputFileLink,s); // link input files
    inputFileMesh >> nNodes >> nCells; // reading number of grid points and number of physical cells
    std::cout << "nCells=" << nCells << std::endl;
    vector3D* xN=new vector3D[nNodes]; // grid point array    
    physicalElement* e=new physicalElement[nCells]; // physical cells array    
    double volume=readMesh(inputFileMesh,inputFileLink,nNodes,xN,nCells,e,glb,&cc);
    inputFileMesh.close(); inputFileLink.close();
    std::cout << "volume=" << volume << std::endl;
// initialization    
    initialConditions(caseName,nCells,e,glb);
    matrix H;
// start time marching simulation       
    double residual, maxCFL, eCFL;
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
    double t_start = omp_get_wtime();
    bool dmpR, dmpH=false;
    if (glb.ctr[0]==0) {printOut(caseName,e,nCells,0,glb);}
// start iterations in time
    for (int i=glb.ctr[0]+1; i<(glb.ctr[1]+1); i++)
    {
        residual=0.; maxCFL=0.;
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
            printOut(caseName,e,nCells,i,glb);
        }
        if (dmpR)
        {
            for (long iC=0; iC<nCells; iC++)
            {
                residual+=std::abs(e[iC].getResidual());
                eCFL=e[iC].CFL(glb.dt);
                if (eCFL>maxCFL) {maxCFL=eCFL;}
            }
            residual=residual/volume;
            outputFileRes << i << " " << residual << " " << maxCFL << std::endl;
            std::cout << "Res.=" << residual << "  CFL=" << maxCFL << std::endl;
            dmpR=false;
        }
        if (i % glb.ctr[4]==0)
        {
            H=e[0].getHist();
            for (long iC=1; iC<nCells; iC++)
            {
               H+=e[iC].getHist();
            }
            H=H/volume;
            outputFileHist << i*glb.dt;
            for (int j=0; j<H.nC(); j++)
            {
                outputFileHist << " " << H.get(j);
            }
            outputFileHist << std::endl;
        }
    }
    delete[] e; delete[] xN;
    double t_end = omp_get_wtime();
    std::cout << "Elapsed time: " << t_end - t_start << " s\n";
    //
    MPI_Finalize();
    return 0;
}
