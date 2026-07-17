//inputOutput.cpp
#include "inputOutput.h"
void initialConditions(std::string caseName, long nCells, physicalElement e[], const global& glb, int rank)  // set initial conditions
{
    if (glb.ctr[0]==0)
    {
        for (long iC=0; iC<nCells; iC++)
        {
            e[iC].setIniCond(caseName);
        }
    }
    else
    {
        double aM; std::string st=std::to_string(glb.ctr[0]); std::string process=std::to_string(rank);
        std::string s="./"+caseName+"/output/result_pr"+process+"_it"+st+".dat"; std::ifstream fin(s); chk(fin,s); 
	skipLine(fin, 1); int Nm;
        for (long iC=0; iC<nCells; iC++)
        {
            skipLine(fin, 1);
            Nm=nModes(glb.sch[0]);
            for (int iM=0; iM<Nm; iM++)
            {
                for (int eq=0; eq<nEq; eq++)
                {
                    fin >> aM;
                    e[iC].setAM(iM,eq,aM);
                }
                skipLine(fin,1);
            }
        }
    }
}  
void readPhysics(std::ifstream &fin, global& glb)
{
    fin >> glb.phy[0]; skipLine(fin, 1); // reference Mach Number
    fin >> glb.phy[1]; skipLine(fin, 1); // reference Reynolds Number
    fin >> glb.phy[2]; skipLine(fin, 1); // reference Froude Number
    fin >> glb.phy[3]; skipLine(fin, 1); // reference Prandtl Number
    fin >> glb.phy[4]; skipLine(fin, 1);// specific heats ratio
    fin >> glb.phy[5]; skipLine(fin, 1);  // Sutherland temperature ratio
}
void readRun(std::ifstream &fin, global& glb)
{
    fin >> glb.sch[0]; skipLine(fin, 1);        // polynomials max order
    fin >> glb.sch[1]; skipLine(fin, 1);      // kind of LES model
    fin >> glb.sch[2]; skipLine(fin, 1);      // kind of convective intercell numerical flux
    fin >> glb.ctr[0]; skipLine(fin, 1);// initial number of time steps
    fin >> glb.ctr[1]; skipLine(fin, 1);  // final number of time steps
    fin >> glb.ctr[2]; skipLine(fin, 1);  // solution damping period
    fin >> glb.ctr[3]; skipLine(fin, 1);  // residual damping period
    fin >> glb.ctr[4]; skipLine(fin, 1); // integrated quantities damping period
    fin >> glb.dt; skipLine(fin, 1);       // time step
    fin >> glb.sch[3];                        // source terms flag

}
double readMesh(std::ifstream &fin1, std::ifstream &fin2, long nNodes, vector3D xN[], long nCells, physicalElement e[], const global& glb, computationalElement *cc, int rank)
{
    double V=0.;
    long iN=-1, iC=-1, iVer[4];
    while (iN<nNodes-1) // reading of the grid point list (the first part of mesh file)
    {
        fin1 >> iN >> xN[iN][0] >> xN[iN][1] >> xN[iN][2]; //grid point reading from mesh file (each point is labeled by a sequential number and is defined by three coordinates)
    }
    longMatrix l(4,4); l.zero();
    long i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15;
    //reading link data:
    while (iC<nCells-1) // reading of the cell list (second part of the mesh file) and linkings (link file)
    {
        fin1 >> iC >> iVer[0] >> iVer[1] >> iVer[2] >> iVer[3];    //reading cell data:
        // each raw of cell list has the sequential index of the cell followed by the indices of the points corresponding to its four vertices (in order from vertex 0 to 3)
        fin2 >> i0 >> i1 >> i2 >> i3 >> i4 >> i5 >> i6 >> i7 >> i8 >> i9 >> i10 >> i11 >> i12 >> i13 >> i14 >> i15; //reading link data:
        // there are four foursomes of data, one for each face (from face 0 to 3), where the first data is the rank of the linked process,
		// the second data is the index of the linked cell, the third is the index of its linked face, and the last one is the index of the face vertex linked to the 0 vertex
		// (unless the celle is connected to itself at the same face: in this case the last index is the BC)
        l.set(0,0,i0); l.set(0,1,i1); l.set(0,2,i2); l.set(0,3,i3); // setting link matrix
        l.set(1,0,i4); l.set(1,1,i5); l.set(1,2,i6); l.set(1,3,i7); // setting link matrix
        l.set(2,0,i8); l.set(2,1,i9); l.set(2,2,i10); l.set(2,3,i11); // setting link matrix
        l.set(3,0,i12); l.set(3,1,i13); l.set(3,2,i14); l.set(3,3,i15); // setting link matrix
        // the link matrix is a 4x4 matrix with a row for each cell face where the first element is the indices of the connected processor
        // the second one is the index of the linked cell, the third one is the indices of the linked side
        // and the fourth is the index of the point linked to 0 point (or the boundary contition)
        e[iC].init(iC,cc,&xN[0],iVer,l,glb); // physical cell initialization
        for (int iS=0; iS<4; iS++)
        {
            if ((l.get(iS,0)==rank)&(l.get(iS,1)==iC)&(l.get(iS,2)==iS)) {e[iC].setBC(iS);}            
        }
        V+=e[iC].Jacobian();
    }
    return 4.*V/3.;  
}
void printOut(std::string caseName, physicalElement e[], long nCells, int it, const global& glb, int rank)
{
    void subCellNodes(std::ofstream &fl, int p, int n0, int n1, int n2, int n3);
    std::string iterations=std::to_string(it); std::string process=std::to_string(rank);
    int Nm=nModes(glb.sch[0]);
    std::ofstream outputFileSol("./"+caseName+"/output/result_pr"+process+"_it"+iterations+".dat");
    outputFileSol << std::setprecision(16);
    outputFileSol << nCells*Nm << " " << glb.dt*it << std::endl;
    for (long i=0; i<nCells; i++)
    {
        outputFileSol << glb.sch[0] << std::endl;
        for (int k=0; k<Nm; k++)
        {
            outputFileSol << e[i].getAM(k,0) << " " << e[i].getAM(k,1) << " " << e[i].getAM(k,2) << " " << e[i].getAM(k,3) << " " << e[i].getAM(k,4) << std::endl;            
        }
    }
}
void subCellNodes(std::ofstream &fl, int p, int n0, int n1, int n2, int n3)
{
    fl << 4 << " " << p+n0 << " " << p+n1 << " " << p+n2 << " " << p+n3 << std::endl;
}
