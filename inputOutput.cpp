//inputOutput.cpp
#include "inputOutput.h"
void initialConditions(std::string caseName, int nCells, physicalElement e[], global* g, boundaryCondition BC[], int rank)  // set initial conditions
{
    if ((*g).ctr[0]==0)
    {
        for (int iC=0; iC<nCells; iC++)
        {
            e[iC].setIniCond(caseName,BC);
        }
    }
    else
    {
        double* aM; std::string st=std::to_string((*g).ctr[0]); std::string process=std::to_string(rank);
        std::string s="./"+caseName+"/output/result_pr"+process+"_it"+st+".dat"; std::ifstream fin(s,std::ios::binary); chk(fin,s); 
        fin.seekg(sizeof(int)+sizeof(double), std::ios::cur);
	int N;
        for (int iC=0; iC<nCells; iC++)
        {
            fin.read(reinterpret_cast<char*>(&N), sizeof(int));
            aM=e[iC].toAM();
            for (int iA=0; iA<nModes(N)*nEq; iA++)
            {
                    fin.read(reinterpret_cast<char*>(aM+iA), sizeof(double));
            }
        }
    }
    if ((*g).ctr[10]>0)
    {
        for (int iC=0; iC<nCells; iC++)
        {
            e[iC].addIniPert(caseName,rank);
        }
    }

}  
void readPhysics(std::ifstream &fin, global* g)
{
    fin >> (*g).phy[0]; skipLine(fin, 1); // reference Mach Number
    fin >> (*g).phy[1]; skipLine(fin, 1); // reference Reynolds Number
    fin >> (*g).phy[2]; skipLine(fin, 1); // reference Froude Number
    fin >> (*g).phy[3]; skipLine(fin, 1); // reference Prandtl Number
    fin >> (*g).phy[4]; skipLine(fin, 1); // specific heats ratio
    fin >> (*g).phy[5];                   // Sutherland temperature ratio
    (*g).phy[6]=(*g).phy[4]*sq((*g).phy[0]); // gamma*Ma^2
    (*g).phy[7]=(*g).phy[4]/((*g).phy[4]-1.); // gamma/(gamma-1)
}
void readRelative(std::ifstream &fin, global* g)
{
    fin >> (*g).ctr[8]; skipLine(fin, 1);  // rotation axis index
    fin >> (*g).phy[8]; skipLine(fin, 1); //  revolution Strouhal number
    fin >> (*g).phy[9]; skipLine(fin, 1); //  amplitude of the rotational oscillation [Deg]
    fin >> (*g).phy[10]; skipLine(fin, 1); // Strouhal of the rotational oscillation
    fin >> (*g).phy[11]; skipLine(fin, 1); // initial phase[Deg] of the rotational oscillation (sinus)
    fin >> (*g).ctr[9]; skipLine(fin, 1);  // translation axis index
    fin >> (*g).phy[12]; skipLine(fin, 1); // amplitude of the translational oscillation
    fin >> (*g).phy[13]; skipLine(fin, 1); // Strouhal of the translational oscillation
    fin >> (*g).phy[14];                   // initial phase [Deg] of the translational oscillation (sinus)
}
void readRun(std::ifstream &fin, global* g)
{
    fin >> (*g).sch[0]; skipLine(fin, 1);  // polynomials max order
    fin >> (*g).sch[1]; skipLine(fin, 1);  // kind of LES model
    fin >> (*g).sch[2]; skipLine(fin, 1);  // kind of convective intercell numerical flux
    fin >> (*g).ctr[0]; skipLine(fin, 1);  // initial number of time steps
    fin >> (*g).ctr[1]; skipLine(fin, 1);  // final number of time steps
    fin >> (*g).ctr[2]; skipLine(fin, 1);  // solution damping period
    fin >> (*g).ctr[3]; skipLine(fin, 1);  // residual damping period
    fin >> (*g).ctr[4]; skipLine(fin, 1); // integrated quantities damping period
    fin >> (*g).ctr[5]; skipLine(fin, 1); // damping period of vtu/pvtu files
    fin >> (*g).snc[0]; skipLine(fin, 1); // time step
    fin >> (*g).sch[3]; skipLine(fin, 1); // source terms flag
    fin >> (*g).ctr[6]; skipLine(fin, 1); // number of computed forces/moments
    fin >> (*g).sch[4]; skipLine(fin, 1); // relative dynamics flag
    fin >> (*g).ctr[10];                  // initial perturbation flag
}
double readMesh(std::ifstream &fin1, std::ifstream &fin2, int nNodes, vector3D xN[], int nCells, physicalElement e[], global* g, computationalElement *cc, boundaryCondition BC[], int rank)
{
    double V=0.;
    int iN=-1, iC=-1, iVer[4];
    while (iN<nNodes-1) // reading of the grid point list (the first part of mesh file)
    {
        fin1 >> iN >> xN[iN][0] >> xN[iN][1] >> xN[iN][2]; //grid point reading from mesh file (each point is labeled by a sequential number and is defined by three coordinates)
    }
    intMatrix l(4,4);
    int i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15;
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
        e[iC].init(iC,cc,&xN[0],iVer,l,g,BC,rank); // physical cell initialization
        V+=e[iC].Jacobian();
    }
    return 4.*V/3.;  
}
void printOut(std::string caseName, physicalElement e[], int nCells, int it, global* g, int rank)
{
    std::string iterations=std::to_string(it); std::string process=std::to_string(rank);
    int Nm=nModes((*g).sch[0]);
    std::ofstream outputFileSol("./"+caseName+"/output/result_pr"+process+"_it"+iterations+".dat",std::ios::binary);
    int nData=nCells*Nm; double time=(*g).snc[0]*it;
    outputFileSol.write(reinterpret_cast<const char*>(&nData), sizeof(int));
    outputFileSol.write(reinterpret_cast<const char*>(&time), sizeof(double));
    double* aM; int N=(*g).sch[0];
    for (int i=0; i<nCells; i++)
    {
        outputFileSol.write(reinterpret_cast<const char*>(&N), sizeof(int));
        aM=e[i].toAM();
        for (int iA=0; iA<Nm*nEq; iA++)
        {
            outputFileSol.write(reinterpret_cast<const char*>(aM+iA), sizeof(double));
        }
    }
}
void makeVtu(std::string caseName, computationalElement* cc, physicalElement e[], int nCells, int it, global* g, int rank, int NRank)
{
    std::string iterations=std::to_string(it); std::string process=std::to_string(rank);
    iterations=std::to_string(it);
    std::fstream density;
    std::string fileNameD="./"+caseName+"/output/density_pr"+process+".dat";
    std::fstream energy;
    std::string fileNameE="./"+caseName+"/output/energy_pr"+process+".dat";
    std::fstream velocity;
    std::string fileNameV="./"+caseName+"/output/velocity_pr"+process+".dat";
    density.clear();
    density.open(fileNameD.c_str(),std::fstream::out);
    density << std::setprecision(12);
    energy.clear();
    energy.open(fileNameE.c_str(),std::fstream::out);
    energy << std::setprecision(12);
    velocity.clear();
    velocity.open(fileNameV.c_str(),std::fstream::out|std::fstream::out);
    velocity << std::setprecision(12);
    std::ofstream fileVisual;
    std::string fileName="./"+caseName+"/output/paraview_pr"+process+"_it"+iterations+".vtu";
    fileVisual.open(fileName.c_str(),std::fstream::out);
    fileVisual << std::setprecision(12);
    fileVisual << "<?xml version=\"1.0\"?>" << std::endl;
    fileVisual << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
    fileVisual << "  <UnstructuredGrid>" << std::endl;
    int Nm=nModes((*g).sch[0]);
    int Np=Nm, NpTot=0;
    fileVisual << "    <Piece NumberOfPoints=\"" << Np*nCells << "\" NumberOfCells=\"" << nCells << "\">" << std::endl;
    fileVisual << "      <Points>" << std::endl;
    fileVisual << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
    double rho, u, v, w, E;
    matrix AM(Np,nEq), qq(Np,nEq), q(nEq);
    vector3D xP, xA;
    for (int iC=0; iC<nCells; iC++)
    {
        for (int i=0; i<AM.size(); i++) {AM.set(i,*(e[iC].toAM()+i));}
        qq=(*cc).getPHIN()*AM;
        for (int k=0; k<Np; k++)
        {
            q=qq.row(k);
            xP=e[iC].getX((*cc).getRqN(k));
            if ((*g).sch[4]>0)
            {
                xA=absoluteX(xP,g);
                fileVisual << "        " <<  xA[0] << " " << xA[1] << " " << xA[2] << std::endl;
                absoluteQ(&q,xP,g);
            }
            else
            {
                fileVisual << "        " <<  xP[0] << " " << xP[1] << " " << xP[2] << std::endl;
            }
            rho=q.get(0);
            density << rho << std::endl;
            energy << q.get(4) << std::endl;
            velocity << q.get(1)/rho << " " << q.get(2)/rho << " " << q.get(3)/rho << std::endl;
        }
        NpTot+=Np;                        
    }
    density.close();
    density.clear();
    energy.close();
    energy.clear();
    velocity.close();
    velocity.clear();
    fileVisual << "        </DataArray>" << std::endl;
    fileVisual << "      </Points>" << std::endl;
    fileVisual << "      <Cells>" << std::endl;
    fileVisual << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
    int p=0;
    for (int iC=0; iC<nCells; iC++)
    {    
        switch((*g).sch[0])
        {
            case 1:
                fileVisual << "        " << p << " " << p+1 << " " << p+2 << " " << p+3 << std::endl;
                p+=4;
            break;
            case 2:
                 fileVisual << "        " << p << " " << p+1 << " " << p+2 << " " << p+3 << " " << p+4 << 
                                      " " << p+5 << " "  << p+6 << " " << p+7 << " " << p+8 << " " << p+9 << std::endl;
                 p+=10;
            break;
            case 3:
                 fileVisual << "        " << p << " " << p+1 << " " << p+2 << " " << p+3 << " " << p+4 << 
                                      " " << p+5 << " "  << p+6 << " " << p+7 << " " << p+8 << " " << p+9 <<
                                      " " << p+10 << " "  << p+11 << " " << p+12 << " " << p+13 << " " << p+14 <<
                                      " " << p+15 << " "  << p+16 << " " << p+17 << " " << p+18 << " " << p+19 << std::endl;
                 p+=20;
            break;
        }
    }
    fileVisual << "        </DataArray>" << std::endl;
    fileVisual << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl;
    for (int iC=0; iC<nCells; iC++)
    {    
        fileVisual << "        " << (iC+1)*Np << std::endl;
    }
    fileVisual << "        </DataArray>" << std::endl;
    fileVisual << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
    for (int i=0; i<nCells; i++)
    {
        fileVisual << 71 << std::endl;
    }
    fileVisual << "        </DataArray>" << std::endl;
    fileVisual << "      </Cells>" << std::endl;
    fileVisual << "      <PointData Scalars=\"density\" Vectors=\"velocity\">" << std::endl;
//
    fileVisual << "        <DataArray type=\"Float64\" Name=\"density\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
    density.open(fileNameD.c_str(),std::fstream::in);
    for (int k=0; k<NpTot; k++)
    {
        density >> rho;
        fileVisual << "        " << rho << std::endl;
    }
    density.close();
    fileVisual << "        </DataArray>" << std::endl;
//
    fileVisual << "        <DataArray type=\"Float64\" Name=\"energy\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
    energy.open(fileNameE.c_str(),std::fstream::in);
    for (int k=0; k<NpTot; k++)
    {
        energy >> E;
        fileVisual << "        " << E << std::endl;
    }
    energy.close();
    fileVisual << "        </DataArray>" << std::endl;
//
    fileVisual << "        <DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
    velocity.open(fileNameV.c_str(),std::fstream::out|std::fstream::in);
    for (int k=0; k<NpTot; k++)
    {
        velocity >> u >> v >> w;
        fileVisual << "        " << u << " " << v << " " << w << std::endl;
    }
    velocity.close();
    fileVisual << "        </DataArray>" << std::endl;
    fileVisual << "      </PointData>" << std::endl;
    fileVisual << "      <CellData>" << std::endl;
    fileVisual << "      </CellData>" << std::endl;
    fileVisual << "    </Piece>" << std::endl;
    fileVisual << "  </UnstructuredGrid>" << std::endl;
    fileVisual << "</VTKFile>" << std::endl;
    fileVisual.close();
    if (rank==0)
    {
        std::string s="./"+caseName+"/output/paraview"+"_it"+iterations+".pvtu"; std::ofstream glbVisual(s); // global file for visualization
        glbVisual << "<?xml version=\"1.0\"?>" << std::endl;
        glbVisual << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
        glbVisual << "  <PUnstructuredGrid GhostLevel=\"0\">" << std::endl;
        glbVisual << "    <PPoints>" << std::endl;
        glbVisual << "      <PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>" << std::endl;
        glbVisual << "    </PPoints>" << std::endl;
        glbVisual << "    <PCells>" << std::endl;
        glbVisual << "      <PDataArray type=\"Int32\" Name=\"connectivity\"/>" << std::endl;
        glbVisual << "      <PDataArray type=\"Int32\" Name=\"offsets\"/>" << std::endl;
        glbVisual << "      <PDataArray type=\"UInt8\" Name=\"types\"/>" << std::endl;
        glbVisual << "    </PCells>" << std::endl;
        glbVisual << "    <PPointData>" << std::endl;
        glbVisual << "      <PDataArray type=\"Float64\" Name=\"density\"/>" << std::endl;
        glbVisual << "      <PDataArray type=\"Float64\" Name=\"energy\"/>" << std::endl;
        glbVisual << "      <PDataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\"/>" << std::endl;
        glbVisual << "    </PPointData>" << std::endl;
        for (int iPr=0; iPr<NRank; iPr++)
        {
            fileName="paraview_pr"+std::to_string(iPr)+"_it"+iterations+".vtu";
            glbVisual << "    <Piece Source=\""+fileName+"\"/>" << std::endl;
        }
        glbVisual << "  </PUnstructuredGrid>" << std::endl;
        glbVisual << "</VTKFile>" << std::endl;        
        glbVisual.close();
    }
    remove(fileNameD.c_str());    
    remove(fileNameE.c_str());    
    remove(fileNameV.c_str());
}

