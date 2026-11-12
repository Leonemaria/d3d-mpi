//computationalElement.cpp
#include "computationalElement.h"

// constructor of computationalElement
computationalElement::computationalElement()
{
}
computationalElement::computationalElement(const global& glb)
{
//
    // the computational tetrahedral element (according to Hesthaven&Waburton 2008, p.409)
    src=glb.sch[3];
    N=glb.sch[0]; Nm=nModes(N); // polynomial order and number of nodal points
    int NF=N/2; NF=N;
    int NmF=nModes(NF);
    vector3D vr[4]; // vertices of right-angled tetrahedron
    vr[0].set(-1, -1, -1);
    vr[1].set(1, -1,-1);
    vr[2].set(-1, 1, -1);
    vr[3].set(-1, -1, 1);
    int iInd[Nm];
    int jInd[Nm];
    int kInd[Nm];
    int m=0;
    for (int i=0; i<N+1; i++)
    {
        for (int j=0; j<N+1-i; j++)
        {
            for (int k=0; k<N+1-i-j; k++)
            {
                iInd[m]=i; jInd[m]=j; kInd[m]=k; // i, j and k indices of the m-th mode
                m++;
            }
        }
    }
    int om, oM, mM;
    int iIs, jIs, kIs;
    for (int n=0; n<Nm; n++)
    {
        oM=iInd[0]+jInd[0]+kInd[0]; mM=0;
        for (int mm=1; mm<(Nm-n); mm++)
        {
            om=iInd[mm]+jInd[mm]+kInd[mm];
            if (om>oM) {oM=om; mM=mm;}
        }
        iIs=iInd[Nm-n-1]; jIs=jInd[Nm-n-1]; kIs=kInd[Nm-n-1];
        iInd[Nm-n-1]=iInd[mM]; jInd[Nm-n-1]=jInd[mM]; kInd[Nm-n-1]=kInd[mM];
        iInd[mM]=iIs; jInd[mM]=jIs; kInd[mM]=kIs;
    }
    matrix gP=gaussPointsAndWeights3D(N*2); // reads the points and weights for Gauss volume quadrature of order oq
    Npq=gP.nR();
    rq=new vector3D[Npq];
    for (int i=0; i<Npq; i++)
    {
        rq[i][0]=gP.get(i,0); rq[i][1]=gP.get(i,1); rq[i][2]=gP.get(i,2);
    }
    w3D=gP.col(3);
    D_r.dim(Nm,Npq); D_s.dim(Nm,Npq); D_t.dim(Nm,Npq);
    E.dim(Nm,Npq); PHI.dim(Npq,Nm);
    for (int i=0; i<Npq; i++)
    {
        for (int j=0; j<Nm; j++)
        {
            PHI.set(i,j,psi3D(iInd[j],jInd[j],kInd[j],rq[i]));
            E.set(j,i,w3D.get(i)*psi3D(iInd[j],jInd[j],kInd[j],rq[i]));
            D_r.set(j,i,w3D.get(i)*dPsi3Ddr(iInd[j],jInd[j],kInd[j],rq[i]));
            D_s.set(j,i,w3D.get(i)*dPsi3Dds(iInd[j],jInd[j],kInd[j],rq[i]));
            D_t.set(j,i,w3D.get(i)*dPsi3Ddt(iInd[j],jInd[j],kInd[j],rq[i]));
        }
    }
//
    gP=gaussPointsAndWeights2D(N*2);  // reads the points and weights for Gauss face quadrature of order oq2
    Npq2=gP.nR(); vector3D rq2[Npq2*4];
    w2D.dim(Npq2*4);
    extP.dim(Npq2,3); double toll=0.00000001;
    double lambda1,lambda2,lambda3; // baricentric coordinates
    double lambda1J, lambda2J, lambda3J;
    for (int i=0; i<Npq2; i++)
    {
        lambda1=gP.get(i,0); lambda2=gP.get(i,1); lambda3=gP.get(i,2);
        rq2[i]=lambda2*vr[1]+lambda3*vr[0]+lambda1*vr[2];
        w2D.set(i,gP.get(i,3));
        rq2[i+Npq2]=lambda2*vr[0]+lambda3*vr[1]+lambda1*vr[3];
        w2D.set(i+Npq2,gP.get(i,3));
        rq2[i+2*Npq2]=lambda2*vr[1]+lambda3*vr[2]+lambda1*vr[3];
        w2D.set(i+2*Npq2,gP.get(i,3));
        rq2[i+3*Npq2]=lambda2*vr[2]+lambda3*vr[0]+lambda1*vr[3];
        w2D.set(i+3*Npq2,gP.get(i,3));
        for (int j=0; j<Npq2; j++)
        {
            lambda1J=gP.get(j,0); lambda2J=gP.get(j,1); lambda3J=gP.get(j,2);
            if ((std::abs(lambda1-lambda3J)+std::abs(lambda2-lambda2J))<toll) {extP.set(i,0,j);}
            if ((std::abs(lambda1-lambda1J)+std::abs(lambda2-lambda3J))<toll) {extP.set(i,1,j);}
            if ((std::abs(lambda1-lambda2J)+std::abs(lambda3-lambda3J))<toll) {extP.set(i,2,j);}
        }
    }
    for (int iS=0; iS<4; iS++)
    {
        PHIS[iS].dim(Npq2,Nm);
        ES[iS].dim(Nm,Npq2);
        for (int i=0; i<Npq2; i++)
        {
            for (int j=0; j<Nm; j++)
            {
               PHIS[iS].set(i,j,psi3D(iInd[j],jInd[j],kInd[j],rq2[iS*Npq2+i]));
               ES[iS].set(j,i,w2D.get(i)*psi3D(iInd[j],jInd[j],kInd[j],rq2[iS*Npq2+i]));
            }
        }
    }
//    
}
// destructor of computationalElement
computationalElement::~computationalElement()
{
    delete[] rq; rq=nullptr;
}
// methods
int computationalElement::extIndex(int i, int j)
{
    return extP.get(i,j); // returs the index (on the face) of the external corresponding point 
}
vector3D computationalElement::getRq(int k)
{
    return rq[k];
}
matrix computationalElement::getPHI()
{
    return PHI;
}
matrix computationalElement::getPHIS(int iS)
{
    return PHIS[iS];
}
matrix computationalElement::getE()
{
    return E;
}
matrix computationalElement::wQuad(int d)
{
    if (d==2) {return w2D.part(0,Npq2);}
    if (d==3) {return w3D;}
}
matrix computationalElement::RHS_V(matrix f[], matrix* B, vector3D r_x, vector3D r_y, vector3D r_z)
{  
// Right Hand Side computation (with volume integral contributions only)
    if(src>0) { return (r_x[0]*D_r+r_x[1]*D_s+r_x[2]*D_t)*f[0]+(r_y[0]*D_r+r_y[1]*D_s+r_y[2]*D_t)*f[1]+(r_z[0]*D_r+r_z[1]*D_s+r_z[2]*D_t)*f[2]+E*(*B);}
    else {return (r_x[0]*D_r+r_x[1]*D_s+r_x[2]*D_t)*f[0]+(r_y[0]*D_r+r_y[1]*D_s+r_y[2]*D_t)*f[1]+(r_z[0]*D_r+r_z[1]*D_s+r_z[2]*D_t)*f[2];}
}
matrix computationalElement::RHS_S(matrix fS[])
{
// Right Hand Side computation (surface contributions)
    return ES[0]*fS[0]+ES[1]*fS[1]+ES[2]*fS[2]+ES[3]*fS[3];
}
void computationalElement::gradSolution(matrix* A_x, matrix* A_y, matrix* A_z, matrix f, matrix fS[], vector3D r_x, vector3D r_y, vector3D r_z,vector3D n[])
{
    matrix EfS[4]; for (int iS=0; iS<4; iS++) {EfS[iS]=ES[iS]*fS[iS];}
    (*A_x)=-(r_x[0]*D_r+r_x[1]*D_s+r_x[2]*D_t)*f+EfS[0]*n[0][0]+EfS[1]*n[1][0]+EfS[2]*n[2][0]+EfS[3]*n[3][0];
    (*A_y)=-(r_y[0]*D_r+r_y[1]*D_s+r_y[2]*D_t)*f+EfS[0]*n[0][1]+EfS[1]*n[1][1]+EfS[2]*n[2][1]+EfS[3]*n[3][1];
    (*A_z)=-(r_z[0]*D_r+r_z[1]*D_s+r_z[2]*D_t)*f+EfS[0]*n[0][2]+EfS[1]*n[1][2]+EfS[2]*n[2][2]+EfS[3]*n[3][2];
}
void computationalElement::RK_V(double d, int m, matrix* KA, matrix* A, matrix* A_0, matrix f[], matrix* B, vector3D r_x, vector3D r_y, vector3D r_z)
{
// five step fourth order Runge Kutta scheme SSPRK(5,4) (Gottlieb&Ketcheson&Shu p.23)
// RK substep (with volume intergral contributions only)
    switch (m)
    {
        case 0:
            (*A_0)=(*A); (*A)+=0.391752226571890*d*RHS_V(f,B,r_x,r_y,r_z);
        break;
        case 1:
            (*A)=0.444370493651235*(*A_0)+0.555629506348765*(*A)+0.368410593050371*d*RHS_V(f,B,r_x,r_y,r_z);
        break;
        case 2:
            (*KA)=0.517231671970585*(*A); (*A)=0.620101851488403*(*A_0)+0.379898148511597*(*A)+0.251891774271694*d*RHS_V(f,B,r_x,r_y,r_z);
        break;
        case 3:
            (*A)=0.821920045606868*(*A)+0.544974750228521*d*RHS_V(f,B,r_x,r_y,r_z);
        break;
        case 4:
            (*A)=(*KA)+0.386708617503269*(*A)+0.226007483236906*d*RHS_V(f,B,r_x,r_y,r_z);
         break;
    }
}
void computationalElement::RK_S(double d, int m, matrix* KA, matrix* A, matrix* A_0, matrix fS[])
{
// five step fourth order Runge Kutta scheme (addition of surface integral)
// RK substep (surface contributions)
    switch (m)
    {
        case 0:
            (*A)-=0.391752226571890*d*RHS_S(fS);
        break;
        case 1:
            (*A)-=0.368410593050371*d*RHS_S(fS);
        break;
        case 2:
            (*A)-=0.251891774271694*d*RHS_S(fS);
        break;
        case 3:
            (*A)-=0.544974750228521*d*RHS_S(fS); (*KA)+=0.116872329662213*(*A); (*A)+=0.178079954393132*(*A_0);
        break;
        case 4:
            (*A)-=0.226007483236906*d*RHS_S(fS);
         break;
    }
}
int computationalElement::time(int i, int m, double d) // RK stage time for non autononous cases (Liu and Song 2025, eq.2.2 and the following line)
{
    double t; i--;
    switch(m)
    {
        case(1):
            t=d*(i+0.39175222657189);
        break;
        case(2):
            t=d*(i+0.58607968931154);
        break;
        case(3):
            t=d*(i+0.47454236312140);
        break;
        case(4):
            t=d*(i+0.93501063096765);
        break;
    }
    return t;        
}
