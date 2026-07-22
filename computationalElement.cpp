//computationalElement.cpp
#include "computationalElement.h"

// constructor of computationalElement
computationalElement::computationalElement()
{
}
computationalElement::computationalElement(const global& glb)
{
//
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
    matrix EF(NmF,Npq); matrix PHIF(Npq,NmF);
    for (int i=0; i<Npq; i++)
    {
        for (int j=0; j<NmF; j++)
        {
            PHIF.set(i,j,psi3D(iInd[j],jInd[j],kInd[j],rq[i]));
            EF.set(j,i,w3D.get(i)*psi3D(iInd[j],jInd[j],kInd[j],rq[i]));
        }
    }
    F=PHIF*EF; // filtering matrix for points inside the volume
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
    PHI2.dim(4*Npq2,Nm);
    E2.dim(Nm,4*Npq2);
    for (int i=0; i<Npq2*4; i++)
    {
        for (int j=0; j<Nm; j++)
        {
           PHI2.set(i,j,psi3D(iInd[j],jInd[j],kInd[j],rq2[i]));
           E2.set(j,i,w2D.get(i)*psi3D(iInd[j],jInd[j],kInd[j],rq2[i]));
        }
    }
    matrix PHIF2(4*Npq2,NmF);
    for (int i=0; i<Npq2*4; i++)
    {
        for (int j=0; j<NmF; j++)
        {
            PHIF2.set(i,j,psi3D(iInd[j],jInd[j],kInd[j],rq2[i]));
        }
    }
    F2=PHIF2*EF; //filtering matrix for points on the sides
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
matrix computationalElement::getPHI2()
{
    return PHI2;
}
matrix computationalElement::getE()
{
    return E;
}
matrix computationalElement::getF()
{
    return F;
}
matrix computationalElement::getF2()
{
    return F2;
}
matrix computationalElement::wQuad()
{
    return w3D;
}
matrix computationalElement::subStep_a(double d, matrix f[], matrix* B, vector3D r_x, vector3D r_y, vector3D r_z)
{  
    if(src>0) { return d*((r_x[0]*D_r+r_x[1]*D_s+r_x[2]*D_t)*f[0]+(r_y[0]*D_r+r_y[1]*D_s+r_y[2]*D_t)*f[1]+(r_z[0]*D_r+r_z[1]*D_s+r_z[2]*D_t)*f[2]+E*(*B));}
    else {return d*((r_x[0]*D_r+r_x[1]*D_s+r_x[2]*D_t)*f[0]+(r_y[0]*D_r+r_y[1]*D_s+r_y[2]*D_t)*f[1]+(r_z[0]*D_r+r_z[1]*D_s+r_z[2]*D_t)*f[2]);}
}
matrix computationalElement::subStep_b(double d, matrix* fS)
{  
    return -d*E2*(*fS);
}
void computationalElement::step_I(matrix* A_x, matrix* A_y, matrix* A_z, matrix f, matrix fS[], vector3D r_x, vector3D r_y, vector3D r_z)
{
    (*A_x)=-(r_x[0]*D_r+r_x[1]*D_s+r_x[2]*D_t)*f+E2*fS[0];
    (*A_y)=-(r_y[0]*D_r+r_y[1]*D_s+r_y[2]*D_t)*f+E2*fS[1];
    (*A_z)=-(r_z[0]*D_r+r_z[1]*D_s+r_z[2]*D_t)*f+E2*fS[2];
}
void computationalElement::step_IIa(double d, int m, matrix* KA, matrix* A, matrix* A_0, matrix f[], matrix* B, vector3D r_x, vector3D r_y, vector3D r_z)
{
// five step fourth order Runge Kutta scheme SSPRK(5,4) (Gottlieb&Ketcheson&Shu pag.23)
    switch (m)
    {
        case 0:
            (*A_0)=(*A); (*A)+=subStep_a(0.391752226571890*d,f,B,r_x,r_y,r_z);
        break;
        case 1:
            (*A)=0.444370493651235*(*A_0)+0.555629506348765*(*A)+subStep_a(0.368410593050371*d,f,B,r_x,r_y,r_z);
        break;
        case 2:
            (*KA)=0.517231671970585*(*A); (*A)=0.620101851488403*(*A_0)+0.379898148511597*(*A)+subStep_a(0.251891774271694*d,f,B,r_x,r_y,r_z);
        break;
        case 3:
            (*A)=0.821920045606868*(*A)+subStep_a(0.544974750228521*d,f,B,r_x,r_y,r_z);
        break;
        case 4:
            (*A)=(*KA)+0.386708617503269*(*A)+subStep_a(0.226007483236906*d,f,B,r_x,r_y,r_z);
         break;
    }
}void computationalElement::step_IIb(double d, int m, matrix* KA, matrix* A, matrix* A_0, matrix* fS)
{
// five step fourth order Runge Kutta scheme (addition of surface integral)
    switch (m)
    {
        case 0:
            (*A)+=subStep_b(0.391752226571890*d,fS);
        break;
        case 1:
            (*A)+=subStep_b(0.368410593050371*d,fS);
        break;
        case 2:
            (*A)+=subStep_b(0.251891774271694*d,fS);
        break;
        case 3:
            (*A)+=subStep_b(0.544974750228521*d,fS); (*KA)+=0.116872329662213*(*A); (*A)+=0.178079954393132*(*A_0);
        break;
        case 4:
            (*A)+=subStep_b(0.226007483236906*d,fS);
         break;
    }
}
