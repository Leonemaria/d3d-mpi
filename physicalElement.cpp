//physicalElement.cpp
#include "physicalElement.h"

// constructor of physicalElement
physicalElement::physicalElement()
{
}
// methods
void physicalElement::init(int my, computationalElement *c, vector3D x[], long iV[], intMatrix l, const global& glb)
{
    phy=glb.phy; sch=glb.sch;
    mySelf=my; N=glb.ctr.N; Nm=nModes(N); //N is the polynomial order. Nm is the number of modes
    int NF=N/2; NmF=nModes(NF);
    cE=c; // cE is the computational element pointer
    Npq=((*cE).getPHI()).nR(); //Npq is the number of quadrature points for volume integral
    Npq2=((*cE).getPHI2()).nR()/4; //Npq2 is the number of quadrature points for side (face) integral
    x0=x;
    join=l;
    for (int i=0; i<4; i++) {iVer[i]=iV[i];}    
    for (int i=0; i<4; i++) {BS[i]=false;}    
    vector3D edge01=(*(x0+iV[1])-*(x0+iV[0]));
    vector3D edge02=(*(x0+iV[2])-*(x0+iV[0]));
    vector3D edge03=(*(x0+iV[3])-*(x0+iV[0]));
    vector3D edge12=(*(x0+iV[2])-*(x0+iV[1]));
    vector3D edge13=(*(x0+iV[3])-*(x0+iV[1]));
    x_r=edge01/2.;
    x_s=edge02/2.;
    x_t=edge03/2.;
    J=x_r[0]*(x_s[1]*x_t[2]-x_s[2]*x_t[1])-x_r[1]*(x_s[0]*x_t[2]-x_s[2]*x_t[0])+x_r[2]*(x_s[0]*x_t[1]-x_s[1]*x_t[0]);
    JS[0]=(edge01.cross(edge02)).norm()/(2.*J);
    n[0]=-(edge01.cross(edge02)).normalized();
    JS[1]=(edge01.cross(edge03)).norm()/(2.*J);
    n[1]=(edge01.cross(edge03)).normalized();
    JS[2]=(edge12.cross(edge13)).norm()/(2.*J);
    n[2]=(edge12.cross(edge13)).normalized();
    JS[3]=(edge02.cross(edge03)).norm()/(2.*J);
    n[3]=-(edge02.cross(edge03)).normalized();
    h=std::pow(J*4./3.,1./3.); // lenght for CFL computation (Bakosi 2024)
    d=std::pow(J*4./3./Nm,1./3.); // grid LES scale (Abba' et al. 2015)
    dF=std::pow(J*4./3./NmF,1./3.); // grid LES larger scale
    r_x[0]=(x_s[1]*x_t[2]-x_s[2]*x_t[1])/J;
    r_y[0]=-(x_s[0]*x_t[2]-x_s[2]*x_t[0])/J;
    r_z[0]=(x_s[0]*x_t[1]-x_s[1]*x_t[0])/J;
    r_x[1]=-(x_r[1]*x_t[2]-x_r[2]*x_t[1])/J;
    r_y[1]=(x_r[0]*x_t[2]-x_r[2]*x_t[0])/J;
    r_z[1]=-(x_r[0]*x_t[1]-x_r[1]*x_t[0])/J;
    r_x[2]=(x_r[1]*x_s[2]-x_r[2]*x_s[1])/J;
    r_y[2]=-(x_r[0]*x_s[2]-x_r[2]*x_s[0])/J;
    r_z[2]=(x_r[0]*x_s[1]-x_r[1]*x_s[0])/J;
    flxq[0].dim(Npq,nEq); flxq[1].dim(Npq,nEq); flxq[2].dim(Npq,nEq); // dimensioning of the numerical fluxes matrices (one for each dimension)
    flxS.dim(4*Npq2,nEq);  // dimensioning of the face physical fluxes matrices (convective and viscous)
    qS.dim(4*Npq2,nEq); // dimensioning of the matrices containing the conservative variables interpolated over the face quadrature points
    qAuxS.dim(4*Npq2,nAux); // dimensioning of the matrices containing the auxiliary variables interpolated over the face quadrature points
    //
    if (sch.CIF==1)
    {
        vector3D ey, ez, eY(0.,1.,0.), eZ(0.,0.,1.);
        for (int i=0; i<4; i++)
        {
            if (std::abs(n[i][1])>std::abs(n[i][2]))
            {
                ey=(eZ.cross(n[i])).normalized();
                ez=(n[i].cross(ey)).normalized();
            }
            else
            {
                ez=(n[i].cross(eY)).normalized();
                ey=(ez.cross(n[i])).normalized();
            }
            rot[i].dim(3,3); invRot[i].dim(3,3);
            rot[i].set(0,0,n[i][0]); rot[i].set(0,1,n[i][1]); rot[i].set(0,2,n[i][2]);
            rot[i].set(1,0,ey[0]); rot[i].set(1,1,ey[1]); rot[i].set(1,2,ey[2]);
            rot[i].set(2,0,ez[0]); rot[i].set(2,1,ez[1]); rot[i].set(2,2,ez[2]);
            invRot[i]=rot[i].inv();
            rot[i]=rot[i].T();
            invRot[i]=invRot[i].T();   
        }
    }
    if (sch.src) {B.dim(Npq,nEq);}
    
}
double physicalElement::CFL(double dt)
{
    matrix qq=(*cE).getPHI()*A;
    double U=0.;
    for (int i=0; i<Npq; i++)
    {
        U+=std::sqrt(sq(qq.get(i,1))+sq(qq.get(i,2))+sq(qq.get(i,3)))/qq.get(i,0)+soundSpeed(qq.row(i),phy);
    }
    U=U/Npq;
    return U*dt/h;
}
void physicalElement::convFlux(matrix cF[], double rho, double rhoU, double rhoV, double rhoW, double E)
{
    double u=rhoU/rho, v=rhoV/rho, w=rhoW/rho, p=pressure(rho,rhoU,rhoV,rhoW,E,phy);
    cF[0].set(0,rhoU); cF[1].set(0,rhoV); cF[2].set(0,rhoW);
    cF[0].set(1,rhoU*u+p/phy.gaM2); cF[1].set(1,rhoU*v); cF[2].set(1,rhoU*w);
    cF[0].set(2,rhoV*u); cF[1].set(2,rhoV*v+p/phy.gaM2); cF[2].set(2,rhoV*w);
    cF[0].set(3,rhoW*u); cF[1].set(3,rhoW*v); cF[2].set(3,rhoW*w+p/phy.gaM2);
    cF[0].set(4,u*(E+p)); cF[1].set(4,v*(E+p)); cF[2].set(4,w*(E+p));
}
void physicalElement::convFluxes(matrix* qq)
// compute the convective flux
{
    matrix cF[3]; cF[0].dim(1,nEq); cF[1].dim(1,nEq); cF[2].dim(1,nEq);
    for (int i=0; i<Npq; i++)
    {
        convFlux(cF,(*qq).get(i,0),(*qq).get(i,1),(*qq).get(i,2),(*qq).get(i,3),(*qq).get(i,4));
        flxq[0].add(i,0,cF[0]); flxq[1].add(i,0,cF[1]); flxq[2].add(i,0,cF[2]);        
    }
}
void physicalElement::convFluxes(int iS)
// compute the convective flux
{
    matrix cF[3]; cF[0].dim(1,nEq); cF[1].dim(1,nEq); cF[2].dim(1,nEq);
    for (int i=iS*Npq2; i<(iS+1)*Npq2; i++)
    {
        convFlux(cF,qS.get(i,0),qS.get(i,1),qS.get(i,2),qS.get(i,3),qS.get(i,4));
        flxS.add(i,0,cF[0]*n[iS][0]+cF[1]*n[iS][1]+cF[2]*n[iS][2]);
    }
}
double physicalElement::getAM(int i, int eq)
{
    return A.get(i,eq);
}
matrix physicalElement::getHist()
{
    return H;
}
double physicalElement::getQAS(int i, int eq)
{
    return qAuxS.get(i,eq);
}
matrix physicalElement::getQAS(int i)
{
    return qAuxS.row(i);
}
double physicalElement::getQS(int i, int eq)
{
    return qS.get(i,eq);
}
matrix physicalElement::getQS(int i)
{
    return qS.row(i);
}
double physicalElement::getResidual()
{
    return res;
}
matrix physicalElement::getFlxS(int i)
{
    return flxS.row(i);
}
vector3D physicalElement::getX(vector3D rP)
{
    return (-(rP[0]+rP[1]+rP[2]+1.)*(*(x0+iVer[0]))+(rP[0]+1.)*(*(x0+iVer[1]))+(rP[1]+1.)*(*(x0+iVer[2]))+(rP[2]+1.)*(*(x0+iVer[3])))/2.;
}
matrix physicalElement::HLL(int iS, matrix qInt, matrix qExt)
// HLL (Mengaldo et al. 2014): approximated Riemann solver for intercell/boundary convective fluxes
{
    matrix flxHLL(1,5);
    double pInt=pressure(&qInt,phy), pExt=pressure(&qExt,phy);
    double cInt=soundSpeed(&qInt,phy), cExt=soundSpeed(&qExt,phy);
    qExt.set(0,1,qExt.part(0,1,1,3)*rot[iS]); qInt.set(0,1,qInt.part(0,1,1,3)*rot[iS]);
    double rhoI=qInt.get(0), rhoUI=qInt.get(1), rhoE=qExt.get(0), rhoUE=qExt.get(1);
    double uI=rhoUI/rhoI, vI=qInt.get(2)/rhoI, wI=qInt.get(3)/rhoI, EI=qInt.get(4);
    double uE=rhoUE/rhoE, vE=qExt.get(2)/rhoE, wE=qExt.get(3)/rhoE, EE=qExt.get(4);
    double SInt=std::min(0.,uI-cInt), SExt=std::max(0.,uE+cExt); // wave propagation speeds (Toro 1997)
    double dS=SExt-SInt; double C0=SExt/dS, C1=-SInt/dS, C2=SInt*SExt/dS;
    flxHLL.set(0,C0*qInt.get(1)+C1*qExt.get(1)+C2*(rhoE-rhoI));
    flxHLL.set(1,C0*(rhoUI*uI+pInt/phy.gaM2)+C1*(rhoUE*uE+pExt/phy.gaM2)+C2*(rhoUE-rhoUI));
    flxHLL.set(2,C0*rhoUI*vI+C1*rhoUE*vE+C2*(rhoE*vE-rhoI*vI));
    flxHLL.set(3,C0*rhoUI*wI+C1*rhoUE*wE+C2*(rhoE*wE-rhoI*wI));
    flxHLL.set(4,C0*uI*(EI+pInt)+C1*uE*(EE+pExt)+C2*(EE-EI));
    flxHLL.set(0,1,flxHLL.part(0,1,1,3)*invRot[iS]);
    return flxHLL; // convective numerical fluxes times side Jacobian

}
double physicalElement::integral(matrix m)
{
    return J*(m.T()*(*cE).wQuad()).get(0);
}
matrix physicalElement::integralM(matrix m)
{
    matrix mm(1,m.nC());
    for (int i=0; i<m.nC(); i++)
    {
        mm.set(i,integral(m.col(i)));
    }
    return mm;
}
double physicalElement::Jacobian()
{
    return J;
}
matrix physicalElement::LaxFriedrichs(int iS, matrix qInt, matrix qExt, matrix fMean)
// Lax-Friedrichs numerical flux computation
// NB: If fMean includes viscous flux this numerical flux will include viscous contribution too, according to BR1 scheme (Bassi&Rebay 1997)
{
    double UnInt=(qInt.get(1)*n[iS][0]+qInt.get(2)*n[iS][1]+qInt.get(3)*n[iS][2])/qInt.get(0);
    double UnExt=(qExt.get(1)*n[iS][0]+qExt.get(2)*n[iS][1]+qExt.get(3)*n[iS][2])/qExt.get(0);
    double cInt=soundSpeed(&qInt,phy), cExt=soundSpeed(&qExt,phy);
    double C=std::max(abs(UnInt)+cInt, abs(UnExt)+cExt)/2.;
    return fMean+C*(qInt-qExt);
}
int physicalElement::ord()
{
    return N;
}
void physicalElement::setAM(int i, int eq, double a)
{
    A.set(i,eq,a);
}
void physicalElement::setBC(int iS)
{
    BS[iS]=true;
}
void physicalElement::setIniCond(std::string caseName, const physics& phy)
{
    vector3D xP;
    matrix qq(Npq,nEq);
    for (int iP=0; iP<Npq; iP++)
    {
        xP=getX((*cE).getRq(iP));
        iniFunc(caseName,iP,xP,&qq,phy);
    }
    A=(*cE).getE()*qq;
}
void physicalElement::step_0(boundaryCondition BC[])
// computes conservative and auxiliary (primitive) variables on side quadrature points
// the conservative variables are stored in the matrix qS while the auxiliary variables are stored in the matrix qAuxS
{
    double rho;
    qS=(*cE).getPHI2()*A; // projects nodal solution on quadrature points: qS is the matrix containing side values (in quadrature points)
    for (int i=0; i<Npq2*4; i++)
    {
        rho=qS.get(i,0);
        qAuxS.set(i,0,qS.get(i,1)/rho); qAuxS.set(i,1,qS.get(i,2)/rho); qAuxS.set(i,2,qS.get(i,3)/rho); // set the side values (in i-th quadrature point) of the three first auxiliary variables (u, v and w) in qAuxS matrix
        qAuxS.set(i,3,pressure(qS.row(i),phy)/rho); // set the side values (in i-th quadrature point) of the fourth auxiliary variables (T) in qAuxS matrix
    }
    for (int iS=0; iS<4; iS++)
    {
        if (BS[iS])
        {
            int kind=BC[join.get(iS,2)].getKind();
            switch (kind)
            {
                case 21: // weak-Riemann no-slip isothermal condition (Mengaldo et al. 2014)
                    for (int i=iS*Npq2; i<(iS+1)*Npq2; i++)
                    {
                        for (int eq=0; eq<nAux; eq++) {qAuxS.set(i,eq,BC[join.get(iS,2)].getQ(eq+nEq+1));}
                    }
                break;
                case 22: // weak-Riemann no-slip adiabatic condition (Mengaldo et al. 2014)
                         // for this case the temperature is not set to a specific value (so it is kept as the internal value)
                    for (int i=iS*Npq2; i<(iS+1)*Npq2; i++)
                    {
                        for (int eq=0; eq<nAux-1; eq++) {qAuxS.set(i,eq,BC[join.get(iS,2)].getQ(eq+nEq+1));}
                    }
                break;
            }
        }
    }
}
void physicalElement::step_I(std::string nameCase, physicalElement e[], boundaryCondition BC[], bool* dmpH)
// computes the auxiliary variable gradients and physical fluxes on internal and side quadrature points
// the internal fluxes (convective-viscous) are stored in the matrix array fluxq[3]: the i-th matrix contains the i-th component of fluxes
// the side fluxes are stored
{
    int nGrad=nAux;
    if (sch.LES>1) {nGrad++;}
    matrix qMed(1,nGrad);
    double rho, rhoU, rhoV, rhoW, E, p;
    matrix numFlx[3]; numFlx[0].dim(4*Npq2,nGrad); numFlx[1].dim(4*Npq2,nGrad); numFlx[2].dim(4*Npq2,nGrad);
    flxS.zero();
    for (int iS=0; iS<4; iS++)
    {
        if (BS[iS])
        {
            for (int i=iS*Npq2; i<(iS+1)*Npq2; i++)
            {
                qMed.set(0,0,qAuxS.row(i)); // No condition                       
                if (sch.LES>1) {qMed.set(0,nGrad-1,qS.get(i,0));}
                numFlx[0].set(i,0,JS[iS]*qMed*n[iS][0]);
                numFlx[1].set(i,0,JS[iS]*qMed*n[iS][1]);
                numFlx[2].set(i,0,JS[iS]*qMed*n[iS][2]);
            }
        }
        else
        {
            int eJ=join.get(iS,0); // element connected to iS-th side
            int iExt;
            for (int i=iS*Npq2; i<(iS+1)*Npq2; i++)
            {
                iExt=(*cE).extIndex(i%Npq2,join.get(iS,2))+join.get(iS,1)*Npq2; // i index of the point on the connected element side
                qMed.set(0,0,0.5*(qAuxS.row(i)+e[eJ].getQAS(iExt)));
                if (sch.LES>1) {qMed.set(0,nGrad-1,0.5*(qS.get(i,0)+e[eJ].getQS(iExt,0)));}
                numFlx[0].set(i,0,JS[iS]*qMed*n[iS][0]);
                numFlx[1].set(i,0,JS[iS]*qMed*n[iS][1]);
                numFlx[2].set(i,0,JS[iS]*qMed*n[iS][2]);
            }
        }
    }
    matrix qq=(*cE).getPHI()*A;
    matrix qqAux(Npq,nGrad);
    for (int i=0; i<Npq; i++)
    {
     // computes the auxiliary variables on volume quadrature points
        rho=qq.get(i,0); rhoU=qq.get(i,1); rhoV=qq.get(i,2); rhoW=qq.get(i,3);
        E=qq.get(i,4); p=pressure(rho,rhoU,rhoV,rhoW,E,phy);
        qqAux.set(i,0,rhoU/rho); qqAux.set(i,1,rhoV/rho); qqAux.set(i,2,rhoW/rho); qqAux.set(i,3,p/rho); // auxiliary variables
        if (sch.LES>1) {qqAux.set(i,4,rho);}
    }
    // computes the auxiliary variable gradients
    matrix A_x, A_y, A_z;
    (*cE).step_I(&A_x,&A_y,&A_z,qqAux,numFlx,r_x,r_y,r_z); //  computation of mode amplitudes of auxiliary variable gradients
    matrix qq_x=(*cE).getPHI()*A_x; matrix qq_y=(*cE).getPHI()*A_y; matrix qq_z=(*cE).getPHI()*A_z; // compute the aux. variable gradients on volume quadrature points
    matrix q_xS=(*cE).getPHI2()*A_x; matrix q_yS=(*cE).getPHI2()*A_y; matrix q_zS=(*cE).getPHI2()*A_z; // computes auxiliary variable gradients on side quadrature points
    //
    matrix qL, qF;
    if (sch.LES>1) {qL=qLES(&qq,&qqAux,&qq_x,&qq_y,&qq_z,sch.LES); qF=(*cE).getF()*qL;} // compute the LES quantities
    if(*dmpH) // compute the possible other(statistic) variables
    {
        matrix var=varHist(nameCase,cE,d,dF,&qq,&qqAux,&qq_x,&qq_y,&qq_z,&qF,phy,sch);
        H=integralM(var);
    }     
    // compute the viscous part (including sgs) of physical fluxes on quadrature points
    viscFluxes(&qq,&qqAux,&qq_x,&qq_y,&qq_z,&qF); // viscous/sgs flux on internal quadrature points
    if (sch.LES>1) {qF=(*cE).getF2()*qL;}
    for (int iS=0; iS<4; iS++)
    {
        if (BS[iS])
        {
            boundaryFluxes(&q_xS,&q_yS,&q_zS,iS,&BC[join.get(iS,2)]);
        }
        else
        {
            viscFluxes(&q_xS,&q_yS,&q_zS,&qF,iS); // viscous/sgs flux on side quadrature points
            if (sch.CIF==0) {convFluxes(iS);}  // add convective flux on side quadrature points
        }
    }
    //
    // compute the convective part of physical fluxes on quadrature points
    convFluxes(&qq); // add convective flux on volume quadrature points
    if (sch.src) {srcFunc(nameCase,&B,&qq,phy);}   
 
}
void physicalElement::step_II(double dt, int m, physicalElement e[], bool dmpR)
{
    matrix numFlx(4*Npq2,nEq);
    for (int iS=0; iS<4; iS++)
    {
        if (BS[iS])
        {
            for (int i=iS*Npq2; i<(iS+1)*Npq2; i++)
            {
                numFlx.set(i,0,JS[iS]*flxS.row(i)); // BC flux times face Jacobian
            }
        }
        else
        {
            physicalElement* eJ=&e[join.get(iS,0)]; // element connected to iS-th side
            int iExt;
            switch (sch.CIF)
            {
                case 0: // Lax-Friedrichs+BR1
                    for (int i=iS*Npq2; i<(iS+1)*Npq2; i++)
                    {
                        iExt=(*cE).extIndex(i%Npq2,join.get(iS,2))+join.get(iS,1)*Npq2; // index of the point on the connected element side
                        numFlx.set(i,0,JS[iS]*LaxFriedrichs(iS,qS.row(i),(*eJ).getQS(iExt),0.5*(flxS.row(i)-(*eJ).getFlxS(iExt)))); // numerical flux times face Jacobian
                    }
                break;
                case 1: // HLL+BR1
                    for (int i=iS*Npq2; i<(iS+1)*Npq2; i++)
                    {
                        iExt=(*cE).extIndex(i%Npq2,join.get(iS,2))+join.get(iS,1)*Npq2; // index of the point on the connected element side
                        numFlx.set(i,0,JS[iS]*(0.5*(flxS.row(i)-(*eJ).getFlxS(iExt))+HLL(iS,qS.row(i),(*eJ).getQS(iExt)))); //  numerical flux times face Jacobian
                    }
                break;
            }
        }
    }
    (*cE).step_II(dt,m,&KA,&A,&A_0,flxq,&numFlx,&B,r_x,r_y,r_z); //computation of mode amplitude of condervative variables
    if (dmpR) {res=integral((*cE).getPHI()*(A.col(0)-A_0.col(0)));}
}
void physicalElement::viscousFlux(matrix vF[], double u, double v, double w, symTensor tau, vector3D heat)
{
    vF[0].set(0,0.); vF[1].set(0,0.); vF[2].set(0,0.);
    vF[0].set(1,-tau.get(0,0)); vF[1].set(1,-tau.get(0,1)); vF[2].set(1,-tau.get(0,2));
    vF[0].set(2,-tau.get(1,0)); vF[1].set(2,-tau.get(1,1)); vF[2].set(2,-tau.get(1,2));
    vF[0].set(3,-tau.get(2,0)); vF[1].set(3,-tau.get(2,1)); vF[2].set(3,-tau.get(2,2));
    vF[0].set(4,-phy.gaM2*(u*tau.get(0,0)+v*tau.get(0,1)+w*tau.get(0,2))+heat[0]);
    vF[1].set(4,-phy.gaM2*(u*tau.get(1,0)+v*tau.get(1,1)+w*tau.get(1,2))+heat[1]);
    vF[2].set(4,-phy.gaM2*(u*tau.get(2,0)+v*tau.get(2,1)+w*tau.get(2,2))+heat[2]);
}
void physicalElement::viscousFlux(matrix vF[], double u, double v, double w, symTensor tau, vector3D heat, vector3D taK)
{
    vF[0].set(0,0.); vF[1].set(0,0.); vF[2].set(0,0.);
    vF[0].set(1,-tau.get(0,0)); vF[1].set(1,-tau.get(0,1)); vF[2].set(1,-tau.get(0,2));
    vF[0].set(2,-tau.get(1,0)); vF[1].set(2,-tau.get(1,1)); vF[2].set(2,-tau.get(1,2));
    vF[0].set(3,-tau.get(2,0)); vF[1].set(3,-tau.get(2,1)); vF[2].set(3,-tau.get(2,2));
    vF[0].set(4,-phy.gaM2*(taK[0]+u*tau.get(0,0)+v*tau.get(0,1)+w*tau.get(0,2))+heat[0]);
    vF[1].set(4,-phy.gaM2*(taK[1]+u*tau.get(1,0)+v*tau.get(1,1)+w*tau.get(1,2))+heat[1]);
    vF[2].set(4,-phy.gaM2*(taK[2]+u*tau.get(2,0)+v*tau.get(2,1)+w*tau.get(2,2))+heat[2]);
}
void physicalElement::viscFluxes(matrix* qq, matrix* qA, matrix* dx, matrix* dy, matrix* dz, matrix* qF)
/// computation of viscous fluxes in the three internal flx matrices.
// (flx[0] for the x-fluxes, flx[1] for the y-fluxes and flx[2] for the z-fluxes)
{
    switch (sch.LES)
    {
        case 0: // ILES/DNS
        {
            double u, v, w, T, mu, k, ggm1=phy.gam/(phy.gam-1.);
            symTensor S, tau; vector3D gradT, heat;
            matrix vF[3]; vF[0].dim(1,nEq); vF[1].dim(1,nEq); vF[2].dim(1,nEq);
            for (int i=0; i<Npq; i++)
            {
                u=(*qA).get(i,0); v=(*qA).get(i,1); w=(*qA).get(i,2); T=(*qA).get(i,3);
                S=strainRate((*dx).row(i),(*dy).row(i),(*dz).row(i));
                mu=Sutherland(T,phy); k=mu/phy.Pr;
                gradT.set((*dx).get(i,3),(*dy).get(i,3),(*dz).get(i,3));
                tau=mu*S.noTrace();
                heat=-k*ggm1*gradT;
                viscousFlux(vF,u,v,w,tau,heat);
                flxq[0].set(i,0,vF[0]); flxq[1].set(i,0,vF[1]); flxq[2].set(i,0,vF[2]);        
            }
        }
        break;
        case 1: // Smagorinsky model
        {
            double u, v, w, T, mu, k, ggm1=phy.gam/(phy.gam-1.), mut, kt;
            symTensor S, tau; vector3D gradT, heat;
            matrix vF[3]; vF[0].dim(1,nEq); vF[1].dim(1,nEq); vF[2].dim(1,nEq);
            for (int i=0; i<Npq; i++)
            {
                u=(*qA).get(i,0); v=(*qA).get(i,1); w=(*qA).get(i,2); T=(*qA).get(i,3);
                S=strainRate((*dx).row(i),(*dy).row(i),(*dz).row(i));
                mu=Sutherland(T,phy); k=mu/phy.Pr;
                gradT.set((*dx).get(i,3),(*dy).get(i,3),(*dz).get(i,3));
                mut=sq(0.1*d)*S.Frobenius()/std::sqrt(2.);
                tau=(mu+mut)*S.noTrace();
                kt=mut/0.85;
                heat=-(k+kt)*ggm1*gradT;
                viscousFlux(vF,u,v,w,tau,heat);
                flxq[0].set(i,0,vF[0]); flxq[1].set(i,0,vF[1]); flxq[2].set(i,0,vF[2]);        
            }
        }
        break;
        case 2: // dynamic (isotropic) Germano model
        {
            double u, v, w, T, mu, k, ggm1=phy.gam/(phy.gam-1.), mut, kt;
            double CS, CI, CQ, CJ, s, rhoD2S, isoTau, dk;
            symTensor S, tau; vector3D gradT, heat;
            matrix vF[3]; vF[0].dim(1,nEq); vF[1].dim(1,nEq); vF[2].dim(1,nEq);
            vector3D taK;
            for (int i=0; i<Npq; i++)
            {
                u=(*qA).get(i,0); v=(*qA).get(i,1); w=(*qA).get(i,2); T=(*qA).get(i,3);
                S=strainRate((*dx).row(i),(*dy).row(i),(*dz).row(i));
                mu=Sutherland(T,phy); k=mu/phy.Pr;
                gradT.set((*dx).get(i,3),(*dy).get(i,3),(*dz).get(i,3));
                s=S.Frobenius()/std::sqrt(2.);
                rhoD2S=(*qq).get(i,0)*sq(d)*s;
                coeffGermano(d,dF,(*qF).row(i),&CS,&CI,&CQ,&CJ);
                CS=std::min(0.2,CS); // clipping on CS
                mut=std::max(-mu,CS*rhoD2S); // clipping on mut
                tau=(mu+mut)*S.noTrace(); // deviatory part of stress tensor
                CI=std::clamp(CI,0.,0.1); // clipping on CI
                isoTau=-CI*rhoD2S*s; // isotropic stress component
                tau.add(0,0,isoTau); tau.add(1,1,isoTau); tau.add(2,2,isoTau);
                CQ=std::min(0.2/0.85,CQ); // clipping on CQ
                kt=std::max(-mu/phy.Pr,CQ*rhoD2S); // clipping on kt
                heat=-(k+kt)*ggm1*gradT;
                dk=(*dx).get(0)*u+(*dx).get(1)*v+(*dx).get(2)*w;
                CJ=std::max(0.,CJ);
                taK[0]=CJ*rhoD2S*dk/2.;                        
                dk=(*dy).get(0)*u+(*dy).get(1)*v+(*dy).get(2)*w;
                taK[1]=CJ*rhoD2S*dk/2.;                         
                dk=(*dz).get(0)*u+(*dz).get(1)*v+(*dz).get(2)*w;
                taK[2]=CJ*rhoD2S*dk/2.;
                viscousFlux(vF,u,v,w,tau,heat,taK);
                flxq[0].set(i,0,vF[0]); flxq[1].set(i,0,vF[1]); flxq[2].set(i,0,vF[2]);        
            }
        }
        break;
    }
}
void physicalElement::viscFluxes(matrix* dx, matrix* dy, matrix* dz, matrix* qF, int iS)
// computation of face normal viscous flux
{
    switch(sch.LES)
    {
        case(0): // ILES/DNS
        {
            double u, v, w, T, mu, k, ggm1=phy.gam/(phy.gam-1.);
            symTensor S, tau; vector3D gradT, heat;
            matrix vF[3]; vF[0].dim(1,nEq); vF[1].dim(1,nEq); vF[2].dim(1,nEq);
            for (int i=iS*Npq2; i<(iS+1)*Npq2; i++)
            {
                u=qAuxS.get(i,0); v=qAuxS.get(i,1); w=qAuxS.get(i,2); T=qAuxS.get(i,3);
                S=strainRate((*dx).row(i),(*dy).row(i),(*dz).row(i));
                mu=Sutherland(T,phy); k=mu/phy.Pr;
                gradT.set((*dx).get(i,3),(*dy).get(i,3),(*dz).get(i,3));
                tau=mu*S.noTrace();
                heat=-k*ggm1*gradT;
                viscousFlux(vF,u,v,w,tau,heat);
                flxS.set(i,0,vF[0]*n[iS][0]+vF[1]*n[iS][1]+vF[2]*n[iS][2]);
            }
        }
        break;
        case 1: // Smagorinsky model
        {
            double u, v, w, T, mu, k, ggm1=phy.gam/(phy.gam-1.), mut, kt;
            symTensor S, tau; vector3D gradT, heat;
            matrix vF[3]; vF[0].dim(1,nEq); vF[1].dim(1,nEq); vF[2].dim(1,nEq);
            for (int i=iS*Npq2; i<(iS+1)*Npq2; i++)
            {
                u=qAuxS.get(i,0); v=qAuxS.get(i,1); w=qAuxS.get(i,2); T=qAuxS.get(i,3);
                S=strainRate((*dx).row(i),(*dy).row(i),(*dz).row(i));
                mu=Sutherland(T,phy); k=mu/phy.Pr;
                gradT.set((*dx).get(i,3),(*dy).get(i,3),(*dz).get(i,3));
                mut=sq(0.1*d)*S.Frobenius()/std::sqrt(2.);
                tau=(mu+mut)*S.noTrace();
                kt=mut/0.85;
                heat=-(k+kt)*ggm1*gradT;
                viscousFlux(vF,u,v,w,tau,heat);
                flxS.set(i,0,vF[0]*n[iS][0]+vF[1]*n[iS][1]+vF[2]*n[iS][2]);
            }
        }
        break;
        case 2: // dynamic (isotropic) Germano model
        {
            double u, v, w, T, mu, k, ggm1=phy.gam/(phy.gam-1.), mut, kt;
            double CS, CI, CQ, CJ, s, rhoD2S, isoTau, dk;
            symTensor S, tau; vector3D gradT, heat;
            matrix vF[3]; vF[0].dim(1,nEq); vF[1].dim(1,nEq); vF[2].dim(1,nEq);
            vector3D taK;
            for (int i=iS*Npq2; i<(iS+1)*Npq2; i++)
            {
                u=qAuxS.get(i,0); v=qAuxS.get(i,1); w=qAuxS.get(i,2); T=qAuxS.get(i,3);
                S=strainRate((*dx).row(i),(*dy).row(i),(*dz).row(i));
                mu=Sutherland(T,phy); k=mu/phy.Pr;
                gradT.set((*dx).get(i,3),(*dy).get(i,3),(*dz).get(i,3));
                s=S.Frobenius()/std::sqrt(2.);
                rhoD2S=qS.get(i,0)*sq(d)*s;
                coeffGermano(d,dF,(*qF).row(i),&CS,&CI,&CQ,&CJ);
                CS=std::min(0.2,CS); // clipping on CS
                mut=std::max(-mu,CS*rhoD2S); // clipping on mut
                tau=(mu+mut)*S.noTrace(); // deviatory part of stress tensor
                CI=std::clamp(CI,0.,0.1); // clipping on CI
                isoTau=-CI*rhoD2S*s; // isotropic stress component
                tau.add(0,0,isoTau); tau.add(1,1,isoTau); tau.add(2,2,isoTau);
                CQ=std::min(0.2/0.85,CQ); // clipping on CQ
                kt=std::max(-mu/phy.Pr,CQ*rhoD2S); // clipping on kt
                heat=-(k+kt)*ggm1*gradT;
                dk=(*dx).get(0)*u+(*dx).get(1)*v+(*dx).get(2)*w;
                CJ=std::max(0.,CJ);
                taK[0]=CJ*rhoD2S*dk/2.;                        
                dk=(*dy).get(0)*u+(*dy).get(1)*v+(*dy).get(2)*w;
                taK[1]=CJ*rhoD2S*dk/2.;                         
                dk=(*dz).get(0)*u+(*dz).get(1)*v+(*dz).get(2)*w;
                taK[2]=CJ*rhoD2S*dk/2.;
                viscousFlux(vF,u,v,w,tau,heat,taK);
                flxS.set(i,0,vF[0]*n[iS][0]+vF[1]*n[iS][1]+vF[2]*n[iS][2]);
            }
        }
        break;
    }
}
void physicalElement::boundaryFluxes(matrix* dx, matrix* dy, matrix* dz, int iS, boundaryCondition* BC)
// computation of a face normal viscous flux
{
    int kind=(*BC).getKind(); // the boundary type
    switch (kind)
    {
        case 0: // No condition
        {
            double mu, k, u, v, w, T, ggm1=phy.gam/(phy.gam-1.);
            symTensor S, tau; vector3D gradT, heat;
            matrix qInt(1,nEq); matrix flxB;
            matrix bF[3]; bF[0].dim(1,nEq); bF[1].dim(1,nEq); bF[2].dim(1,nEq);
            for (int i=iS*Npq2; i<(iS+1)*Npq2; i++)
            {
                qInt=qS.row(i);
                convFlux(bF,qInt.get(0),qInt.get(1),qInt.get(2),qInt.get(3),qInt.get(4));
                flxB=bF[0]*n[iS][0]+bF[1]*n[iS][1]+bF[2]*n[iS][2];
                u=qAuxS.get(i,0); v=qAuxS.get(i,1); w=qAuxS.get(i,2); T=qAuxS.get(i,3);
                S=strainRate((*dx).row(i),(*dy).row(i),(*dz).row(i));
                mu=Sutherland(T,phy); k=mu/phy.Pr;
                gradT.set((*dx).get(i,3),(*dy).get(i,3),(*dz).get(i,3));
                tau=mu*S.noTrace();
                heat=-k*ggm1*gradT;
                viscousFlux(bF,u,v,w,tau,heat);
                flxB.add(0,0,bF[0]*n[iS][0]+bF[1]*n[iS][1]+bF[2]*n[iS][2]);
                flxS.set(i,0,flxB);
            }
        }
        break;
        case 21: // weak-Riemann no-slip isothermal condition (Mengaldo et al. 2014)
        {
            double mu, k, u, v, w, T, ggm1=phy.gam/(phy.gam-1.);
            symTensor S, tau; vector3D gradT, heat;
            matrix qInt(1,nEq); matrix qExt(1,nEq); matrix flxB;
            matrix bF[3]; bF[0].dim(1,nEq); bF[1].dim(1,nEq); bF[2].dim(1,nEq);
            for (int i=iS*Npq2; i<(iS+1)*Npq2; i++)
            {
                qInt=qS.row(i);
                qExt.set(0,0,qInt.get(0)); // external density equal to internal value
                qExt.set(0,1,2.*(*BC).getQ(6)*qExt.get(0)-qInt.get(1)); // x momentum external value (extrapolating from inetrnal and boundary values)
                qExt.set(0,2,2.*(*BC).getQ(7)*qExt.get(0)-qInt.get(2)); // y momentum external value (extrapolating from inetrnal and boundary values)
                qExt.set(0,3,2.*(*BC).getQ(8)*qExt.get(0)-qInt.get(3)); // z momentum external value (extrapolating from inetrnal and boundary values)
                qExt.set(0,4,qInt.get(4)); // external energy equal to internl value
                convFlux(bF,qInt.get(0),qInt.get(1),qInt.get(2),qInt.get(3),qInt.get(4));
                flxB=(bF[0]*n[iS][0]+bF[1]*n[iS][1]+bF[2]*n[iS][2])/2.;
                convFlux(bF,qExt.get(0),qExt.get(1),qExt.get(2),qExt.get(3),qExt.get(4));
                flxB.add(0,0,(bF[0]*n[iS][0]+bF[1]*n[iS][1]+bF[2]*n[iS][2])/2.);
                u=qAuxS.get(i,0); v=qAuxS.get(i,1); w=qAuxS.get(i,2); T=qAuxS.get(i,3);
                S=strainRate((*dx).row(i),(*dy).row(i),(*dz).row(i));
                mu=Sutherland(T,phy); k=mu/phy.Pr;
                gradT.set((*dx).get(i,3),(*dy).get(i,3),(*dz).get(i,3));
                tau=mu*S.noTrace();
                heat=-k*ggm1*gradT;
                viscousFlux(bF,u,v,w,tau,heat);
                flxB.add(0,0,bF[0]*n[iS][0]+bF[1]*n[iS][1]+bF[2]*n[iS][2]);
                flxS.set(i,0,LaxFriedrichs(iS,qInt,qExt,flxB));
            }
        }
        break;
        case 22: // weak-Riemann no-slip adiabatic condition (Mengaldo et al. 2014)
        {
            double mu, k, u, v, w, T, ggm1=phy.gam/(phy.gam-1.);
            symTensor S, tau; vector3D gradT, heat;
            matrix qInt(1,nEq); matrix qExt(1,nEq); matrix flxB;
            matrix bF[3]; bF[0].dim(1,nEq); bF[1].dim(1,nEq); bF[2].dim(1,nEq);
            double q_n;
            for (int i=iS*Npq2; i<(iS+1)*Npq2; i++)
            {
                qInt=qS.row(i);
                qExt.set(0,0,qInt.get(0)); // external density equal to internal value
                qExt.set(0,1,2.*(*BC).getQ(6)*qExt.get(0)-qInt.get(1)); // x momentum external value (extrapolating from inetrnal and boundary values)
                qExt.set(0,2,2.*(*BC).getQ(7)*qExt.get(0)-qInt.get(2)); // y momentum external value (extrapolating from inetrnal and boundary values)
                qExt.set(0,3,2.*(*BC).getQ(8)*qExt.get(0)-qInt.get(3)); // z momentum external value (extrapolating from inetrnal and boundary values)
                qExt.set(0,4,qInt.get(4)); // external energy equal to internl value
                convFlux(bF,qInt.get(0),qInt.get(1),qInt.get(2),qInt.get(3),qInt.get(4));
                flxB=(bF[0]*n[iS][0]+bF[1]*n[iS][1]+bF[2]*n[iS][2])/2.;
                convFlux(bF,qExt.get(0),qExt.get(1),qExt.get(2),qExt.get(3),qExt.get(4));
                flxB.add(0,0,(bF[0]*n[iS][0]+bF[1]*n[iS][1]+bF[2]*n[iS][2])/2.);
                u=qAuxS.get(i,0); v=qAuxS.get(i,1); w=qAuxS.get(i,2); T=qAuxS.get(i,3);
                S=strainRate((*dx).row(i),(*dy).row(i),(*dz).row(i));
                mu=Sutherland(T,phy); k=mu/phy.Pr;
                gradT.set((*dx).get(i,3),(*dy).get(i,3),(*dz).get(i,3));
                q_n=gradT[0]*n[iS][0]+gradT[1]*n[iS][1]+gradT[2]*n[iS][2];
                gradT[0]-=q_n*n[iS][0]; gradT[1]-=q_n*n[iS][1]; gradT[2]-=q_n*n[iS][2]; // dT/dn=0 (thus zero heat at the wall)
                tau=mu*S.noTrace();
                heat=-k*ggm1*gradT;
                viscousFlux(bF,u,v,w,tau,heat);
                flxB.add(0,0,bF[0]*n[iS][0]+bF[1]*n[iS][1]+bF[2]*n[iS][2]);
                flxS.set(i,0,LaxFriedrichs(iS,qInt,qExt,flxB));
            }
        }
        break;
    }
}

