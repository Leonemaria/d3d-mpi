//physicalRelations.cpp
#include "physicalRelations.h"
double energy(double rho, double rhoU, double rhoV, double rhoW, double p, global* g)
{
    double gaM2=(*g).phy[4]*sq((*g).phy[0]);
    return p/((*g).phy[4]-1.)+0.5*gaM2*(sq(rhoU)+sq(rhoV)+sq(rhoW))/rho;
}
void now(int i, int m, global* g)
{
    double time, amp, angFr, phi;
    switch(m)
    {
        case(1):
            time=(*g).snc[0]*(i+0.39175222657189);
        break;
        case(2):
            time=(*g).snc[0]*(i+0.58607968931154);
        break;
        case(3):
            time=(*g).snc[0]*(i+0.47454236312140);
        break;
        case(4):
            time=(*g).snc[0]*(i+0.93501063096765);
        break;
    }
    if ((*g).sch[4]>0)
    {
        amp=radPerDeg*(*g).phy[9]; angFr=2*pi*(*g).phy[10]; phi=radPerDeg*(*g).phy[11];
        (*g).snc[2]=2.*pi*(*g).phy[8]*time+amp*std::sin(angFr*time+phi);
        (*g).snc[3]=2.*pi*(*g).phy[8]+angFr*amp*std::cos(angFr*time+phi);
        (*g).snc[4]=-sq(angFr)*amp*std::sin(angFr*time+phi);
        amp=radPerDeg*(*g).phy[12]; angFr=2*pi*(*g).phy[13]; phi=radPerDeg*(*g).phy[14];
        (*g).snc[5]=amp*std::sin(angFr*time+phi);
        (*g).snc[6]=angFr*amp*std::cos(angFr*time+phi);
        (*g).snc[7]=-sq(angFr)*amp*std::sin(angFr*time+phi);
    }
    (*g).snc[1]=time;
}
void relativeQ(matrix* q, vector3D x, global* g) // corrects the momentum setting relative velocities
{
    double r, rho=(*q).get(0);
    vector3D v(0.,0.,0.), omega(0.,0.,0.);
    if ((*g).ctr[8]>0)
    {
        omega[(*g).ctr[8]-1]=(*g).snc[3];
        v=-omega.cross(x);
    }
    if ((*g).ctr[9]>0)
    {
        v[(*g).ctr[9]-1]-=(*g).snc[6]*cos((*g).snc[2]);
    }
    (*q).add(1,v[0]*rho);
    (*q).add(2,v[1]*rho);
    (*q).add(3,v[2]*rho);
}
void absoluteQ(matrix* q, vector3D x, global* g) // corrects the momentum setting "absolute" velocities
{
    double r, rho=(*q).get(0);
    vector3D v(0.,0.,0.), omega(0.,0.,0.);
    if ((*g).ctr[8]>0)
    {
        omega[(*g).ctr[8]-1]=(*g).snc[3];
        v=omega.cross(x);
    }
    if ((*g).ctr[9]>0)
    {
        v[(*g).ctr[9]-1]+=(*g).snc[6]*cos((*g).snc[2]);
    }
    (*q).add(1,v[0]*rho);
    (*q).add(2,v[1]*rho);
    (*q).add(3,v[2]*rho);
}
vector3D absoluteX(vector3D x, global* g)
{
    vector3D X; double th=(*g).snc[2];
    switch((*g).ctr[8])
    {
        case(1):
            X[1]=x[1]*cos(th)-x[2]*sin(th);
            X[2]=x[1]*sin(th)+x[2]*cos(th);
        break;
        case(2):
            X[2]=x[2]*cos(th)-x[0]*sin(th);
            X[0]=x[2]*sin(th)+x[0]*cos(th);
        break;
        case(3):
            X[0]=x[0]*cos(th)-x[1]*sin(th);
            X[1]=x[0]*sin(th)+x[1]*cos(th);
        break;
    }
    if ((*g).ctr[9]>0) {X[(*g).ctr[9]-1]+=(*g).snc[5];}
    return X;
}
vector3D fictitious(vector3D x, double u, double v, double w, global* g)  // this returns the fictitious force/mass
{
    vector3D a(0.,0.,0.), va(u,v,w), vr(0.,0.,0.), omega(0.,0.,0.), dOmdt(0.,0.,0.);
    if ((*g).ctr[8]>0)
    {
        omega[(*g).ctr[8]-1]=(*g).snc[3];
        dOmdt[(*g).ctr[8]-1]=(*g).snc[4];
        vr=va-omega.cross(x);
        a-=dOmdt.cross(x)+omega.cross(omega.cross(x))+2.*omega.cross(vr);
    }
    if ((*g).ctr[9]>0)
    {
        a[(*g).ctr[9]-1]-=(*g).snc[7];
    }
    return a;
}
double pressure(double rho, double rhoU, double rhoV, double rhoW, double E, global* g)
{
    double gaM2=(*g).phy[4]*sq((*g).phy[0]);
    return ((*g).phy[4]-1.)*(E-0.5*gaM2*(sq(rhoU)+sq(rhoV)+sq(rhoW))/rho);
}
double pressure(matrix* q, global* g)
{
    double gaM2=(*g).phy[4]*sq((*g).phy[0]);
    return ((*g).phy[4]-1.)*((*q).get(4)-0.5*gaM2*(sq((*q).get(1))+sq((*q).get(2))+sq((*q).get(3)))/(*q).get(0));
}
double pressure(matrix q, global* g)
{
    double gaM2=(*g).phy[4]*sq((*g).phy[0]);
    return ((*g).phy[4]-1.)*(q.get(4)-0.5*gaM2*(sq(q.get(1))+sq(q.get(2))+sq(q.get(3)))/q.get(0));
}
matrix qAuxiliary(matrix q, global* g)
{
    matrix qAS(1,nAux);
    double rho=q.get(0), rhoU=q.get(1), rhoV=q.get(2), rhoW=q.get(3);
    qAS.set(0,rhoU/rho); qAS.set(1,rhoV/rho); qAS.set(2,rhoW/rho);
    qAS.set(3,pressure(rho,rhoU,rhoV,rhoW,q.get(4),g)/rho);
    return qAS;   
}
double soundSpeed(matrix* q, global* g)
{
    double p=pressure(q,g);
    return std::sqrt(p/(*q).get(0))/(*g).phy[0];
}
double soundSpeed(matrix q, global* g)
{
    double p=pressure(q,g);
    return std::sqrt(p/q.get(0))/(*g).phy[0];
}
double Sutherland(double T, global* g)
{
    return std::pow(T,1.5)*(1.+(*g).phy[5])/(T+(*g).phy[5])/(*g).phy[1];  // mu/Re=(t^(1.5))*(1+S)/(T+S)/Re
}
symTensor strainRate(matrix dx, matrix dy, matrix dz) // strain rate tensor (defined as in Abba' et al. 2015 !!! thus no 1/2)
{
    symTensor Sd;
    Sd.set(0,0,2.*dx.get(0));
    Sd.set(1,0,dx.get(1)+dy.get(0)); Sd.set(1,1,2.*dy.get(1));
    Sd.set(2,0,dx.get(2)+dz.get(0)); Sd.set(2,1,dy.get(2)+dz.get(1)); Sd.set(2,2,2.*dz.get(2));
    return Sd; 
}
double vortex(int k, double G, double rc, double r) // tangential velocity/r of a 2D vortex
{
    double v=0.;
    switch(k) // switch on the kind of vortex
    {
        case(0):
            v=G/(2.*pi*std::sqrt(std::pow(rc,4)+std::pow(r,4)));
        break;
    }
    return v;
}
