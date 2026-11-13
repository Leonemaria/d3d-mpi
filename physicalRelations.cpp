//physicalRelations.cpp
#include "physicalRelations.h"
double energy(double rho, double rhoU, double rhoV, double rhoW, double p, global* g)
{
    double gaM2=(*g).phy[4]*sq((*g).phy[0]);
    return p/((*g).phy[4]-1.)+0.5*gaM2*(sq(rhoU)+sq(rhoV)+sq(rhoW))/rho;
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
