//physicalRelations.cpp
#include "physicalRelations.h"
double energy(double rho, double rhoU, double rhoV, double rhoW, double p, double gam, double Ma)
{
    return p/(gam-1.)+0.5*gam*sq(Ma)*(sq(rhoU)+sq(rhoV)+sq(rhoW))/rho;
}
double pressure(double rho, double rhoU, double rhoV, double rhoW, double E, double gam, double Ma)
{
    return (gam-1.)*(E-0.5*gam*sq(Ma)*(sq(rhoU)+sq(rhoV)+sq(rhoW))/rho);
}
double pressure(matrix* q, double gam, double Ma)
{
    return (gam-1.)*((*q).get(4)-0.5*gam*sq(Ma)*(sq((*q).get(1))+sq((*q).get(2))+sq((*q).get(3)))/(*q).get(0));
}
double pressure(matrix q, double gam, double Ma)
{
    return (gam-1.)*(q.get(4)-0.5*gam*sq(Ma)*(sq(q.get(1))+sq(q.get(2))+sq(q.get(3)))/q.get(0));
}
double soundSpeed(double rho, double p, double Ma)
{
    return std::sqrt(p/rho)/Ma;
}
double soundSpeed(double rho, double rhoU, double rhoV, double rhoW, double E, double gam, double Ma)
{
    double p=pressure(rho,rhoU,rhoV,rhoW,E,gam,Ma);
    return std::sqrt(p/rho)/Ma;
}
double soundSpeed(matrix* q, double gam, double Ma)
{
    double p=pressure(q,gam,Ma);
    return std::sqrt(p/(*q).get(0))/Ma;
}
double soundSpeed(matrix q, double gam, double Ma)
{
    double p=pressure(q,gam,Ma);
    return std::sqrt(p/q.get(0))/Ma;
}
double Sutherland(double T, double S, double Re)
{
    return std::pow(T,1.5)*(1.+S)/(T+S)/Re;
}
symTensor strainRate(matrix dx, matrix dy, matrix dz) // strain rate tensor (defined as in Abba' et al. 2015 !!! thus no 1/2)
{
    symTensor Sd;
    Sd.set(0,0,2.*dx.get(0));
    Sd.set(1,0,dx.get(1)+dy.get(0)); Sd.set(1,1,2.*dy.get(1));
    Sd.set(2,0,dx.get(2)+dz.get(0)); Sd.set(2,1,dy.get(2)+dz.get(1)); Sd.set(2,2,2.*dz.get(2));
    return Sd; 
}
