//cases.cpp
#include "cases.h"
void iniFunc(std::string caseName, int i, vector3D x, matrix* q, double gam, double Ma)
{
    double rho=1., u=0., v=0., w=0., E=1., p=1.;
    if (caseName=="test")
    {
        u=1.5-1.5*sq(x[1]-1.);
        u=0.5-0.5*sq(x[0]-.25);
        E=energy(rho,rho*u,rho*v,rho*w,p,gam,Ma);
    }
    if (caseName=="poiseuille") // Poiseuille plane channel flow with channel height h=1
    {
        u=2./3.;
        E=energy(rho,rho*u,rho*v,rho*w,p,gam,Ma);
    }
    if (caseName=="adiabaticPoiseuille") // Poiseuille plane channel flow with channel height h=1
    {
        u=2./3.;
        E=energy(rho,rho*u,rho*v,rho*w,p,gam,Ma);
    }
    if (caseName=="taylorgreen")
    {
        p=1.+gam*sq(Ma)*(cos(2.*x[1])+cos(2.*x[0]))*(cos(2.*x[2])+2.)/16.;
        rho=p;
        v=sin(x[1])*cos(x[0])*cos(x[2]);
        u=-cos(x[1])*sin(x[0])*cos(x[2]);
        E=energy(rho,rho*u,rho*v,rho*w,p,gam,Ma);
    }
    (*q).set(i,0,rho);
    (*q).set(i,1,rho*u);
    (*q).set(i,2,rho*v);
    (*q).set(i,3,rho*w);
    (*q).set(i,4,E);
}
void srcFunc(std::string caseName, matrix* B, matrix* q, double gam, double Ma, double Re)
{
    int Np=(*q).nR();
    if (caseName=="adiabaticPoiseuille") // Adiabatic Poiseuille plane channel flow with channel height h=1
    {
        double p_x_gM2=8./Re;
        for (int i=0; i<Np; i++)
        {
            (*B).set(i,0,0.);
            (*B).set(i,1,p_x_gM2);
            (*B).set(i,2,0.);
            (*B).set(i,3,0.);
//            (*B).set(i,4,(*q).get(i,1)*p_x_gM2*gam*sq(Ma)*(1.-gam/(gam-1)/(Re*0.71))/(*q).get(i,0));
            (*B).set(i,4,(*q).get(i,1)*p_x_gM2*gam*sq(Ma)/(*q).get(i,0));
        }
    }
    if (caseName=="poiseuille") // Isothermal Poiseuille plane channel flow with channel height h=1
    {
        double p_x_gM2=8./Re;
        for (int i=0; i<Np; i++)
        {
            (*B).set(i,0,0.);
            (*B).set(i,1,p_x_gM2);
            (*B).set(i,2,0.);
            (*B).set(i,3,0.);
            (*B).set(i,4,(*q).get(i,1)*p_x_gM2*gam*sq(Ma)/(*q).get(i,0));
        }
    }
}
matrix varHist(std::string caseName, computationalElement* cE, double d, double dF, matrix *q, matrix* qA, matrix* q_x, matrix* q_y, matrix* q_z, matrix* qF,
               double S, double Re, int LES)
{
    double rho, u, v, w, T, mu, mut=0., s;
    symTensor Sr, tau;
    int Np=(*q).nR();
    matrix var;
    if (caseName=="taylorgreen")
    {
       var.dim(Np,2);
       for (int i=0; i<Np; i++)
        {
            rho=(*q).get(i,0); u=(*qA).get(i,0); v=(*qA).get(i,1); w=(*qA).get(i,2), T=(*qA).get(i,3);
            Sr=strainRate((*q_x).row(i),(*q_y).row(i),(*q_z).row(i));
            mu=Sutherland(T,S,Re);
            switch(LES)
            {
                case 0:
                    tau=mu*Sr.noTrace();
                break;
                case 1:
                    mut=sq(0.1*d)*Sr.Frobenius()/std::sqrt(2.);
                    tau=(mu+mut)*Sr.noTrace();
                break;
                case 2:
                {
                    double CS, CI, CQ, CJ;
                    s=Sr.Frobenius()/std::sqrt(2.);
                    double rhoD2S=(*q).get(i,0)*sq(d)*s;
                    coeffGermano(d,dF,(*qF).row(i),&CS,&CI,&CQ,&CJ);
                    CS=std::min(0.2,CS); // clipping on CS
                    double mut=std::max(-mu,CS*rhoD2S); // clipping on mut
                    tau=(mu+mut)*Sr.noTrace();
                    CI=std::clamp(CI,0.,0.1);
                    double isoTau=-CI*rhoD2S*s;
                    tau.add(0,0,isoTau); tau.add(1,1,isoTau); tau.add(2,2,isoTau);
                }
                break;
            }
            var.set(i,0,0.5*rho*(sq(u)+sq(v)+sq(w))-tau.trace()/2.);
            var.set(i,1,tau.doubleContr(Sr));
        }
    }
    if (caseName=="adiabaticPoiseuille")
    {
       var.dim(Np,1);
       for (int i=0; i<Np; i++)
        {
            u=(*qA).get(i,0);
            var.set(i,0,u);
        }
    }
    if (caseName=="poiseuille")
    {
       var.dim(Np,1);
       for (int i=0; i<Np; i++)
        {
            u=(*qA).get(i,0);
            var.set(i,0,u);
        }
    }
    if (caseName=="test")
    {
       var.dim(Np,1);
       for (int i=0; i<Np; i++)
        {
            u=(*qA).get(i,0);
            var.set(i,0,u);
        }
    }
    return var;
}

