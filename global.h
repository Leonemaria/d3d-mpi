#ifndef GLOBAL_H
#define GLOBAL_H

// global.h
//
struct physics
{
    double Ma;   // reference Mach number Ma
    double Re;   // reference Reynolds number Re
    double Fr;   // reference Froude number Fr
    double Pr;   // reference Prandtl number Pr
    double gam;  // air specific heats ratio gamma
    double gaM2; // gamma*Ma^2
    double S;    // Suherland temperature ratio
};
struct control
{
    int rank;     // MPI process rank
    int N;        // polynomials max order
    int startStep;// initial number of time steps
    int endStep;  // final number of time steps
    int nStSave;  // solution damping period
    int nStDRes;  // residual damping period
    int nStDHist; // integrated quantities damping period
    double dt;    // time step
};
struct scheme
{
    int LES;      // kind of LES model
    int CIF;      // kind of convective intercell numerical flux
    bool src;     // source terms flag
};
struct global
{
    double dt;
    physics phy;
    control ctr;
    scheme  sch;
};
#endif
