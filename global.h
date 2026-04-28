#ifndef GLOBAL_H
#define GLOBAL_H

// global.h
//
struct global
{
    double dt; // time step
    double phy[6]; // array of phyisical data
        // phy[0]: reference Mach number Ma        
        // phy[1]: reference Reynolds number Re            
        // phy[2]: reference Froude number Fr            
        // phy[3]: reference Prandtl number Pr            
        // phy[4]: air specific heats ratio gamma            
        // phy[5]: Suherland temperature ratio         
    int ctr[5]; // array of run controls
     // ctr[0]: initial number of time steps
     // ctr[1]: final number of time steps
     // ctr[2]: solution damping period
     // ctr[3]: residual damping period
     // ctr[4]: integrated quantities damping period
    int sch[4]; // array of numerical scheme definers
     // sch[0]: polynomials max order
     // sch[1]: kind of LES model
     // sch[2]: kind of convective intercell numerical flux
     // sch[3]: source terms flag
};
#endif
