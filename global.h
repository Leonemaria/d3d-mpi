#ifndef GLOBAL_H
#define GLOBAL_H
// global.h
//
struct global
{
    double dt; // time step
    double phy[16]; // array of phyisical data
        // phy[0]: reference Mach number Ma        
        // phy[1]: reference Reynolds number Re            
        // phy[2]: reference Froude number Fr            
        // phy[3]: reference Prandtl number Pr            
        // phy[4]: air specific heats ratio gamma            
        // phy[5]: Suherland temperature ratio         
        // phy[6]: gamma*Ma^2
        // phy[7]: gamma/(gamma-1)
        // phy[8]: rotational speed (mean value)
        // phy[9]: amplitude of the rotational oscillation
        // phy[10]:frequency of the rotational oscillation
        // phy[11]:initial phase (of a sinus) of the rotational oscillation
        // phy[12]:amplitude of the traslational oscillation
        // phy[14]:frequency of the traslational oscillation
        // phy[15]:initial phase (of a sinus) of the traslational oscillation
    int ctr[7]; // array of run controls
     // ctr[0]: initial number of time steps
     // ctr[1]: final number of time steps
     // ctr[2]: solution damping period
     // ctr[3]: residual damping period
     // ctr[4]: traced quantities damping period
     // ctr[5]: number of computed forces/moments
     // ctr[6]: number of boundary conditions
    int sch[5]; // array of numerical scheme definers
     // sch[0]: polynomials max order
     // sch[1]: kind of LES model
     // sch[2]: kind of convective intercell numerical flux
     // sch[3]: source terms flag
     // sch[4]: relative dynamics flag
};
#endif
