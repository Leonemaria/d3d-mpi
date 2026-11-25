#ifndef GLOBAL_H
#define GLOBAL_H
// global.h
//
struct global
{
    double snc[8]={}; // array for timing and syncronized relative kinematic values 
        // timing:
        // snc[0]: time step
        // snc[1]: Runge-Kutta time
        // relative motion:
        // snc[2]: present value of angular displacement
        // snc[3]: present value of angular speed
        // snc[4]: present value of angular acceleration
        // snc[5]: present value of translational displacement
        // snc[6]: present value of translational velocity
        // snc[7]: present value of translational acceleration
    double phy[15]={}; // array of phyisical data
        // phy[0]: reference Mach number Ma        
        // phy[1]: reference Reynolds number Re            
        // phy[2]: reference Froude number Fr            
        // phy[3]: reference Prandtl number Pr            
        // phy[4]: air specific heats ratio gamma            
        // phy[5]: Suherland temperature ratio         
        // phy[6]: gamma*Ma^2
        // phy[7]: gamma/(gamma-1)
        // phy[8]: revolution strouhal number speed
        // phy[9]: amplitude of the rotational oscillation [Deg]
        // phy[10]:Strouhal of the rotational oscillation
        // phy[11]:initial phase[Deg] of the rotational oscillation (sinus)
        // phy[12]:amplitude of the traslational oscillation
        // phy[13]:Strouhal of the traslational oscillation
        // phy[14]:initial phase[Deg] of the traslational oscillation (sinus)
    int ctr[11]={}; // array of run controls
     // ctr[0]: initial number of time steps
     // ctr[1]: final number of time steps
     // ctr[2]: solution damping period
     // ctr[3]: residual damping period
     // ctr[4]: traced quantities damping period
     // ctr[5]: damping period of vtu/pvtu files (for Paraview)
     // ctr[6]: number of computed forces/moments
     // ctr[7]: number of boundary conditions
     
     // ctr[8]: index of the rotation axis
     // ctr[9]: index of the translation axis
     // ctr[10]:initial perturbation flag
    int sch[5]={}; // array of numerical scheme definers
     // sch[0]: polynomials max order
     // sch[1]: kind of LES model
     // sch[2]: kind of convective intercell numerical flux
     // sch[3]: source terms flag
     // sch[4]: relative dynamics flag
};
#endif
