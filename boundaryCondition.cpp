//boundaryCondition.cpp
#include "boundaryCondition.h"

// constructor of boundaryCondition
boundaryCondition::boundaryCondition()
{
// it is initialized as a Neumann condition
    kind[0]=0; kind[1]=-1;
    for (int eq=0; eq<nEq; eq++)
    {
        qB[eq]=0.;
    }
}
// methods
int boundaryCondition::getKind(int i)
{
    return kind[i];
}
double boundaryCondition::getQ(int eq)
{
    return qB[eq];
}
void boundaryCondition::input(std::ifstream &iFile, double a)
{
    iFile >> kind[0]; skipLine(iFile, 1);
    switch(kind[0])
    {
        case(11): // Weak-Riemann farfield (inlet/outlet) condition
        {
            double alph=a*radPerDeg, u, v;
            iFile >> qB[0]; // rho_inf
            iFile >> qB[5]; // p_inf
            iFile >> u; // u_inf at Alpha=0
            iFile >> v; // v_inf at Alpha=0
            qB[6]=u*cos(alph)-v*sin(alph); // u_inf with Alpha correction
            qB[7]=u*sin(alph)+v*cos(alph); // v_inf with Alpha correction
            iFile >> qB[8]; // w_inf
        }
        break;
        case(21): // Weak-Riemann isothermal no-slip
            iFile >> qB[6]; // u_wall
            iFile >> qB[7]; // v_wall
            iFile >> qB[8]; // w_wall
            iFile >> qB[9]; // T_wall
            iFile >> kind[1]; // index of force
        break;
        case(22): // Weak-Riemann adiabatic no-slip
            iFile >> qB[6]; // u_wall
            iFile >> qB[7]; // v_wall
            iFile >> qB[8]; // w_wall
            iFile >> kind[1]; // index of force
        break;
    }
    skipLine(iFile, 1);
}
double* boundaryCondition::doubleData()
{
    return &qB[0];
}
int* boundaryCondition::intData()
{
    return &kind[0];
}



