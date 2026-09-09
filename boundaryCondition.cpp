//boundaryCondition.cpp
#include "boundaryCondition.h"

// constructor of boundaryCondition
boundaryCondition::boundaryCondition()
{
// it is initialized as a Neumann condition
    kind=0;
    for (int eq=0; eq<nEq; eq++)
    {
        qB[eq]=0.;
    }
}
// methods
int boundaryCondition::getKind()
{
    return kind;
}
double boundaryCondition::getQ(int eq)
{
    return qB[eq];
}
void boundaryCondition::input(std::ifstream &iFile)
{
    iFile >> kind; skipLine(iFile, 1);
    switch(kind)
    {
        case(11): // Weak-Riemann farfield (inlet/outlet) condition
            iFile >> qB[0]; // rho_inf
            iFile >> qB[5]; // p_inf
            iFile >> qB[6]; // u_inf
            iFile >> qB[7]; // v_inf
            iFile >> qB[8]; // w_inf
        break;
        case(21): // Weak-Riemann isothermal no-slip
            iFile >> qB[6]; // u_wall
            iFile >> qB[7]; // v_wall
            iFile >> qB[8]; // w_wall
            iFile >> qB[9]; // T_wall
        break;
        case(22): // Weak-Riemann adiabatic no-slip
            iFile >> qB[6]; // u_wall
            iFile >> qB[7]; // v_wall
            iFile >> qB[8]; // w_wall
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
    return &kind;
}



