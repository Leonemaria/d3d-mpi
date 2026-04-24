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
        case(21): // Weak-Riemann isothermal no-slip
            iFile >> qB[6];
            iFile >> qB[7];
            iFile >> qB[8];
            iFile >> qB[9];
        break;
        case(22): // Weak-Riemann adiabatic no-slip
            iFile >> qB[6];
            iFile >> qB[7];
            iFile >> qB[8];
        break;
    }
    skipLine(iFile, 1);
}
