#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H
// boundaryCondition.h
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "constants.h"
#include "utilities.h"
#include "physicalRelations.h"

// this class defines a boundary condition
//
class boundaryCondition
{
    public:
// constructors/destructor declarations
        boundaryCondition();
// public methods declarations
        int getKind();
        double getQ(int eq);
        void input(std::ifstream &iFile);
    private:
// private variable declarations
        int kind; // number of equation, kind of condition
        double qB[10]; // boundary values (depend on the kind of BC)
};
//
#endif
