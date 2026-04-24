//inputOutput.h
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include "vectors.h"
#include "physicalElement.h"
#include "constants.h"
#include "utilities.h"
#include "global.h"
//
// this is a collection of input/outpu functions
void initialConditions(std::string caseName, long nCells,physicalElement e[], const global& glb);
void readPhysics(std::ifstream &fin, global& glb);   
void readRun(std::ifstream &fin, global& glb);
double readMesh(std::ifstream &fin1, std::ifstream &fin2, long nNodes, vector3D xN[], long nCells, physicalElement e[], const global& glb, computationalElement *cc);
void printOut(std::string caseName, physicalElement e[], long nCells, int it, const global& glb);

