#ifndef PHYSICALRELATIONS_H
#define PHYSICALRELATIONS_H
//utilities.h
#include <cmath>
#include "matrix.h"
#include "tensor.h"
#include "vectors.h"
#include "utilities.h"
#include "constants.h"
#include "computationalElement.h"
#include "global.h"
// this is a collection of usefull physical relations
double energy(double rho, double rhoU, double rhoV, double rhoW, double p, global* g);
double pressure(double rho, double rhoU, double rhoV, double rhoW, double E, global* g);
double pressure(matrix* q, global* g);
double pressure(matrix q, global* g);
matrix qAuxiliary(matrix q, global* g);
double soundSpeed(matrix* q, global* g);
double soundSpeed(matrix q, global* g);
double Sutherland(double T, global* g);
symTensor strainRate(matrix dx, matrix dy, matrix dz);
#endif

