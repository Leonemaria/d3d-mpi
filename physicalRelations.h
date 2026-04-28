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
double energy(double rho, double rhoU, double rhoV, double rhoW, double p, double gam, double Ma);
double pressure(double rho, double rhoU, double rhoV, double rhoW, double E, double gam, double Ma);
double pressure(matrix* q, double gam, double Ma);
double pressure(matrix q, double gam, double Ma);
double soundSpeed(double rho, double p, double Ma);
double soundSpeed(double rho, double rhoU, double rhoV, double rhoW, double E, double gam, double Ma);
double soundSpeed(matrix* q, double gam, double Ma);
double soundSpeed(matrix q, double gam, double Ma);
double Sutherland(double T, double S, double Re);
symTensor strainRate(matrix dx, matrix dy, matrix dz);
#endif

