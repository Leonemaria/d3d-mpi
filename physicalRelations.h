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
double energy(double rho, double rhoU, double rhoV, double rhoW, double p, const physics& phy);
double pressure(double rho, double rhoU, double rhoV, double rhoW, double E, const physics& phy);
double pressure(matrix* q, const physics& phy);
double pressure(matrix q, const physics& phy);
double soundSpeed(double rho, double p, double Ma);
double soundSpeed(double rho, double rhoU, double rhoV, double rhoW, double E, const physics& phy);
double soundSpeed(matrix* q, const physics& phy);
double soundSpeed(matrix q, const physics& phy);
double Sutherland(double T, const physics& phy);
symTensor strainRate(matrix dx, matrix dy, matrix dz);
#endif

