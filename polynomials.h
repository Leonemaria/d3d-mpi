#ifndef POLYNOMIALS_H
#define POLYNOMIALS_H
//polynomials.h
#include <iostream>
#include <fstream>
#include <cmath>
#include "matrix.h"
#include "vectors.h"
#include "utilities.h"
// this is a collection of functions related to polynomial basis
//typedef double (*iniPar)(vector2D, int);
double dJacobi(int n, int A, int B, double x);
matrix gaussPointsAndWeights2D(int p);
matrix gaussPointsAndWeights3D(int p);
double jacobi(int n, int A, int B, double x);
vector2D l_to_xy(double l1, double l2, double l3, vector2D v[]);
double psi3D(int i, int j, int k, vector3D rP);
double dPsi3Ddr(int i, int j, int k, vector3D rP);
double dPsi3Dds(int i, int j, int k, vector3D rP);
double dPsi3Ddt(int i, int j, int k, vector3D rP);
int nModes(int N);
double power(double b, int i);
#endif

