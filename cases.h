//cases.h
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include "constants.h"
#include "vectors.h"
#include "matrix.h"
#include "tensor.h"
#include "utilities.h"
#include "physicalRelations.h"
#include "global.h"
#include "LES.h"
//
// case distribution functions
void iniFunc(std::string caseName, int i, vector3D x, matrix* qP, const physics& phy);
void srcFunc(std::string caseName, matrix* B, matrix* q, const physics& phy);
matrix varHist(std::string caseName, computationalElement* cE,double d, double dF, matrix *q, matrix* qA, matrix* q_x, matrix*q_y, matrix* q_z, matrix* qF,
               const physics& phy, const scheme& sch);

