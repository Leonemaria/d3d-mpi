#ifndef LES_H
#define LES_H
//utilities.h
#include <cmath>
#include "matrix.h"
#include "tensor.h"
#include "vectors.h"
#include "utilities.h"
#include "constants.h"
#include "computationalElement.h"
#include "physicalRelations.h"
// this is a collection of usefull physical relations
void coeffGermano(double d, double df, matrix qF, double *CS, double *CI, double *CQ, double *CJ);
matrix qLES(matrix *q, matrix *qA, matrix *dx, matrix *dy, matrix *dz, int LES);
#endif

