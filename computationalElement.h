#ifndef COMPUTATIONALELEMENT_H
#define COMPUTATIONALELEMENT_H
// computationalElement.h
#include<cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "constants.h"
#include "matrix.h"
#include "utilities.h"
#include "polynomials.h"
#include "vectors.h"
#include "global.h"
//
class computationalElement
{
    public:
// constructors/destructor declarations
        computationalElement(const global& glb);
        ~computationalElement();
// public methods declarations
        vector3D getRq(int k); // return the computational coordinates vector [r,s] of the k-th quadrature point
        matrix getPHI();
        matrix getPHI2();
        matrix getE();
        matrix getF();
        matrix getF2();
        matrix wQuad();
        void step_I(matrix* A_x, matrix* A_y, matrix* A_z, matrix f, matrix fS[], vector3D r_x, vector3D r_y, vector3D r_z);
        void step_II(double dt, int m, matrix* KA, matrix* A, matrix* A_0, matrix f[], matrix* fS, matrix* B, vector3D r_x, vector3D r_y, vector3D r_z);
        matrix subStep(double dt, matrix f[], matrix* fS, matrix* B, vector3D r_x, vector3D r_y, vector3D r_z);
        int extIndex(int i, int j);
    private:
// private variable declarations
        bool src;
        int N, Nm; //element polynomial order, number of modes
        int Npq, Npq2; // number of volume quadrature points, number of face quadrature points
        intMatrix extP;
        vector3D* rq=nullptr; // array of quadrature point coordinates [r,s,t]
        matrix w2D; // column matrix of triangle quadrature point weights
        matrix w3D; // column matrix of volume quadrature point weights
        matrix PHI, PHI2, E, E2, F, F2;
        matrix D_r, D_s, D_t;
};
//
#endif
