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
        computationalElement();
        computationalElement(global* g);
        ~computationalElement();
// public methods declarations
        vector3D getRq(int i); // return the computational coordinates vector [r,s] of the i-th quadrature point in the volume
        vector3D getRqS(int iS, int i); // return the computational coordinates vector [r,s] of the i-th quadrature point of the k-th face
        vector3D getRqN(int i); // return the computational coordinates vector [r,s] of the i-th nodal (for graphics) point of the cell
        matrix getPHI();
        matrix getPHIS(int iS);
        matrix getPHIN();
        matrix getE();
        matrix wQuad(int d);
        void gradSolution(matrix* A_x, matrix* A_y, matrix* A_z, matrix f, matrix fS[], vector3D r_x, vector3D r_y, vector3D r_z, vector3D n[]);
        void RK_V(double dt, int m, matrix* KA, matrix* A, matrix* A_0, matrix f[], matrix* B, vector3D r_x, vector3D r_y, vector3D r_z);
        void RK_S(double dt, int m, matrix* KA, matrix* A, matrix* A_0, matrix fS[]);
        matrix RHS_V(matrix f[], matrix* B, vector3D r_x, vector3D r_y, vector3D r_z);
        matrix RHS_S(matrix fS[]);
        int extIndex(int i, int j);
    private:
// private variable declarations
        int src=0;
        int N, Nm; //element polynomial order, number of modes
        int Npq, Npq2; // number of volume quadrature points, number of face quadrature points
        intMatrix extP;
        vector3D* rq=nullptr; // array of volume quadrature point coordinates [r,s,t]
        vector3D* rq2=nullptr; // array of face quadrature point coordinates [r,s,t]
        vector3D* rqN=nullptr; // array of graphic nodal point coordinates [r,s,t]
        matrix w2D; // column matrix of triangle quadrature point weights
        matrix w3D; // column matrix of volume quadrature point weights
        matrix PHI, E;
        matrix PHIS[4], ES[4];
        matrix PHIN;
        matrix D_r, D_s, D_t;
};
//
#endif
