#ifndef PHYSICALELEMENT_H
#define PHYSICALELEMENT_H
// physicalElement.h
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "constants.h"
#include "matrix.h"
#include "tensor.h"
#include "utilities.h"
#include "computationalElement.h"
#include "boundaryCondition.h"
#include "physicalRelations.h"
#include <chrono>
#include "global.h"
#include "vectors.h"
#include "cases.h"
#include "LES.h"

// this class defines the physical element (triangle)
//
class physicalElement
{
    public:
// constructors/destructor declarations
        physicalElement();
// public methods declarations
        void init(int my, computationalElement *c, vector3D x[], long iV[], intMatrix l, const global& g);
        double CFL(double dt);
        void convFlux(matrix cF[], double rho, double rhoU, double rhoV, double rhoW, double E);
        void convFluxes(matrix* qq);
        void convFluxes(int iS);
        double getAM(int i, int eq);
        matrix getHist();
        double getQAS(int i, int eq);
        matrix getQAS(int i);
        double getQP(int i, int eq);
        double getQS(int i, int eq);
        matrix getQS(int i);
        double getResidual();
        matrix getFlxS(int i);
        vector3D getX(vector3D rP);
        matrix HLL(int iS, matrix qInt, matrix qExt);
        double integral(matrix m);
        matrix integralM(matrix m);
        double Jacobian();
        matrix LaxFriedrichs(int iS, matrix qInt, matrix qExt, matrix fMean);
        int ord();
        void print();
        void setAM(int i, int eq, double a);
        void setBC(int iS);
        void setIniCond(std::string caseName, const physics& phy);
        void step_0(boundaryCondition BC[]);
        void step_I(std::string nameCase, physicalElement e[], boundaryCondition BC[], bool* dmpH);
        void step_II(double dt, int m, physicalElement e[], bool dmpR);
        void viscousFlux(matrix vF[], double u, double v, double w, symTensor tau, vector3D heat);
        void viscousFlux(matrix vF[], double u, double v, double w, symTensor tau, vector3D heat, vector3D taK);
        void viscFluxes(matrix* qq, matrix* qA, matrix* dx, matrix* dy, matrix* dz, matrix* qF);
        void viscFluxes(matrix* dx, matrix* dy, matrix* dz, matrix* qF, int iS);
        void boundaryFluxes(matrix* dx, matrix* dy, matrix* dz, int iS, boundaryCondition* BC);

    private:
// private variable declarations
        bool BS[4];
        int mySelf, N, Nm, NmF, Npq, Npq2;
        int iVer[4];
        double J, JS[4], d, dF, h, res;
        matrix rot[4], invRot[4];
        vector3D* x0;
        matrix flxq[3];
        matrix flxS;
        matrix H;
        intMatrix join;
        vector3D x_r, x_s, x_t, r_x, r_y, r_z;
        matrix qS, qAuxS;
        matrix A, KA, A_0, B;
        vector3D n[4];
        computationalElement* cE; // pointer to the computational element
        physics phy; scheme sch;
};
//
#endif
