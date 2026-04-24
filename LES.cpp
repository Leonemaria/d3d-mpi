//LES.cpp
#include "LES.h"
void coeffGermano(double d, double dF, matrix qF, double *CS, double *CI, double *CQ, double *CJ)
{
    double rhoF=qF.get(0), dsq=sq(d), dFsq=sq(dF);
    symTensor SF;
    SF.set(0,0,2.*(qF.get(39)+qF.get(54)-qF.get(1)*qF.get(48)/rhoF)/rhoF);
    SF.set(1,0,(qF.get(40)+qF.get(55)-qF.get(1)*qF.get(49)/rhoF+qF.get(42)+qF.get(57)-qF.get(2)*qF.get(48)/rhoF)/rhoF);
    SF.set(1,1,2.*(qF.get(43)+qF.get(58)-qF.get(2)*qF.get(49)/rhoF)/rhoF);
    SF.set(2,0,(qF.get(41)+qF.get(56)-qF.get(1)*qF.get(50)/rhoF+qF.get(45)+qF.get(60)-qF.get(3)*qF.get(48)/rhoF)/rhoF);
    SF.set(2,1,(qF.get(44)+qF.get(59)-qF.get(2)*qF.get(50)/rhoF+qF.get(46)+qF.get(61)-qF.get(3)*qF.get(49)/rhoF)/rhoF);
    SF.set(2,2,2.*(qF.get(47)+qF.get(62)-qF.get(3)*qF.get(50)/rhoF)/rhoF);
    double tr=SF.trace()/3.;
    double sFsq=SF.autoContr()/2.; double sF=std::sqrt(sFsq);
    symTensor L(qF.get(5)-sq(qF.get(1))/rhoF,qF.get(6)-(qF.get(1)*qF.get(2))/rhoF,qF.get(8)-sq(qF.get(2))/rhoF,
             qF.get(7)-qF.get(1)*qF.get(3)/rhoF,qF.get(9)-qF.get(2)*qF.get(3)/rhoF,qF.get(10)-sq(qF.get(3))/rhoF); // Leonard stress tensor
    symTensor R(dsq*qF.get(24)-dFsq*sF*(SF.get(0,0)-tr),dsq*qF.get(25)-dFsq*sF*SF.get(1,0),dsq*qF.get(26)-dFsq*sF*(SF.get(1,1)-tr),
    dsq*qF.get(27)-dFsq*sF*SF.get(2,0),dsq*qF.get(28)-dFsq*sF*SF.get(2,1),dsq*qF.get(29)-dFsq*sF*(SF.get(2,2)-tr));
    double denominator=R.autoContr();
    if (std::abs(denominator)<0.0000000001) {*CS=0.;}
    else{*CS=(L.noTrace()).doubleContr(R)/denominator;}
    denominator=dFsq*sFsq-dsq*qF.get(23);
    if (std::abs(denominator)<0.0000000001) {*CI=0.;}
    else{*CI=L.trace()/denominator;};
    vector3D LQ(qF.get(14)-qF.get(1)*qF.get(4)/rhoF,qF.get(15)-qF.get(2)*qF.get(4)/rhoF,qF.get(16)-qF.get(3)*qF.get(4)/rhoF);
    vector3D RQ;
    RQ[0]=dsq*qF.get(33)-dFsq*sF*(qF.get(30)+qF.get(51)-qF.get(4)*qF.get(48));
    RQ[1]=dsq*qF.get(34)-dFsq*sF*(qF.get(31)+qF.get(52)-qF.get(4)*qF.get(49));
    RQ[2]=dsq*qF.get(35)-dFsq*sF*(qF.get(32)+qF.get(53)-qF.get(4)*qF.get(50));
    denominator=RQ.dot(RQ);
    if (std::abs(denominator)<0.0000000001) {*CQ=0.;}
    else {*CQ=LQ.dot(RQ)/denominator;}
    double U2=(sq(qF.get(1))+sq(qF.get(2))+sq(qF.get(3)))/sq(rhoF);
    vector3D LJ(qF.get(11)-qF.get(1)*U2,qF.get(12)-qF.get(2)*U2,qF.get(13)-qF.get(3)*U2);
    vector3D RJ;
    RJ[0]=dsq*qF.get(36)-dFsq*sF*(qF.get(1)*(qF.get(39)+qF.get(54)-qF.get(1)*qF.get(48)/rhoF)+
                                  qF.get(2)*(qF.get(42)+qF.get(57)-qF.get(2)*qF.get(48)/rhoF)+
                                  qF.get(3)*(qF.get(45)+qF.get(60)-qF.get(3)*qF.get(48)/rhoF))/rhoF;
    RJ[1]=dsq*qF.get(37)-dFsq*sF*(qF.get(1)*(qF.get(40)+qF.get(55)-qF.get(1)*qF.get(49)/rhoF)+
                                  qF.get(2)*(qF.get(43)+qF.get(58)-qF.get(2)*qF.get(49)/rhoF)+
                                  qF.get(3)*(qF.get(46)+qF.get(61)-qF.get(3)*qF.get(49)/rhoF))/rhoF;
    RJ[2]=dsq*qF.get(38)-dFsq*sF*(qF.get(1)*(qF.get(41)+qF.get(56)-qF.get(1)*qF.get(50)/rhoF)+
                                  qF.get(2)*(qF.get(44)+qF.get(59)-qF.get(2)*qF.get(50)/rhoF)+
                                  qF.get(3)*(qF.get(47)+qF.get(62)-qF.get(3)*qF.get(50)/rhoF))/rhoF;
    denominator=RJ.dot(RJ);
    if (std::abs(denominator)<0.0000000001) {*CJ=0.;}
    else{*CJ=LJ.dot(RJ)/denominator;}
}
matrix qLES(matrix *q, matrix *qA, matrix *dx, matrix *dy, matrix *dz, int LES)
{
    int nr=(*q).nR(); matrix qL; symTensor t, td; double rho, u, v, w, T, g, s;

    switch(LES)
    {
        case(2):
            qL.dim(nr,63);
            for (int i=0; i<nr; i++)
            {
                rho=(*q).get(i,0);
                T=(*qA).get(i,3);
                t=strainRate((*dx).row(i),(*dy).row(i),(*dz).row(i));
                qL.set(i,0,rho);                          // rho     
                qL.set(i,1,(*q).get(i,1));                // rhou
                qL.set(i,2,(*q).get(i,2));                // rhov
                qL.set(i,3,(*q).get(i,3));                // rhow
                qL.set(i,4,rho*T);                        // rhoT
                qL.set(i,5,(*q).get(i,1)*(*qA).get(i,0)); // rhouu
                qL.set(i,6,(*q).get(i,1)*(*qA).get(i,1)); // rhouv=rhovu
                qL.set(i,7,(*q).get(i,1)*(*qA).get(i,2)); // rhouw=rhowu
                qL.set(i,8,(*q).get(i,2)*(*qA).get(i,1)); // rhovv
                qL.set(i,9,(*q).get(i,2)*(*qA).get(i,2)); // rhovw=rhowv
                qL.set(i,10,(*q).get(i,3)*(*qA).get(i,2));// rhoww
                g=sq((*qA).get(i,0))+sq((*qA).get(i,1))+sq((*qA).get(i,2));
                qL.set(i,11,(*q).get(i,1)*g);             // rhou(uu+vv+ww)
                qL.set(i,12,(*q).get(i,2)*g);             // rhov(uu+vv+ww)
                qL.set(i,13,(*q).get(i,3)*g);             // rhow(uu+vv+ww)
                qL.set(i,14,(*q).get(i,1)*T);             // rhouT
                qL.set(i,15,(*q).get(i,2)*T);             // rhovT
                qL.set(i,16,(*q).get(i,3)*T);             // rhowT
                qL.set(i,17,rho*t.get(0,0));              // rhoS00
                qL.set(i,18,rho*t.get(1,0));              // rhoS10=rhoSd10
                qL.set(i,19,rho*t.get(1,1));              // rhoS11
                qL.set(i,20,rho*t.get(2,0));              // rhoS20=rhoSd20
                qL.set(i,21,rho*t.get(2,1));              // rhoS21=rhoSd21
                qL.set(i,22,rho*t.get(2,2));              // rhoS22
                s=t.autoContr()/2.;
                qL.set(i,23,rho*s);                       // rho|S|^2
                s=std::sqrt(s);
                g=t.trace()/3.;
                qL.set(i,24,rho*s*(t.get(0,0)-g));        // rho|S|Sd00
                qL.set(i,25,rho*s*t.get(1,0));            // rho|S|Sd10
                qL.set(i,26,rho*s*(t.get(1,1)-g));        // rho|S|Sd11
                qL.set(i,27,rho*s*t.get(2,0));            // rho|S|Sd20
                qL.set(i,28,rho*s*t.get(2,1));            // rho|S|Sd21
                qL.set(i,29,rho*s*(t.get(2,2)-g));        // rho|S|Sd22
                qL.set(i,30,rho*(*dx).get(i,3));          // rhodT/dx
                qL.set(i,31,rho*(*dy).get(i,3));          // rhodT/dy
                qL.set(i,32,rho*(*dz).get(i,3));          // rhodT/dz
                qL.set(i,33,rho*s*(*dx).get(i,3));        // rho|S|dT/dx
                qL.set(i,34,rho*s*(*dy).get(i,3));        // rho|S|dT/dy
                qL.set(i,35,rho*s*(*dz).get(i,3));        // rho|S|dT/dz
                g=(*qA).get(i,0)*(*dx).get(i,0)+(*qA).get(i,1)*(*dx).get(i,1)+(*qA).get(i,2)*(*dx).get(i,2);
                qL.set(i,36,rho*s*g);                     // rho|S|d(uu+vv+ww)/dx/2          
                g=(*qA).get(i,0)*(*dy).get(i,0)+(*qA).get(i,1)*(*dy).get(i,1)+(*qA).get(i,2)*(*dy).get(i,2);
                qL.set(i,37,rho*s*g);                     // rho|S|d(uu+vv+ww)/dy/2          
                g=(*qA).get(i,0)*(*dz).get(i,0)+(*qA).get(i,1)*(*dz).get(i,1)+(*qA).get(i,2)*(*dz).get(i,2);
                qL.set(i,38,rho*s*g);                     // rho|S|d(uu+vv+ww)/dz/2                
                qL.set(i,39,rho*(*dx).get(i,0));          // rhodu/dx
                qL.set(i,40,rho*(*dy).get(i,0));          // rhodu/dy
                qL.set(i,41,rho*(*dz).get(i,0));          // rhodu/dz
                qL.set(i,42,rho*(*dx).get(i,1));          // rhodv/dx
                qL.set(i,43,rho*(*dy).get(i,1));          // rhodv/dy
                qL.set(i,44,rho*(*dz).get(i,1));          // rhodv/dz
                qL.set(i,45,rho*(*dx).get(i,2));          // rhodw/dx
                qL.set(i,46,rho*(*dy).get(i,2));          // rhodw/dy
                qL.set(i,47,rho*(*dz).get(i,2));          // rhodw/dz
                
                qL.set(i,48,(*dx).get(i,4));              // drho/dx
                qL.set(i,49,(*dy).get(i,4));              // drho/dy
                qL.set(i,50,(*dz).get(i,4));              // drho/dz
                
                u=(*qA).get(0); v=(*qA).get(1); w=(*qA).get(2);
                qL.set(i,51,T*(*dx).get(i,4));            // Tdrho/dx
                qL.set(i,52,T*(*dy).get(i,4));            // Tdrho/dy
                qL.set(i,53,T*(*dz).get(i,4));            // Tdrho/dz
                qL.set(i,54,u*(*dx).get(i,4));            // udrho/dx
                qL.set(i,55,u*(*dy).get(i,4));            // udrho/dy
                qL.set(i,56,u*(*dz).get(i,4));            // udrho/dz
                qL.set(i,57,v*(*dx).get(i,4));            // vdrho/dx
                qL.set(i,58,v*(*dy).get(i,4));            // vdrho/dy
                qL.set(i,59,v*(*dz).get(i,4));            // vdrho/dz
                qL.set(i,60,w*(*dx).get(i,4));            // wdrho/dx
                qL.set(i,61,w*(*dy).get(i,4));            // wdrho/dy
                qL.set(i,62,w*(*dz).get(i,4));            // wdrho/dz
                //
            }
        break;
    }
    return qL;
}
