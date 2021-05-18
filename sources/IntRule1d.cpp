/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <iostream> 
#include "IntRule1d.h"
#include <vector>
#include <math.h>
#include <cmath>
#include "tpanic.h"
using namespace std;

#define PI 3.141592654

IntRule1d::IntRule1d(){}

IntRule1d::IntRule1d(int order) : IntRule(order) {
    SetOrder(order);
}

void IntRule1d::SetOrder(int order) {
    if(order<0||order>=6) DebugStop();

    int npoints = (order+1)/2 + (order+1)%2;
    fOrder = order;
    if(npoints==this->NPoints())return;

    fWeights.resize(npoints);
    fPoints.resize(npoints,1);
    // taken from https://pomax.github.io/bezierinfo/legendre-gauss.html
    switch(npoints){
        case 1:
            fWeights[0] = 2;        fPoints(0, 0) = 0.0;                 
            break;
        case 2:
            fWeights[0] = 1;        fPoints(0, 0) = -0.57735026918962573;
            fWeights[1] = 1;        fPoints(1, 0) = +0.57735026918962573;
            break;
        case 3:
            fWeights[0] =  0.55555555555555558;     fPoints(0, 0) = -0.7745966692414834; 
            fWeights[1] =  0.88888888888888884;     fPoints(1, 0) =  0.0;                
            fWeights[2] =  0.55555555555555558;     fPoints(2, 0) = +0.7745966692414834; 
            break;
        case 4:
            fWeights[0] =  0.34785484513745385;     fPoints(0, 0) = -0.86113631159405257;
            fWeights[1] =  0.65214515486254609;     fPoints(1, 0) = -0.33998104358485626;
            fWeights[2] =  0.65214515486254609;     fPoints(2, 0) =  0.33998104358485626;
            fWeights[3] =  0.34785484513745385;     fPoints(3, 0) =  0.86113631159405257;
            break;
        case 5:
            fWeights[0] =  0.23692688505618908;     fPoints(0, 0) = -0.90617984593866396;
            fWeights[1] =  0.23692688505618908;     fPoints(1, 0) =  0.90617984593866396;
            fWeights[2] =  0.47862867049936647;     fPoints(2, 0) = -0.53846931010568311;
            fWeights[3] =  0.47862867049936647;     fPoints(3, 0) =  0.53846931010568311;
            fWeights[4] =  0.56888888888888889;     fPoints(4, 0) =  0.0;                
            break;
        case 6:
            fWeights[0] =  0.17132449237917036;     fPoints(0, 0) = -0.93246951420315205;
            fWeights[1] =  0.17132449237917036;     fPoints(1, 0) =  0.93246951420315205;
            fWeights[2] =  0.36076157304813861;     fPoints(2, 0) = -0.66120938646626448;
            fWeights[3] =  0.36076157304813861;     fPoints(3, 0) =  0.66120938646626448;
            fWeights[4] =  0.46791393457269104;     fPoints(4, 0) = -0.2386191860831969; 
            fWeights[5] =  0.46791393457269104;     fPoints(5, 0) =  0.2386191860831969; 
            break;
        case 7:
            fWeights[0] =  0.1294849661688697;      fPoints(0, 0) = -0.94910791234275849;
            fWeights[1] =  0.1294849661688697;      fPoints(1, 0) =  0.94910791234275849;
            fWeights[2] =  0.27970539148927664;     fPoints(2, 0) = -0.74153118559939446;
            fWeights[3] =  0.27970539148927664;     fPoints(3, 0) =  0.74153118559939446;
            fWeights[4] =  0.38183005050511892;     fPoints(4, 0) = -0.40584515137739718;
            fWeights[5] =  0.38183005050511892;     fPoints(5, 0) =  0.40584515137739718;
            fWeights[6] =  0.4179591836734694;      fPoints(6, 0) =  0.0;                
            break;
        case 8:
            fWeights[0] =  0.10122853629037626;     fPoints(0, 0) = -0.96028985649753618;
            fWeights[1] =  0.10122853629037626;     fPoints(1, 0) =  0.96028985649753618;
            fWeights[2] =  0.22238103445337448;     fPoints(2, 0) = -0.79666647741362673;
            fWeights[3] =  0.22238103445337448;     fPoints(3, 0) =  0.79666647741362673;
            fWeights[4] =  0.31370664587788727;     fPoints(4, 0) = -0.52553240991632899;
            fWeights[5] =  0.31370664587788727;     fPoints(5, 0) =  0.52553240991632899;
            fWeights[6] =  0.36268378337836199;     fPoints(6, 0) = -0.18343464249564981;
            fWeights[7] =  0.36268378337836199;     fPoints(7, 0) =  0.18343464249564981;
            break;
        default:DebugStop();
    }
}

void IntRule1d::gauleg(const double x1, const double x2, VecDouble &co, VecDouble &w){
    int n = w.size();

    double EPS = 1.0e-14;
    int m, j, i;
    double z1, z, xm, xl, pp, p3, p2, p1;    
    
    m = (n + 1) / 2;
    xm = 0.5 * (x2 + x1);
    xl = 0.5 * (x2 - x1);
    for (i = 0; i < m; i++) {
        z = cos(PI * (i + 0.75) / (n + 0.5));
        do {
            p1 = 1.0;
            p2 = 0.0;
            for (j = 0; j < n; j++) {
                p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1);
            }
            pp = n * (z * p1 - p2) / (z * z - 1.0);
            z1 = z;
            z = z1 - p1 / pp;
        } while (fabs(z - z1) > EPS);
        co[i] = xm - xl*z;
        co[n - 1 - i] = xm + xl*z;
        w[i] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
        w[n - 1 - i] = w[i];
    }
    
}

