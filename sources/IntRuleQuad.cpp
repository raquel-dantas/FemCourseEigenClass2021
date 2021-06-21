/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

///\cond
#include <iostream> 
///\endcond
#include "IntRule1d.h"
#include "IntRuleQuad.h"

IntRuleQuad::IntRuleQuad(){}

IntRuleQuad::IntRuleQuad(int order) {
    // if(order<0||order>=6) DebugStop();
    SetOrder(order);
}

void IntRuleQuad::SetOrder(int order) {
    if(order<0||order>=6) DebugStop();

    int linepoints = (order+1)/2 + (order+1)%2;
    int npoints = linepoints*linepoints;
    fOrder = order;
    if(npoints==this->NPoints())return;

    fWeights.resize(npoints);
    fPoints.resize(npoints,2);

    IntRule1d rule_line(order);
    
    double wx = 0.;
    double wy = 0.;
    VecDouble x(1);
    VecDouble y(1);

    for(int i = 0; i < linepoints; i++){
        rule_line.Point(i,x,wx);
        for(int j = 0; j < linepoints; j++){
            rule_line.Point(j,y,wy);

            fWeights[linepoints*i + j] = wx * wy;
            fPoints(linepoints*i + j,0) = x[0];
            fPoints(linepoints*i + j,1) = y[0];
        }

    }



}

void IntRuleQuad::gaulegQuad(const double x1, const double x2, VecDouble &co, VecDouble &w) {
    IntRule1d x;
    IntRule1d y;
    
    int n = w.size();   

    VecDouble cox(n);
    VecDouble coy(n);
    VecDouble wx(n);
    VecDouble wy(n);


    x.gauleg(x1, x2, cox, wx);
    y.gauleg(x1, x2, coy, wy);
    
    co.resize(2*n*n);
    w.resize(n * n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            co[j + i * n] = cox[j];
            co[j + i * n + n * n] = coy[i];
            w[n * i + j] = wx[i] * wy[j];
        }
    }
}
