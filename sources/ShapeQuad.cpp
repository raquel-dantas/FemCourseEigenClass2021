//
//  ShapeQuad.cpp
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include "Shape1d.h"
#include "ShapeQuad.h"

/// computes the shape functions in function of the coordinate in parameter space and orders of the shape functions (size of orders is number of sides of the element topology)
void ShapeQuad::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, MatrixDouble &dphi){
    
    for (int i = 0; i < orders.size(); i++)
    {
        if (orders[i] < 0) {
            std::cout << "ShapeQuad::Shape: Invalid dimension for arguments: order\n";
            DebugStop();
        }
    }
    if (orders[0] > 1 || orders[1] > 1 || orders[2] > 1 || orders[3] > 1) {
        std::cout << "ShapeQuad::Shape: Invalid dimension for arguments: order\n";
        DebugStop();
    }

    int nshape = NShapeFunctions(orders);
    phi.resize(nshape);
    dphi.resize(2,nshape);

    if (orders[nshape-1] > 2) {
        std::cout << "ShapeQuad::Shape, only implemented until order = 2" << std::endl;
        DebugStop();
    }
    // Linear
    phi[0] = 1./4.*(1. - xi[0])*(1. - xi[1]);
    phi[1] = 1./4.*(1. + xi[0])*(1. - xi[1]);
    phi[2] = 1./4.*(1. + xi[0])*(1. + xi[1]);
    phi[3] = 1./4.*(1. - xi[0])*(1. + xi[1]);

    // dphi(i,j) represents the ith derivative of function j
    dphi(0,0) = -1./4.*(1. - xi[1]);
    dphi(1,0) = -1./4.*(1. - xi[0]);
    dphi(0,1) = 1./4.*(1. - xi[1]);
    dphi(1,1) = -1./4.*(1. + xi[0]);
    dphi(0,2) = 1./4.*(1. + xi[1]);
    dphi(1,2) = 1./4.*(1. + xi[0]);
    dphi(0,3) = -1./4.*(1. + xi[1]);
    dphi(1,3) = 1./4.*(1. - xi[0]);

    // Quadratic
    if(nshape > nCorners)
    {
        int count = nCorners;
        for(int jside = nCorners; jside < nSides-1; jside++)
        {
            if(orders[jside] >= 2){
                double a = xi[jside%2];
                double b = xi[(jside+1)%2];
                double sign = ((jside+1)%nCorners > 1 ? +1. : -1.);
                phi[count] = 0.5*(1. - a*a)*(1 + sign*b);
                dphi(  jside  %2,count) = -a*(1. + sign*b);
                dphi((jside+1)%2,count) = sign*0.5*(1. - a*a);
                count++;
            }
        }
    }
    if(orders[nSides - 1] >= 2){
        phi[8] = (1. - xi[0]*xi[0])*(1. - xi[1]*xi[1]);
        dphi(0,8) = -2.*xi[0]*(1. - xi[1]*xi[1]);
        dphi(1,8) = -2.*xi[1]*(1. - xi[0]*xi[0]);
    }
}

/// returns the number of shape functions associated with a side
int ShapeQuad::NShapeFunctions(int side, int order){
    if(order < 1 || order >2) DebugStop();
    if(side<4)
        return 1;//0 a 4
    else if(side<8)
        return (order-1);//6 a 14
    else if(side==8)
        return ((order-1)*(order-1));
    
    std::cout << "ShapeQuad::NShapeFunctions, bad parameter side " << side << std::endl;
    DebugStop();
    
    return 0;
}

/// returns the total number of shape functions
int ShapeQuad::NShapeFunctions(VecInt &orders){
    
    int res=4;
    for(int in=4; in<orders.size(); in++) {
        res += NShapeFunctions(in, orders[in]);
    }
    
    return res;
}
