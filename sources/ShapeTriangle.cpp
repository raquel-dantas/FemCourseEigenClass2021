//
//  ShapeTriangle.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include "ShapeTriangle.h"
#include "Shape1d.h"

/// computes the shape functions in function of the coordinate in parameter space and orders of the shape functions (size of orders is number of sides of the element topology)
void ShapeTriangle::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, MatrixDouble &dphi){
    
    for (int i = 0; i < orders.size(); i++)
    {
        if (orders[i] < 0) {
            std::cout << "ShapeTriangle::Shape: Invalid dimension for arguments: order\n";
            DebugStop();
        }
    }
    if (orders[0] > 1 || orders[1] > 1 || orders[2] > 1) {
        std::cout << "ShapeTriangle::Shape: Invalid dimension for arguments: order\n";
        DebugStop();
    }

    int nshape = NShapeFunctions(orders);
    phi.resize(nshape);
    dphi.resize(2,nshape);

    if (orders[nshape-1] > 2) {
        std::cout << "ShapeTriangle::Shape, only implemented until order = 2" << std::endl;
        DebugStop();
    }
    
    // Linear order
    phi[0] =  1.-xi[0]-xi[1];
    phi[1] =  xi[0];
    phi[2] =  xi[1];

    // dphi(i,j) represents the ith derivative of function j
    dphi(0,0) = -1.;
    dphi(1,0) = -1.;
    dphi(0,1) =  1.;
    dphi(1,1) =  0.;
    dphi(0,2) =  0.;
    dphi(1,2) =  1.;
    
    // Quadratic
    if(nshape > nCorners)
    {
        int count = nCorners;
        for(int jside = nCorners; jside < nSides - 1; jside++)
        {
            int no_anterior = jside - nCorners;
            int no_posterior = (no_anterior + 1)%nCorners;
            if(orders[jside] >= 2)
            {
                phi[count] = 4*phi[no_anterior]*phi[no_posterior];
                // Product rule for derivatives
                dphi(0,count) = 4.*(dphi(0,no_anterior)*phi[no_posterior] + phi[no_anterior]*dphi(0,no_posterior));
                dphi(1,count) = 4.*(dphi(1,no_anterior)*phi[no_posterior] + phi[no_anterior]*dphi(1,no_posterior));
                count++;
            } 
        }
    }
    
    //Cubic
    if(orders[nSides - 1] >= 3){
        phi[6] = (orders[6] >= 2) * 27 *phi[0] * phi[1] * phi[2];
        dphi(0,6) = (orders[6] >= 2) * -27*xi[1]*(-1 + xi[1] + 2*xi[0]);
        dphi(1,6) = (orders[6] >= 2) * -27*xi[0]*(-1 + 2*xi[1] + xi[0]);
    }
    
}

/// returns the number of shape functions associated with a side
int ShapeTriangle::NShapeFunctions(int side, int order){
    switch(side) {
        case 0:
        case 1:
        case 2:
            return 1;
        case 3:
        case 4:
        case 5:
            return order-1;
        case 6:
            return 0;
    }
    
    DebugStop();
    std::cout << "ShapeTriangle::NShapeFunctions, bad parameter side " << std::endl;
    return 0;
}

/// returns the total number of shape functions
int ShapeTriangle::NShapeFunctions(VecInt &orders){
    
    int res=3;
    for(int in=3; in<orders.size(); in++) {
        res += NShapeFunctions(in, orders[in]);
    }
    
    return res;
    
}
