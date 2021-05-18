/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <iostream> 
#include "IntRuleTriangle.h"
#include "tpanic.h"

IntRuleTriangle::IntRuleTriangle(){}

IntRuleTriangle::IntRuleTriangle(int order) {
    SetOrder(order);
}

void IntRuleTriangle::SetOrder(int order) {
    if(order<0 || order > gMaxOrder()) DebugStop();
                                // 0 1 2 3 4 5
    const int map_order_npoint[6] {1,1,3,4,6,7};
    int npoints = map_order_npoint[order];

    fOrder = order;
    if(npoints == this->NPoints()) return;
    fPoints.resize(npoints,2);
    fWeights.resize(npoints);

    // Taken from https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
    switch (npoints){
        case 1: // order up to 1
            fWeights[0] = 1./2.;    fPoints(0,0) = 1./3.;
                                    fPoints(0,1) = 1./3.;
            break;
        case 3: // order up to 2
            fWeights[0] = 1./6.;    fPoints(0,0) = 1./2.;
                                    fPoints(0,1) = 1./2.;
            fWeights[1] = 1./6.;    fPoints(1,0) = 0.;
                                    fPoints(1,1) = 1./2.;
            fWeights[2] = 1./6.;    fPoints(2,0) = 1./2.;
                                    fPoints(2,1) = 0.;
            break;
        case 4: // order up to 3
            fWeights[0] = -27./96.; fPoints(0,0) = 1./3.;
                                    fPoints(0,1) = 1./3.;
            fWeights[1] =  25./96.; fPoints(1,0) = 1./5.;
                                    fPoints(1,1) = 1./5.;
            fWeights[2] =  25./96.; fPoints(2,0) = 3./5.;
                                    fPoints(2,1) = 1./5.;
            fWeights[3] =  25./96.; fPoints(3,0) = 1./5.;
                                    fPoints(3,1) = 3./5.;
            break;
        case 6: // order up to 4
            fWeights[0] =  0.11169079483900574;     fPoints(0,0) = 0.44594849091596489;  
                                                    fPoints(0,1) = 0.10810301816807023;
            fWeights[1] =  0.11169079483900574;     fPoints(1,0) = 0.10810301816807023;  
                                                    fPoints(1,1) = 0.44594849091596489;
            fWeights[2] =  0.11169079483900574;     fPoints(2,0) = 0.44594849091596489;  
                                                    fPoints(2,1) = 0.44594849091596489;
            fWeights[3] =  0.054975871827660935;    fPoints(3,0) = 0.091576213509770743;  
                                                    fPoints(3,1) = 0.81684757298045851;
            fWeights[4] =  0.054975871827660935;    fPoints(4,0) = 0.81684757298045851;  
                                                    fPoints(4,1) = 0.091576213509770743;
            fWeights[5] =  0.054975871827660935;    fPoints(5,0) = 0.091576213509770743;  
                                                    fPoints(5,1) = 0.091576213509770743;
            break;
        case 7: // order up to 5
            fWeights[0] = 0.06296959027241357;      fPoints(0,0) = 0.10128650732345634;  
                                                    fPoints(0,1) = 0.79742698535308731;
            fWeights[1] = 0.06296959027241357;      fPoints(1,0) = 0.79742698535308731;  
                                                    fPoints(1,1) = 0.10128650732345634;
            fWeights[2] = 0.06296959027241357;      fPoints(2,0) = 0.10128650732345634;  
                                                    fPoints(2,1) = 0.10128650732345634;
            fWeights[3] = 0.066197076394253096;     fPoints(3,0) = 0.47014206410511511;  
                                                    fPoints(3,1) = 0.059715871789769823;
            fWeights[4] = 0.066197076394253096;     fPoints(4,0) = 0.059715871789769823;  
                                                    fPoints(4,1) = 0.47014206410511511;
            fWeights[5] = 0.066197076394253096;     fPoints(5,0) = 0.47014206410511511;  
                                                    fPoints(5,1) = 0.47014206410511511; 
            fWeights[6] = 0.112500000000000000;     fPoints(6,0) = 0.33333333333333333;  
                                                    fPoints(6,1) = 0.33333333333333333;
            break;
        default:
            DebugStop();
            break;
    }

}
