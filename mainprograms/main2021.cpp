

#include <iostream>
#include <math.h>
#include "IntRule.h"
#include "IntRule1d.h"
#include "DataTypes.h"

using std::cout;
using std::endl;
using std::cin;

// f(x) = x
double func(double x){
    return sin(x);
}


int main ()
{
    // Create Quadrature
    IntRule1d regra(5);
    
    regra.Print(cout);

    // Compute integral
    double integral = 0.0;
    for(int i =0; i < regra.NPoints(); i++){
        VecDouble xi(1);
        double wi = 0.0;
        // IntRule::Point gives you the i-th point and weight
        regra.Point(i,xi,wi);
        // Integral is obtained through the sum
        integral += wi*func(xi[0]);
    }
    std::cout<<integral<<std::endl;
    return 0;
}
