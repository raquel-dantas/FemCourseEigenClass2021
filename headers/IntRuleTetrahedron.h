//
//  IntRuleTriangle.h
//  FemSC
//
//  Created by Philippe Devloo on 7/30/15.
//
//

#ifndef __FemSC__IntRuleTetrahedron__
#define __FemSC__IntRuleTetrahedron__

///\cond
#include <stdio.h>
///\endcond
#include "IntRule.h"

/**
@brief Integration rule associated with a tetrahedron
@ingroup integration
*/
class IntRuleTetrahedron : public IntRule
{
  
    public:
    
    // Default Constructor of integration rule for tetrahedron elements
    IntRuleTetrahedron();
    
    // Constructor of integration rule for tetrahedron elements
    IntRuleTetrahedron(int order);

    // Dimension of the integration rule
    virtual int Dimension() const override{
        return 3;
    }
    
    // Return the maximum polynomial order that can be integrated exactly
    static int gMaxOrder()
    {
        return 14;
    }

    // Return the maximum polynomial order that can be integrated exactly
    virtual int MaxOrder() const override
    {
        return gMaxOrder();
    }
    // Method to set polynomial order of the integration rule for tetrahedro elements
    virtual void SetOrder(int order) override;
    
};


#endif /* defined(__FemSC__TIntRule1d__) */
