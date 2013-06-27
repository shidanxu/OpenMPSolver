//
//  envelope.h
//  LPsolver
//
//  Created by Dragos Ciocan on 6/20/13.
//  Copyright (c) 2013 Dragos Ciocan. All rights reserved.
//

#ifndef __LPsolver__envelope__
#define __LPsolver__envelope__

#include <iostream>

#include "subproblem.h"

// Implementation of Andrew's monotone chain 2D convex hull algorithm.
// Asymptotic complexity: O(n log n).
// Practical performance: 0.5-1.0 seconds for n=1000000 on a 1GHz machine.
#include <algorithm>
#include <vector>
using namespace std;

namespace distributed_solver {
    typedef long double coord_t;         // coordinate type
    typedef long double coord2_t;  // must be big enough to hold 2*max(|coordinate|)^2

    struct compare_Constraint_lexicographically {
        bool operator() (const Constraint &a, const Constraint &b) const {
            return ((a.coefficient_ > b.coefficient_) ||
                    (a.coefficient_ == b.coefficient_ && a.price_ > b.price_));
        }
    };

    // Returns a list of points on the convex hull in counter-clockwise order.
    // Note: the last point in the returned list is the same as the first one.
    vector<Constraint> upper_envelope(vector<Constraint>* P, long double numerical_tolerance);
}

#endif /* defined(__LPsolver__envelope__) */
