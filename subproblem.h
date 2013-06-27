//
//  subproblem.h
//  LPsolver
//
//  Created by Dragos Ciocan on 6/20/13.
//  Copyright (c) 2013 Dragos Ciocan. All rights reserved.
//

#ifndef __LPsolver__subproblem__
#define __LPsolver__subproblem__

#include <ext/hash_map>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

namespace distributed_solver {
    class Constraint {
    public:
        long double price_;
        long double coefficient_;
        long double weight_; // Add this to avoid arithmetic precision errors.
        bool is_active_;
        int advertiser_index_;
        Constraint(); // Default constructor.
        Constraint(long double price, long double coefficient, long double weight, int index);
        void set_active(bool value);
    };

    struct compare_Constraint_by_weight
    {
        bool operator() (const Constraint & lhs, const Constraint & rhs) {
            return (lhs.weight_ < rhs.weight_);
        }
    };

    class Subproblem {
    public:
        int num_vars_;
        std::vector<Constraint> constraints_;
        std::vector<std::pair<long double, long double> > envelope_points_;
        std::vector<long double> budget_cutoffs_;
        Subproblem(int num_vars, std::vector<std::pair<long double, long double> >* coefficients, std::vector<int>* advertiser_index);
        void SolveSubproblem(int iteration, int index);
        void SolveSubproblemConvexHull(int iteration, int index);
        void SolveSubproblemConvexHullOptimized(int iteration, int index);
    };
}

#endif /* defined(__LPsolver__subproblem__) */
