//
//  global_problem.h
//  LPsolver
//
//  Created by Dragos Ciocan on 6/20/13.
//  Copyright (c) 2013 Dragos Ciocan. All rights reserved.
//

#ifndef __LPsolver__global_problem__
#define __LPsolver__global_problem__

#include <cmath>
#include <iostream>
#include <vector>

#include "subproblem.h"

using namespace std;

namespace distributed_solver {
    class GlobalProblem {
    public:
        int num_partitions_;
        long double budget_;
        int num_iterations_;

        vector<pair<int, long double> > budget_allocation_;
        vector<Subproblem> subproblems_;
        vector<long double> slacks_;

        vector<vector<pair<int, pair<long double, long double> > > >* solution_;

        vector<pair<pair<int, long double>, pair<int, long double> > > primal_changes_;

        GlobalProblem(int num_partitions, long double max_bid, long double advertiser_indegree,
                      long double numerical_accuracy_tolerance,
                      vector<vector<pair<int, pair<long double, long double> > > >* solution);
        void InitializeInstance();
        void InitializeBudgetAllocation();
        void ConstructPrimal(int iteration);

    private:
        void FindOptimalBudgetAllocation();
        void ConstructSubproblemPrimal(int subproblem_index, long double budget_allocation, int opt_region);
        long double numerical_accuracy_tolerance_;
        long double primal_assignment_test_;
        long double FindOptimalBudgetAllocationBinSearch(long double lower, long double upper);
        long double CalculateAllocationDelta(long double critical_ratio, vector<long double>* budget_usage);
        void AllocateCurrentRatio(long double critical_ratio);
    };

    class Slope {
    public:
        long double slope_;
        int subproblem_index_;
        int region_index_;
        Slope(long double slope, int subproblem_index, int region_index);
    };

    struct compare_Slope
    {
        bool operator() (const Slope & lhs, const Slope & rhs) {
            return lhs.slope_ > rhs.slope_;
        }
    };

}

#endif /* defined(__LPsolver__global_problem__) */
