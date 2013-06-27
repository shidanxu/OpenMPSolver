//
//  global_problem.cpp
//  LPsolver
//
//  Created by Dragos Ciocan on 6/20/13.
//  Copyright (c) 2013 Dragos Ciocan. All rights reserved.
//

#include <algorithm>
#include <cmath>

#include "global_problem.h"
#include "instance.h"
#include "time.h"
#include <omp.h>

namespace distributed_solver {
    GlobalProblem::GlobalProblem(int num_partitions, long double max_bid, long double advertiser_indegree,
                                 long double numerical_accuracy_tolerance,
                                 vector<vector<pair<int, pair<long double, long double> > > >* solution) {
        num_partitions_ = num_partitions;
        budget_ = 0;
        solution_ = solution;
        num_iterations_ = 0;
        numerical_accuracy_tolerance_ = numerical_accuracy_tolerance;

        for (int i = 0; i < num_partitions; ++i) {
            primal_changes_.push_back(make_pair(make_pair(-1, 0.0), make_pair(-1, 0.0)));
        }
    }

    void GlobalProblem::FindOptimalBudgetAllocation() {
        long double remaining_budget = budget_;

        // Create list of slopes.
        vector<Slope> slopes;
        for (int i = 0; i < num_partitions_; ++i) {
            // Reset budget allocation.
            budget_allocation_[i].second = 0;

            if (subproblems_[i].envelope_points_.size() > 1) {
                // -1 because at last envelope return on budget is 0.
                for (int j = 0; j < subproblems_[i].envelope_points_.size() - 1; ++j) {
                    slopes.push_back(Slope(subproblems_[i].envelope_points_[j].first, i, j));
                }
            }
        }
        // Sort list of slopes.
        sort(slopes.begin(), slopes.end(), compare_Slope());

        // Allocate budgets in decreasing order of slopes.
        int sp_index;
        for (int k = 0; k < slopes.size(); ++k) {
            sp_index = slopes[k].subproblem_index_;
            long double allocation_increase = min(remaining_budget,
                                                  subproblems_[sp_index].budget_cutoffs_[slopes[k].region_index_ + 1] -
                                                  subproblems_[sp_index].budget_cutoffs_[slopes[k].region_index_]);
            if ((allocation_increase > 0) && (remaining_budget > 0)) {
                budget_allocation_[sp_index] = make_pair(slopes[k].region_index_, budget_allocation_[sp_index].second + allocation_increase);
            }
            remaining_budget = remaining_budget - allocation_increase;
        }
    }

    void GlobalProblem::ConstructPrimal(int iteration) {
        num_iterations_++;

        // Reset primal solution.
        Instance::ResetCurrentPrimal(solution_);

        clock_t t1, t2;
        float diff;
        t1 = clock();
        // Find optimal budget allocations to problems.
        cout << "Solving subproblems \n";

        /*
        for (int i = 0; i < num_partitions_; ++i) {
            //subproblems_[i].SolveSubproblem(iteration, i);
            //subproblems_[i].SolveSubproblemConvexHull(iteration, i);
            subproblems_[i].SolveSubproblemConvexHullOptimized(iteration, i);
        }
         */

        int num_threads = 4;
        int threadId = 0;
        int i, ending;

        omp_set_num_threads(num_threads);

        #pragma omp parallel private ( threadId, i, ending ) \
        shared( iteration )

		{
			threadId = omp_get_thread_num();
			ending = (threadId + 1 ) * num_partitions_ / num_threads;

			for (i = threadId * num_partitions_ / num_threads; i < ending; i++){
				subproblems_[i].SolveSubproblemConvexHullOptimized(iteration, i);
			}
		}

        t2 = clock();
        diff = ((float)t2-(float)t1);
        cout << "subproblems took  " << diff << "\n";


        cout << "Find optimal budget allocation \n";
        t1 = clock();
        FindOptimalBudgetAllocation();
        t2 = clock();
        diff = ((float)t2-(float)t1);
        cout << "budget allocation took  " << diff << "\n";

        // Calculate primal solution for each subproblem.
        cout << "Constructing primal \n";
        t1 = clock();
        long double dual_val = 0;
        primal_assignment_test_ = 0;

        /*
        for (int i = 0; i < num_partitions_; ++i) {
            if (budget_allocation_[i].second > 0) {
                long double u = subproblems_[i].envelope_points_[budget_allocation_[i].first].first;
                long double v = subproblems_[i].envelope_points_[budget_allocation_[i].first].second;
                dual_val += u * budget_allocation_[i].second + v;
                ConstructSubproblemPrimal(i, budget_allocation_[i].second, budget_allocation_[i].first);
            }
        }
        */

        threadId = 0;

        omp_set_num_threads(num_threads);

        #pragma omp parallel private ( threadId, ending, i ) \
        shared( dual_val )
        {
            threadId = omp_get_thread_num();
            cout << omp_get_num_threads();
            ending = (threadId + 1 ) * num_partitions_ / num_threads;
            for (i = threadId * num_partitions_ / num_threads; i < ending; i++){
                long double u = subproblems_[i].envelope_points_[budget_allocation_[i].first].first;
                long double v = subproblems_[i].envelope_points_[budget_allocation_[i].first].second;
                dual_val += u * budget_allocation_[i].second + v;
                ConstructSubproblemPrimal(i, budget_allocation_[i].second, budget_allocation_[i].first);
            }

        }

        Instance::UpdatePrimal(iteration, solution_, primal_changes_);

        t2 = clock();
        diff = ((float)t2-(float)t1);
        cout << "constructing primal took  " << diff << "\n";

        cout << "Dual Value = ";
        cout << dual_val;
        cout << "\n";
    }

    void GlobalProblem::ConstructSubproblemPrimal(int subproblem_index, long double budget_allocation,
                                                  int opt_region) {
        // Figure out opt u, v.
        long double u = subproblems_[subproblem_index].envelope_points_[opt_region].first;
        long double v = subproblems_[subproblem_index].envelope_points_[opt_region].second;

        long double allocation_value = 0;

        // If optimum is u = 0, optimal allocation is greedy wrt price.
        if (u == 0) {
            int max_price_index;
            long double max_price = 0;
            for (int i = 0; i < subproblems_[subproblem_index].num_vars_; ++i) {
                if ((max_price < subproblems_[subproblem_index].constraints_[i].price_) and (subproblems_[subproblem_index].constraints_[i].is_active_)) {
                    max_price_index = i;
                    max_price = subproblems_[subproblem_index].constraints_[i].price_;
                }
            }

            //(*primal_sol)[subproblems_[subproblem_index].advertiser_index_->at(max_price_index)][subproblem_index] = 1;
            //(*solution_)[subproblems_[subproblem_index].constraints_[max_price_index].advertiser_index_][subproblem_index].first = 1;
            //(*primal_changes)[subproblems_[subproblem_index].constraints_[max_price_index].advertiser_index_][subproblem_index] = 1;
            primal_changes_[subproblem_index].first = make_pair(subproblems_[subproblem_index].constraints_[max_price_index].advertiser_index_, 1.0);
            primal_changes_[subproblem_index].second = make_pair(-1, 0.0);
            allocation_value = primal_changes_[subproblem_index].first.second *
                               subproblems_[subproblem_index].constraints_[max_price_index].price_;
            primal_assignment_test_ += allocation_value;
        }

        // If optimum is v = 0, optimal allocation is greedy wrt price/weight ratio.
        if (v == 0) {
            int max_ratio_index;
            long double max_ratio = 0;
            for (int i = 0; i < subproblems_[subproblem_index].num_vars_; ++i) {
                if (max_ratio < (subproblems_[subproblem_index].constraints_[i].price_ /
                                 subproblems_[subproblem_index].constraints_[i].coefficient_) and
                    (subproblems_[subproblem_index].constraints_[i].is_active_)) {
                    max_ratio_index = i;
                    max_ratio = (subproblems_[subproblem_index].constraints_[i].price_ /
                                 subproblems_[subproblem_index].constraints_[i].coefficient_);
                }
            }

            //(*primal_changes)[subproblems_[subproblem_index].constraints_[max_ratio_index].advertiser_index_][subproblem_index] = budget_allocation / subproblems_[subproblem_index].constraints_[max_ratio_index].coefficient_;
            //(*solution_)[subproblems_[subproblem_index].constraints_[max_ratio_index].advertiser_index_][subproblem_index].first = budget_allocation / subproblems_[subproblem_index].constraints_[max_ratio_index].coefficient_;
            //(*primal_sol)[subproblems_[subproblem_index].advertiser_index_->at(max_ratio_index)][subproblem_index] = budget_allocation / subproblems_[subproblem_index].constraints_[max_ratio_index].coefficient_;

            primal_changes_[subproblem_index].first = make_pair(subproblems_[subproblem_index].constraints_[max_ratio_index].advertiser_index_, budget_allocation / subproblems_[subproblem_index].constraints_[max_ratio_index].coefficient_);
            primal_changes_[subproblem_index].second = make_pair(-1, 0.0);
            allocation_value = primal_changes_[subproblem_index].first.second * subproblems_[subproblem_index].constraints_[max_ratio_index].price_;
            primal_assignment_test_ += allocation_value;
        }

        // If opt is u, v > 0, the optimum whp only has two positive allocations, which are the
        // solutions to a system of 2 equations.
        if ((u > 0) and (v > 0)) {
            vector<int> tight_constraint_indices;
            for (int i = 0; i < subproblems_[subproblem_index].num_vars_; ++i) {
                if (subproblems_[subproblem_index].constraints_[i].is_active_) {
                    long double slack = subproblems_[subproblem_index].constraints_[i].price_ -
                    (u * subproblems_[subproblem_index].constraints_[i].coefficient_ + v);
                    if (slack < 0) { slack = (-1) * slack;}
                    if (slack < numerical_accuracy_tolerance_) {
                        tight_constraint_indices.push_back(i);
                    }
                }
            }
            if (tight_constraint_indices.size() > 2) {
                cout << "ERROR, PERTURB PRICES \n";
            }
            if (tight_constraint_indices.size() == 1) {
                /*
                 (*primal_sol)[subproblems_[subproblem_index].advertiser_index_->
                 at(tight_constraint_indices[0])]
                 [subproblem_index] = fmin(budget_allocation / subproblems_[subproblem_index].constraints_[tight_constraint_indices[0]].coefficient_, 1);
                 */
                //(*solution_)[subproblems_[subproblem_index].constraints_[tight_constraint_indices[0]].advertiser_index_][subproblem_index].first = fmin(budget_allocation / subproblems_[subproblem_index].constraints_[tight_constraint_indices[0]].coefficient_, 1);
                //(*primal_changes)[subproblems_[subproblem_index].constraints_[tight_constraint_indices[0]].advertiser_index_][subproblem_index] = fmin(budget_allocation / subproblems_[subproblem_index].constraints_[tight_constraint_indices[0]].coefficient_, 1);
                primal_changes_[subproblem_index].first = make_pair(subproblems_[subproblem_index].constraints_[tight_constraint_indices[0]].advertiser_index_, fmin(budget_allocation / subproblems_[subproblem_index].constraints_[tight_constraint_indices[0]].coefficient_, 1));
                primal_changes_[subproblem_index].second = make_pair(-1, 0.0);
                allocation_value = primal_changes_[subproblem_index].first.second * subproblems_[subproblem_index].constraints_[tight_constraint_indices[0]].coefficient_;
                primal_assignment_test_ += allocation_value;
            }
            if (tight_constraint_indices.size() == 2) {
                int first_index = tight_constraint_indices[0];
                int second_index = tight_constraint_indices[1];

                long double x_1 = (budget_allocation -
                                   subproblems_[subproblem_index].constraints_[second_index].coefficient_) / (subproblems_[subproblem_index].constraints_[first_index].coefficient_ - subproblems_[subproblem_index].constraints_[second_index].coefficient_);
                /*
                 (*primal_sol)[subproblems_[subproblem_index].advertiser_index_->at(first_index)][subproblem_index] = x_1;
                 (*primal_sol)[subproblems_[subproblem_index].advertiser_index_->at(second_index)][subproblem_index] = 1 - x_1;
                 */
                /*
                (*solution_)[subproblems_[subproblem_index].constraints_[first_index].advertiser_index_][subproblem_index].first = x_1;
                (*solution_)[subproblems_[subproblem_index].constraints_[second_index].advertiser_index_][subproblem_index].first = 1 - x_1;
                */
                /*
                (*primal_changes)[subproblems_[subproblem_index].constraints_[first_index].advertiser_index_][subproblem_index] = x_1;
                (*primal_changes)[subproblems_[subproblem_index].constraints_[second_index].advertiser_index_][subproblem_index] = 1 - x_1;
                */
                primal_changes_[subproblem_index].first = make_pair(subproblems_[subproblem_index].constraints_[first_index].advertiser_index_, x_1);
                primal_changes_[subproblem_index].second = make_pair(subproblems_[subproblem_index].constraints_[second_index].advertiser_index_, 1 - x_1);

                 allocation_value = x_1 * subproblems_[subproblem_index].constraints_[first_index].price_ + (1-x_1) * subproblems_[subproblem_index].constraints_[second_index].price_;
                primal_assignment_test_ += allocation_value;

            }
        }
        if (((allocation_value - (u * budget_allocation + v)) > numerical_accuracy_tolerance_) ||
            ((allocation_value - (u * budget_allocation + v)) < -numerical_accuracy_tolerance_)) {
            cout << "**************Error at problem************ " << subproblem_index << " equal to " <<
            allocation_value << " - " << (u * budget_allocation + v) << "\n";
        }
    }

    void GlobalProblem::InitializeBudgetAllocation() {
        for (int i = 0; i < num_partitions_; ++i) {
            budget_allocation_.push_back(make_pair(0, 0.0));
        }
    }

    Slope::Slope(long double slope, int subproblem_index, int region_index) {
        slope_ = slope;
        subproblem_index_ = subproblem_index;
        region_index_ = region_index;
    }

    long double GlobalProblem::FindOptimalBudgetAllocationBinSearch(long double lower, long double upper) {
        vector<long double>* budget_usage = new vector<long double>(num_partitions_);
        long double lower_bound = lower;
        long double upper_bound = upper;
        long double critical_ratio = (lower_bound + upper_bound) / 2.0;
        long double critical_ratio_tolerance = 0.0001;

        //figure out (budget_usage - budget_) for current critical ratio
        long double delta = CalculateAllocationDelta(critical_ratio, budget_usage);

        if (delta < 0) {
            // Critical ratio too high.
            upper_bound = critical_ratio;
        } else if (delta > 0) {
            // Critical ratio too low.
            lower_bound = critical_ratio;
        } else {
            // Budget perfectly assigned, allocate as is.
            AllocateCurrentRatio(critical_ratio);
            return critical_ratio;
        }

        if (upper_bound - lower_bound < critical_ratio_tolerance) {
            // Terminate if critical ratio is closely bounded.
            AllocateCurrentRatio(critical_ratio);
            return critical_ratio;
        }
        else {
            // Recurse.
            return FindOptimalBudgetAllocationBinSearch(lower_bound, upper_bound);
        }
    }

    long double GlobalProblem::CalculateAllocationDelta(long double critical_ratio, vector<long double>* budget_usage) {
        long double total_budget_usage = 0.0;
        for (int i = 0; i < num_partitions_; ++i) {
            (*budget_usage)[i] = 0;
            for (int j = 0; j < subproblems_[i].envelope_points_.size() - 1; ++j) {
                if (subproblems_[i].envelope_points_[j].first >= critical_ratio) {
                    (*budget_usage)[i] += subproblems_[i].budget_cutoffs_[j + 1] - subproblems_[i].budget_cutoffs_[j];
                }
            }
            total_budget_usage += (*budget_usage)[i];
        }
        return total_budget_usage - budget_;
    }

    void GlobalProblem::AllocateCurrentRatio(long double critical_ratio) {
        // This method updates optimal budget_allocation_ based on critical ratio method.
        for (int i = 0; i < num_partitions_; ++i) {
            for (int j = 0; j < subproblems_[i].envelope_points_.size() - 1; ++j) {
                if (subproblems_[i].envelope_points_[j].first >= critical_ratio) {
                    budget_allocation_[i] = make_pair(j,
                                                      budget_allocation_[i].second
                                                      + subproblems_[i].budget_cutoffs_[j + 1]
                                                      - subproblems_[i].budget_cutoffs_[j]);
                }
            }
        }
    }
}
