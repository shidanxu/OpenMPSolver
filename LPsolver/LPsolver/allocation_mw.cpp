//
//  allocation_mw.cpp
//  LPsolver
//
//  Created by Dragos Ciocan on 6/20/13.
//  Copyright (c) 2013 Dragos Ciocan. All rights reserved.
//

#include <cmath>

#include "allocation_mw.h"
#include "convex_hull.h"
#include "instance.h"
#include "time.h"


namespace distributed_solver {
    typedef IloArray<IloNumArray>    NumMatrix;
    typedef IloArray<IloNumVarArray> NumVarMatrix;
    
    AllocationMW::AllocationMW(int num_advertisers, int num_impressions, int num_slots, long double bid_sparsity,
                               long double max_bid, long double epsilon, long double numerical_accuracy_tolerance,
                               vector<vector<pair<int, long double> > >* bids_matrix,
                               vector<vector<pair<int, long double> > >* transpose_bids_matrix,
                               vector<long double>* budgets,
                               vector<vector<pair<int, pair<long double, long double> > > >* solution) :
    global_problem_(num_impressions, max_bid, bid_sparsity * num_impressions, numerical_accuracy_tolerance, solution) {
        
        epsilon_ = epsilon;
        iteration_count_ = 0;
        max_bid_ = max_bid;
        numerical_accuracy_tolerance_ = numerical_accuracy_tolerance;
        
        num_advertisers_ = num_advertisers;
        num_impressions_ = num_impressions;
        num_slots_ = num_slots;
        bid_sparsity_ = bid_sparsity;
        num_shards_ = 10;
        
        budgets_ = budgets;
        bids_matrix_ = bids_matrix;
        transpose_bids_matrix_ = transpose_bids_matrix;
        
        solution_ = solution;
        
        weights_ = vector<long double>(num_advertisers_);
        slacks_ = vector<long double>(num_advertisers_);
        avg_slacks_ = vector<long double>(num_advertisers_);
        for (int i = 0; i < num_advertisers_; ++i) {
            weights_[i] = 1;
            // Perturbation against degeneracy
            // weights_[i] = 1 + ((long double) (rand() + 1) / ((long double) RAND_MAX)) / 10000;
            slacks_[i] = 0;
        }
        
        CalculateAllocationMWWidth();
        CreateGlobalProblem();
        global_problem_.InitializeBudgetAllocation();
    }
    
    void AllocationMW::ComputeWeightedBudget() {
        // Compute budget.
        global_problem_.budget_ = 0;
        for (int i = 0; i < num_advertisers_; ++i) {
            global_problem_.budget_ = global_problem_.budget_ + weights_[i] * (*budgets_)[i];
        }
    }
    
    void AllocationMW::CreateGlobalProblem() {
        ComputeWeightedBudget();
        for (int i = 0; i < global_problem_.num_partitions_; ++i) {
            // Construct subproblems.
            vector<pair<long double, long double> >* coefficients = new vector<pair<long double, long double> >();
            vector<int>* advertiser_index = new vector<int>();
            
            // Go through bids matrix and identify all bids for impression i.
            int subproblem_size = 0;
            
            for (int j = 0; j < (*transpose_bids_matrix_)[i].size(); ++j) {
                coefficients->push_back(make_pair((*transpose_bids_matrix_)[i][j].second,
                                                  weights_[(*transpose_bids_matrix_)[i][j].first]));
                advertiser_index->push_back((*transpose_bids_matrix_)[i][j].first);
                ++subproblem_size;
            }
            
            global_problem_.subproblems_.push_back(Subproblem(subproblem_size, coefficients, advertiser_index));
            (*coefficients).clear();
            (*advertiser_index).clear();
        }
        (*transpose_bids_matrix_).clear();
        // cout << "Size of global problem object is " << (sizeof(global_problem_) * 4.0) / (1024.0 * 1024.0) << "\n";
    }
    
    void AllocationMW::UpdateGlobalProblem() {
        // Update weights.
        for (int i = 0; i < global_problem_.num_partitions_; ++i) {
            for (int j = 0; j < global_problem_.subproblems_[i].num_vars_; ++j) {
                global_problem_.subproblems_[i].constraints_[j].coefficient_ = global_problem_.subproblems_[i].constraints_[j].price_ * weights_[global_problem_.subproblems_[i].constraints_[j].advertiser_index_];
                global_problem_.subproblems_[i].constraints_[j].weight_ = weights_[global_problem_.subproblems_[i].constraints_[j].advertiser_index_];
                global_problem_.subproblems_[i].constraints_[j].is_active_ = true;
            }
        }
        ComputeWeightedBudget();
    }
    
    void AllocationMW::CalculateAllocationMWWidth() {
        width_ = max_bid_ * (num_impressions_ * bid_sparsity_);
        for (int i = 0; i < num_advertisers_; ++i) {
            if ((*budgets_)[i] > width_) {
                width_ = (*budgets_)[i];
            }
        }
    }
    
    void AllocationMW::CalculateSlacks() {
        for (int j = 0; j < num_advertisers_; ++j) {
            slacks_[j] = (-1) * (*budgets_)[j];
            for (int i = 0; i < (*solution_)[j].size(); ++i) {
                slacks_[j] += (*solution_)[j][i].second.second * (*bids_matrix_)[j][i].second;
                //slacks_[j] = slacks_[j] + iter->second.first * (*bids_matrix_)[j].find(iter->first)->second;
            }
        }
    }
    
    void AllocationMW::UpdateWeights() {
        for (int j = 0; j < num_advertisers_; ++j) {
            long double tmp = slacks_[j] / width_;
            if (abs(tmp) > 1) {
                cout << "Slack normalization error\n";
            }
            
            if (tmp >= 0) {
                weights_[j] = weights_[j] * pow((1 + epsilon_), tmp);
            }
            else {
                weights_[j] = weights_[j] * pow((1 - epsilon_), -tmp);
            }
            
            // weights_[j] = weights_[j] * (1 + epsilon_ * tmp);
            // cout << "weight " << j << " set to " << weights_[j] << "\n";
        }
    }
    
    void AllocationMW::SetBudgets(long double scaling_factor) {
        long double average_bid = 0.5; // Need to change this manually depending on how bids are drawn.
        for (int j = 0; j < num_advertisers_; ++j) {
            (*budgets_)[j] = average_bid * (num_impressions_ / num_advertisers_) * scaling_factor;
        }
    }
    
    void AllocationMW::UpdateAvgSlacks(int t) {
        for (int i = 0; i < num_advertisers_; ++i) {
            avg_slacks_[i] = (long double) (t - 1) / t * avg_slacks_[i] + (long double) 1 / t * slacks_[i];
        }
    }
    
    void AllocationMW::ReportWorstInfeasibility(int t) {
        long double max_infeasibility = 0.0;
        int max_infeasibility_index = -1;
        for (int i = 0; i < num_advertisers_; ++i) {
            if ((avg_slacks_[i] > 0) and ((avg_slacks_[i] / (*budgets_)[i]) > max_infeasibility)){
                max_infeasibility = avg_slacks_[i] / (*budgets_)[i];
                max_infeasibility_index = i;
            }
        }
        cout << "At iteration ";
        cout << t;
        cout << ", max infeasiblity was ";
        cout << max_infeasibility;
        cout << " on constraint ";
        cout << max_infeasibility_index;
        cout << "\n";
        
        ofstream report_file;
        report_file.open("/Users/ciocan/Documents/Google/data/infeasibility_report.csv", ios::app);
        report_file << t;
        report_file << ", ";
        report_file << max_infeasibility;
        report_file << ", ";
        report_file << max_infeasibility_index;
        report_file << "\n";
        report_file.close();
    }
    
    void AllocationMW::VerifySolution() {
        IloEnv env;
        IloInt i, j;
        IloModel model(env);
        IloInt nbImpressions = num_impressions_;
        IloInt nbAdvertisers = num_advertisers_;
        
        NumVarMatrix x(env, nbImpressions);
        // NumVarMatrix y(env, nbImpressions);
        
        for(i = 0; i < nbImpressions; i++) {
            x[i] = IloNumVarArray(env, nbAdvertisers, 0.0, 1.0, ILOFLOAT);
            // y[i] = IloNumVarArray(env, nbAdvertisers, 0.0, 1.0, ILOFLOAT);
        }
        
        // Add impression constraints.
        for(i = 0; i < nbImpressions; i++) {
            model.add(IloSum(x[i]) <= 1.0);
        }
        
        // Add weighted contraint.
        IloExpr weighted_constraint(env);
        for(j = 0; j < nbImpressions; j++) {      // demand must meet supply
            for(i = 0; i < nbAdvertisers; i++) {
                weighted_constraint += ((double) (*transpose_bids_matrix_)[j][i].second) *
                                       ((double) weights_[(int) (*transpose_bids_matrix_)[j][i].first]) * x[j][i];
            }
        }
        model.add(weighted_constraint <= ((double) global_problem_.budget_));
        weighted_constraint.end();
        
        IloExpr obj_exp(env);
        for(i = 0; i < nbImpressions; i++) {
            for(j = 0; j < nbAdvertisers; j++) {
                obj_exp += ((double) (*transpose_bids_matrix_)[i][j].second) * x[i][j];
            }
        }
        
        model.add(IloMaximize(env, obj_exp));
        obj_exp.end();
        
        IloCplex cplex(env);
        cplex.setOut(env.getNullStream());
        cplex.extract(model);
        
        // Optimize the problem and obtain solution.
        if ( !cplex.solve() ) {
            env.error() << "Failed to optimize LP" << endl;
            throw(-1);
        }
        
        IloNumArray vals(env);
        
        long double sum_b = 0;
        for (int a = 0; a < num_advertisers_; ++a) {
            //slacks_[a] = 0;
            for (i = 0; i < nbImpressions; i++) {
                if (cplex.getValue(x[i][a]) > 0) {
                    //slacks_[a] = slacks_[a] + cplex.getValue(x[i][a]) * bids_matrix_[a].find(i)->second;
                    sum_b += cplex.getValue(x[i][a]) * (*bids_matrix_)[a][i].second * weights_[a];
                }
            }
        }
        
        cout << "Cplex buget allocation is = ";
        cout << sum_b;
        cout << "\n";
        
        if (cplex.getObjValue() == CalculateGlobalMWProblemOpt()) {
            cout << "Solution checks \n";
        } else {
            cout << "Solution does not check, Cplex opt is ";
            cout << cplex.getObjValue();
            cout << "\n";
        }
        env.end();
    }
    
    long double AllocationMW::CalculateGlobalMWProblemOpt() {
        long double MW_revenue = 0;
        /*
         for (int a = 0; a < (*primal_sol_).size(); ++a) {
         for (__gnu_cxx::hash_map<int, long double>::iterator iter = (*primal_sol_)[a].begin();
         iter != (*primal_sol_)[a].end();
         ++iter) {
         MW_revenue += iter->second * (*bids_matrix_)[a][iter->first];
         }
         }
         */
        for (int a = 0; a < (*solution_).size(); ++a) {
            for (int i = 0; i < (*solution_)[a].size(); ++a) {
                MW_revenue += (*solution_)[a][i].second.first * (*bids_matrix_)[a][i].second;
                //iter->second.first * (*bids_matrix_)[a][iter->first];
            }
        }
        
        cout << "MW rev ";
        cout << MW_revenue;
        cout << "\n";
        return MW_revenue;
    }
    
    void AllocationMW::ReportWeightStats() {
        long double min_weight = 100000;
        long double max_weight = 0;
        for (int i = 0; i < num_advertisers_; ++i) {
            if (weights_[i] < min_weight) {
                min_weight = weights_[i];
            }
            if (weights_[i] > max_weight) {
                max_weight = weights_[i];
            }
        }
        cout << "min weight = ";
        cout << min_weight;
        cout << ", max weight = ";
        cout << max_weight;
        cout << "\n";
    }
    
    void AllocationMW::ComputeCPLEXRevenue() {
        IloEnv env;
        IloInt i, j;
        IloModel model(env);
        IloInt nbImpressions = num_impressions_;
        IloInt nbAdvertisers = num_advertisers_;
        
        NumVarMatrix x(env, nbImpressions);
        
        for(i = 0; i < nbImpressions; i++) {
            x[i] = IloNumVarArray(env, nbAdvertisers, 0.0, 1.0, ILOFLOAT);
        }
        
        // Add impression constraints.
        for(i = 0; i < nbImpressions; i++) {
            model.add(IloSum(x[i]) <= 1.0);
        }
        
        // Add weighted contraint.
        for(j = 0; j < nbAdvertisers; j++) {
            IloExpr curr_adv_constraint(env);
            for (int i = 0; i < (*bids_matrix_)[j].size(); ++i) {
                curr_adv_constraint += ((double) (*bids_matrix_)[j][i].second) *
                                       ((double) weights_[(int) (*bids_matrix_)[j][i].first]) *
                                       x[(int) (int) (*bids_matrix_)[j][i].first][j];
            }
            model.add(curr_adv_constraint <= ((double) (*budgets_)[j]));
        }
        
        IloExpr obj_exp(env);
        for(i = 0; i < nbImpressions; i++) {
            for(j = 0; j < nbAdvertisers; j++) {
                obj_exp += ((double) (*transpose_bids_matrix_)[i][j].second) * x[i][j];
            }
        }
        
        model.add(IloMaximize(env, obj_exp));
        obj_exp.end();
        
        IloCplex cplex(env);
        cplex.setOut(env.getNullStream());
        cplex.extract(model);
        
        // Optimize the problem and obtain solution.
        if ( !cplex.solve() ) {
            env.error() << "Failed to optimize LP" << endl;
            throw(-1);
        }
        
        //cplex.exportModel("/Users/ciocan/Documents/Google/data/example.lp");
        cout << "CPLEX opt is " << cplex.getObjValue() << "\n";
        env.end();
    }
    
    void AllocationMW::RunAllocationMW(int num_iterations) {
        cout << "Running allocation LP MW algorithm \n";
        for (int t = 1; t <= num_iterations; ++t) {
            cout << "Entering iteration " << t << "\n";
            
            // Get primal solution.
            clock_t t1, t2;
            float diff;
            t1 = clock();
            global_problem_.ConstructPrimal(t);
            t2 = clock();
            diff = ((float)t2-(float)t1);
            cout << "execution time of relaxation computation was " << diff << "\n";
            
            t1 = clock();
            //Instance::UpdateAvgPrimal(t, solution_);
            t2 = clock();
            diff = ((float)t2-(float)t1);
            cout << "execution time of average primal update was " << diff << "\n";
            
            // Runs CPLEX for debugging only, should be turned off.
            // VerifySolution();
            
            // Calculate slacks, update averages and recalculate weights.
            t1 = clock();
            CalculateSlacks();
            UpdateAvgSlacks(t);
            t2 = clock();
            diff = ((float)t2-(float)t1);
            cout << "execution time of slack update was " << diff << "\n";
            
            // Calculate slacks, update averages and recalculate weights.
            t1 = clock();
            ReportWorstInfeasibility(t);
            t2 = clock();
            diff = ((float)t2-(float)t1);
            cout << "execution time of worst infeasibility update was " << diff << "\n";
            
            // Calculate slacks, update averages and recalculate weights.
            t1 = clock();
            UpdateWeights();
            t2 = clock();
            diff = ((float)t2-(float)t1);
            cout << "execution time of weight update was " << diff << "\n";
            
            t1 = clock();
            UpdateGlobalProblem();
            t2 = clock();
            diff = ((float)t2-(float)t1);
            cout << "execution time of global problem update was " << diff << "\n";
            
            ReportWeightStats();
        }
    }
}
