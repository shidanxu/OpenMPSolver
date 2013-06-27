//
//  instance.cpp
//  LPsolver
//
//  Created by Dragos Ciocan on 6/20/13.
//  Copyright (c) 2013 Dragos Ciocan. All rights reserved.
//


#include <cmath>
#include "convex_hull.h"
#include "instance.h"

namespace distributed_solver {

    Instance::Instance(int num_advertisers, int num_impressions, int num_slots, long double bid_sparsity,
                       long double epsilon, long double scaling_factor, long double numerical_accuracy_tolerance) {
        epsilon_ = epsilon;
        scaling_factor_ = scaling_factor;
        iteration_count_ = 0;
        numerical_accuracy_tolerance_ = numerical_accuracy_tolerance;

        num_advertisers_ = num_advertisers;
        num_impressions_ = num_impressions;
        num_slots_ = num_slots;
        bid_sparsity_ = bid_sparsity;
        num_shards_ = 10;

        budgets_ = vector<long double>(num_advertisers_);
        SetBudgets();
    }

    void Instance::GenerateInstance() {
        vector<__gnu_cxx::hash_map<int, long double> >* bid_mat_hm = new vector<__gnu_cxx::hash_map<int, long double> >();
        vector<__gnu_cxx::hash_map<int, long double> >* transpose_bid_mat_hm = new vector<__gnu_cxx::hash_map<int, long double> >();
        transpose_bids_matrix_.reserve(num_impressions_);
        bids_matrix_.reserve(num_advertisers_);
        for (int i = 0; i < num_impressions_; ++i) {
            transpose_bid_mat_hm->push_back(*new __gnu_cxx::hash_map<int, long double>());
        }
        srand(1);
        max_bid_ = 0;
        for (int j = 0; j < num_advertisers_; ++j) {
            // Generate bids for advertiser j.
            __gnu_cxx::hash_map<int, long double> bid_row;
            for (int i = 0; i < (bid_sparsity_ * num_impressions_); ++i) {
                int index = rand() % num_impressions_;
                long double bid = (long double) (rand() + 1) / ((long double) RAND_MAX);
                if (max_bid_ < bid) {
                    max_bid_ = bid;
                }
                bid_row[index] = bid;
                (*transpose_bid_mat_hm)[index][j] = bid;
            }
            bid_mat_hm->push_back(bid_row);
        }

        for (int i = 0; i < num_impressions_; ++i) {
            transpose_bids_matrix_.push_back(*new vector<pair<int, long double> >());
            for (__gnu_cxx::hash_map<int, long double>::const_iterator iter = (*transpose_bid_mat_hm)[i].begin();
                 iter != (*transpose_bid_mat_hm)[i].end(); iter++) {
                //cout << iter->first << "  " << iter->second << "\n";
                transpose_bids_matrix_[i].push_back(make_pair(iter->first, iter->second));
            }
        }

        for (int a = 0; a < num_advertisers_; ++a) {
            bids_matrix_.push_back(*new vector<pair<int, long double> >());
            for (__gnu_cxx::hash_map<int, long double>::const_iterator iter = (*bid_mat_hm)[a].begin();
                 iter != (*bid_mat_hm)[a].end(); iter++) {
                bids_matrix_[a].push_back(make_pair(iter->first, iter->second));
            }
        }

        delete bid_mat_hm;
        delete transpose_bid_mat_hm;
        cout << "Generated instance \n";
        // ReportGraphTopology();
    }

    void Instance::RunMultiplicativeWeights(long double num_iterations, long double numerical_accuracy_tolerance) {
        BuildPrimals();
        AllocationMW alloc_mw = AllocationMW(num_advertisers_, num_impressions_, num_slots_,
                                             bid_sparsity_, max_bid_, epsilon_, numerical_accuracy_tolerance_,
                                             &bids_matrix_, &transpose_bids_matrix_, &budgets_, solution_);
        alloc_mw.RunAllocationMW(num_iterations);
    }

    void Instance::SetBudgets() {
        long double average_bid = 0.5; // Need to change this manually depending on how bids are drawn.
        for (int j = 0; j < num_advertisers_; ++j) {
            budgets_[j] = average_bid * (num_impressions_ / num_advertisers_) * scaling_factor_;
        }
    }

    void Instance::UpdatePrimal(int t,
                                vector<vector<pair<int, pair<long double, long double> > > >* solution,
                                const vector<pair<pair<int, long double>, pair<int, long double> > >& primal_changes) {
        __gnu_cxx::hash_map<int, long double>::const_iterator iter;
        for (int a = 0; a < (*solution).size(); ++a) {
            for (int j = 0; j < (*solution)[a].size(); ++j) {
                //iter = (*primal_changes)[a].find((*solution)[a][j].first);
                //if (iter != (*primal_changes)[a].end())  {
                    //cout << "setting " << j << ", " << a << " to " << iter->second << "\n";
                //    (*solution)[a][j].second.first = iter->second;
                //}
                if (primal_changes[j].first.first == a) {
                    (*solution)[a][j].second.first = primal_changes[j].first.second;
                } else if (primal_changes[j].second.first == a) {
                    (*solution)[a][j].second.first = primal_changes[j].second.second;
                }
                (*solution)[a][j].second.second = (long double) (t - 1) / t * (*solution)[a][j].second.second +
                                                  (long double) 1 / t * (*solution)[a][j].second.first;
            }
        }
    }

    void Instance::BuildPrimals() {
        solution_ = new vector<vector<pair<int, pair<long double, long double> > > >();
        for (int a = 0; a < num_advertisers_; ++a) {
            vector<pair<int, pair<long double, long double> > > row;
            for (int j = 0; j < bids_matrix_[a].size(); ++j) {
                row.push_back(make_pair(bids_matrix_[a][j].first, make_pair(0.0, 0.0)));
            }
            solution_->push_back(row);
        }
    }

    void Instance::ResetCurrentPrimal(vector<vector<pair<int, pair<long double, long double> > > >* sol) {
        for (int a = 0; a < (*sol).size(); ++a) {
            for (int j = 0; j < (*sol)[a].size(); ++j) {
                (*sol)[a][j].second.first = 0.0;
            }
        }
    }

    void Instance::ReportGraphTopology() {
        for (int a = 0; a < bids_matrix_.size(); ++a) {
            cout << "Advertiser " << a << " degree is " << bids_matrix_[a].size() << "\n";
        }
        for (int i = 0; i < transpose_bids_matrix_.size(); ++i) {
            cout << "Impression " << i << " degree is " << transpose_bids_matrix_[i].size() << "\n";
        }
    }
}
