#pragma once
#include <stdint.h>
#include <type_traits>
#include <vector>
#include <set>
#include <unordered_map>
#include <map>
#include <sstream>
#include <string>

//#include "lifelib/cpads/include/cpads/mxor.hpp"
//#include "unique_ensurer.hpp"
//#include "lifelib/avxlife/uli.h"

namespace gsod {

struct ObjectWithCost {
    apg::pattern object, object_env, object_env_flipped;
    uint32_t object_cost;
};

struct ProblemState {

    apg::bitworld added_objects; // including starting glider
    uint32_t added_objects_cost; // from added_objects_cost ... evaluation of result

    // fields to be computed from above and origin (reaction_allowed, objects_forbidden)
    uint64_t early_hash;
    uint32_t early_generation, locally_stable_generation, stable_generation;
    uint32_t num_output_gliders;
    //apg::bitworld output_gliders; not needed during the computation could be calculated when processing a(n unclear) solution
    double spanning_tree_cost; // evaluation of progress
    double total_cost(double pessimism) const {
        return spanning_tree_cost * pessimism + added_objects_cost; // todo think about the relative constant
    }

};

struct BeamSearchContainer {

    uint32_t maxsize;
    std::vector<ProblemState> contents;
    std::set<std::pair<double, uint32_t>> pmpq; // poor man's priority queue
    std::unordered_map<uint64_t, uint32_t> hash_to_idx;

    void try_insert(const ProblemState& ps, double pessimism) {

        auto it = hash_to_idx.find(ps.early_hash);
        if (it != hash_to_idx.end()) {//cheaper way to reach the same position
            uint32_t idx = it->second;
            if (ps.added_objects_cost < contents[idx].added_objects_cost) {
                pmpq.erase(std::pair<double, uint32_t>{contents[idx].total_cost(pessimism), idx});
                contents[idx] = ps;
                pmpq.insert(std::pair<double, uint32_t>{contents[idx].total_cost(pessimism), idx});
            }
            return;
        }

        uint32_t idx = contents.size();

        if (idx >= maxsize) {
            auto opair = *(pmpq.rbegin());
            if (ps.total_cost(pessimism) >= opair.first) { return; }
            idx = opair.second;
            pmpq.erase(opair);
            hash_to_idx.erase(contents[idx].early_hash);
            contents[idx] = ps;
        } else {
            contents.push_back(ps);
        }

        pmpq.insert(std::pair<double, uint32_t>{contents[idx].total_cost(pessimism), idx});
        hash_to_idx[ps.early_hash] = idx;
    }

    void add(const BeamSearchContainer &other, double pessimism) {
        for (const auto &x : other.contents) {
            try_insert(x, pessimism);
        }
    }
};

} // namespace gsod
