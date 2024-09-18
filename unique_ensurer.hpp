#pragma once
#include <mutex>
#include "lifelib/cpads/include/cpads/indextable.hpp"
#include "lifelib/pattern2.h"

namespace gsod {

/**
 * We use hh::indextable instead of std::unordered_map because it allows us
 * to store up to 2**32 mappings from 8-byte keys to 4-byte values with only
 * 24-32 bytes of memory per element. (By comparison, std::unordered_map
 * takes around 48 bytes of memory per element.)
 */
struct UniqueEnsurerEntry {
    uint64_t key;
    uint32_t next;
    uint32_t value;
};

static_assert(sizeof(UniqueEnsurerEntry) == 16, "UniqueEnsurerEntry should be 16 bytes");

struct UniqueEnsurerTable : public hh::indextable<UniqueEnsurerEntry, uint64_t, UniqueEnsurerTable> {
    UniqueEnsurerTable() : hh::indextable<UniqueEnsurerEntry, uint64_t, UniqueEnsurerTable>(1) { }
    auto compute_key(const UniqueEnsurerEntry &element) const { return element.key; }
};


/**
 * Ensures that we never duplicate or revisit a state. We do this by storing
 * a global append-only mapping (64-bit hash --> 32-bit glider count) which
 * is updated every time we produce a child state, even if that child state
 * is outside the width of the beam search.
 */
struct UniqueEnsurer {

    using mutex_type = std::mutex;

    mutex_type mtx;

    apg::bitworld origin;
    apg::bitworld reaction_allowed; // includes both reaction and possible objects placement
    apg::bitworld objects_forbidden; // objects cannot be placed here as it would react with pattern history
    apg::bitworld input_glider_lines_forbidden; // pixel in a lane blocks incoming glider in the lane (blocks just one x+y lane and/or one y-x lane) outside reaction_allowed
    apg::bitworld ori_added_objects;
    uint32_t max_output_gliders_allowed = 2;
    uint32_t max_object_types_allowed = 10;
    uint32_t origin_period;
    uint32_t max_branching = 16384; //to experiment with
    uint32_t max_extra_gens = 1024;
    uint32_t pessimism = 70;

    uint32_t best_solution_cost = 999999999;
    std::vector<ProblemState> solutions; // differ by gliders issued, minimizing object cost for each possibility

    UniqueEnsurerTable hash_to_added_object_cost;

    UniqueEnsurer (const UniqueEnsurer &_ue);

    UniqueEnsurer (apg::bitworld _origin, apg::bitworld _reaction_allowed, apg::bitworld _objects_forbidden, apg::bitworld _input_glider_lines_forbidden, apg::bitworld _ori_added_objects) {
        origin = _origin; reaction_allowed = _reaction_allowed; objects_forbidden = _objects_forbidden; input_glider_lines_forbidden = _input_glider_lines_forbidden; ori_added_objects = _ori_added_objects;
    };

    uint32_t get_best_solution_cost() {
        std::unique_lock<mutex_type> lock(mtx);
        return best_solution_cost;
    }

    uint32_t _get_cost_unsafe(uint64_t hash) const {
        uint32_t idx = hash_to_added_object_cost.find(hash);
        if (idx == ((uint32_t) -1)) {
            //std::cerr << "m";
            return ((uint32_t) -1);
        } else {
            //std::cerr << "v(" << idx << ":" << hash_to_added_object_cost[idx].value << ")" ;
            return hash_to_added_object_cost[idx].value;
        }
    }

    bool _leq_cost_unsafe(uint64_t hash, uint32_t cost) const {
        if (best_solution_cost <= cost) {
            //std::cerr << "b";
            return true;
        }
        return (_get_cost_unsafe(hash) <= cost);
    }

    uint32_t size() {
        std::unique_lock<mutex_type> lock(mtx);
        return hash_to_added_object_cost.size();
    }

    bool leq_cost(uint64_t hash, uint32_t cost) {
        std::unique_lock<mutex_type> lock(mtx);
        return _leq_cost_unsafe(hash, cost);
    }

    /*
        negative solution_clarity means there remains something to destroy
        positive is cost of output gliders (each is supposed to be stopped by a blinker (cheapest object)
    */
    bool leq_cost_update(uint64_t hash, uint32_t cost) {
        std::unique_lock<mutex_type> lock(mtx);
        if (_leq_cost_unsafe(hash, cost)) {
            //std::cerr << "h(" << hash << ")";
            return true;
        } else {
            UniqueEnsurerEntry uee;
            uee.key = hash;
            uee.value = cost;
            uee.next = 0;
            bool overwrite = true;
            hash_to_added_object_cost.insert(uee, overwrite);
            //std::cerr << "n";
            return false;
        }
    }

    bool best_cost_update(uint32_t cost, int solution_clarity) {
        std::unique_lock<mutex_type> lock(mtx);
        uint32_t best_cost_estimate = cost + std::abs(solution_clarity);
        if (best_cost_estimate < best_solution_cost) {
            if (solution_clarity >= 0) {
                std::cerr  << std::endl << "\033[32;1mFound solution with cost " << cost << " of added_objects with possble cost " << solution_clarity << " to stop output gliders.\033[0m";
            } else {
                std::cerr  << std::endl << "\033[32;1mSolution cost estimate changed to " << best_cost_estimate << " according to partial " << cost << " of added_objects with estimated cost " << solution_clarity << " to continue with gliders.\033[0m";
            }
            best_solution_cost = best_cost_estimate;
            return true;
        }
        return false;
    }

    void save_last_solution(apg::pattern start, uint32_t cost, uint32_t num_output_gliders, uint32_t expected_extra_cost) {
        std::ofstream out("gSoD_"+std::to_string(solutions.size())+"_"+std::to_string(cost)+"_"+std::to_string(num_output_gliders)+"_"+std::to_string(expected_extra_cost)+".mc");
        start.write_macrocell(out);
        out.close();
    }

    void save_solution(ProblemState ps, apg::pattern start, int solution_clarity) {
        std::unique_lock<mutex_type> lock(mtx);
        if (solution_clarity == 0) {
            if (ps.added_objects_cost == best_solution_cost) {
                std::ofstream out("SoD_best_clean.mc");
                start.write_macrocell(out);
                out.close();
            } else {
                std::cerr << "\33[31;1mClean solution not best (" << ps.added_objects_cost << "!=" << best_solution_cost << ") current best added objects cost!\033[0m" << std::endl;
            }
        }
        solutions.push_back(ps); save_last_solution(start, ps.added_objects_cost, ps.num_output_gliders, (solution_clarity<0) ? -solution_clarity : 0);
    }

    void save_progress(apg::pattern start, uint32_t depth, uint32_t beamIndex, double spanning_tree_cost, uint32_t added_objects_cost) {
        std::ofstream out("SoD_Progress_"+std::to_string(depth)+"_"+std::to_string(beamIndex)+"_"+std::to_string(spanning_tree_cost)+"_"+std::to_string(added_objects_cost)+".mc");
        start.write_macrocell(out);
        out.close();
    }

};
} // namespace gsod
