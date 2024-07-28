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

    apg::bitword origin;
    apg::bitword reaction_allowed; // includes both reaction and possible objects placement
    apg::bitword objects_forbidden; // objects cannot be placed here as it would react with pattern history

    uint32_t max_gliders_allowed = 2;
    uint32_t origin_period;
    uint32_t max_branching = 16384; //to experiment with
    uint32_t max_extra_gens = 1024;

    uint32_t best_solution_cost = 999999999;
    std::vector<ProblemState> solutions; // differ by gliders issued, minimizing object cost for each possibility

    UniqueEnsurerTable hash_to_added_object_cost;

    uint32_t get_best_solution_cost() {
        std::unique_lock<mutex_type> lock(mtx);
        return best_cost_solution;
    }

    uint32_t _get_cost_unsafe(uint64_t hash) const {
        uint32_t idx = hash_to_added_object_cost.find(hash);
        if (idx == ((uint32_t) -1)) {
            return ((uint32_t) -1);
        } else {
            return hash_to_added_object_cost[idx].value;
        }
    }

    bool _leq_cost_unsafe(uint64_t hash, uint32_t cost) const {
        if (best_solution_cost <= cost) { return true; }
        return (_get_cost_unsafe(hash) <= cost);
    }

    uint64_t size() {
        std::unique_lock<mutex_type> lock(mtx);
        return hash_to_added_object_cost.size();
    }

    bool leq_cost(uint64_t hash, uint32_t cost) {
        std::unique_lock<mutex_type> lock(mtx);
        return _leq_cost_unsafe(hash, cost);
    }

    bool leq_cost_update(uint64_t hash, uint32_t cost, bool is_clean_solution) {
        std::unique_lock<mutex_type> lock(mtx);
        if (_leq_cost_unsafe(hash, cost)) {
            return true;
        } else {
            if (is_clean_solution) {
                std::cerr  << std::endl << "\033[32;1mFound clean solution with " << cost << " of added_objects .\033[0m";
                best_solution_cost = cost;
            }
            UniqueEnsurerEntry uee;
            uee.key = hash;
            uee.value = cost;
            uee.next = 0;
            bool overwrite = true;
            hash_to_added_object_cost.insert(uee, overwrite);
            return false;
        }
    }

    void save_last_solution() {
        std::ofstream out("SoD_"+std::to_string(solutions.size())+".mc");
        solutions[solutions.size()-1].starting_pattern(origin).write_macrocell(out);
        out.close();
    }

    void save_solution(ProblemState ps) {
        std::unique_lock<mutex_type> lock(mtx);
        bool has_no_gliders = ps.gliders.totalPopulation() == 0;
        if (has_no_gliders) {
            if (ps.added_objects_cost == best_solution_cost) {
                std::ofstream out("SoD_best_cleean.mc");
                ps.starting_pattern(origin).write_macrocell(out);
                out.close();
            } else {
                std::cerr << "\33[31;1mClean solution not best (" << added_objects_cost << "!=" << best_solution_cost << ") current best added objects cost!\033[0m" << std::endl;
            }
        }
        solutions.push_back(ps); save_last_solution();
    }

    void save_progress(ProblemState ps, uint64_t depth, uint64_t beamIndex) {
        std::ofstream out("SoD_Progress_"+std::to_string(depth)+"_"+std::to_string(beamIndex)+".mc");
        ps.starting_pattern(origin).write_macrocell(out);
        out.close();
    }

};
} // namespace gsod
