#include <unistd.h>
#include <map>
#include <thread>
#include <atomic>
#include <vector>

//#include "lifelib/cpads/include/cpads/random/prng.hpp"
#include "lifelib/spantree.h"
#include "lifelib/bitworld.h"
#include "lifelib/lifetree.h"
#include "lifelib/pattern2.h"
#include "structures.hpp"
#include "unique_ensurer.hpp"

#define BIG_SAVE_FREQ 1 //0

namespace gsod {

template<int N>
std::vector<apg::bitworld> load_file(std::string filename) {
    std::vector<apg::bitworld> vbw;
    apg::lifetree<uint32_t, N> lt(100);
    apg::pattern pat = apg::pattern(&lt, filename);
    for (int i = 0; i < N; i++) { vbw.push_back(pat.flatlayer(i)); }
    return vbw;
}

void save_pattern(apg::pattern patt, std::string filename) {
    std::ofstream out(filename);
    patt.write_macrocell(out);
    out.close();
}

static ObjectWithCost make_ObjectWithCost (apg::pattern _object, uint32_t _object_cost) {
        apg::pattern _object_env = _object + _object[1];
        apg::pattern _object_env_flipped = _object_env.transform("flip", 0, 0);
        return {_object, _object_env, _object_env_flipped, _object_cost};
}

//double plus40(double sqlen) { // square of the euklid distance looks like the best linear estimate of the final cost (*0.7)
//    return sqlen;
//}

struct SODThread {

    apg::lifetree<uint32_t, 1> lt;

    apg::pattern origin; // life pattern to destroy
    apg::pattern reaction_allowed; // includes both reaction and possible objects placement
    apg::pattern objects_forbidden; // objects cannot be placed here as it would react with pattern history
    //I should allow to restrict some input glider lanes+directions.
    apg::pattern influence; //handy constant
    std::vector<ObjectWithCost> objectsSource; //all p at most 2

    //int power_thousands; used to search optimal power of euklid dist for edges

    SODThread() : lt(1000), origin(&lt, "", "b3s23"), reaction_allowed(&lt, "", "b3s23"), objects_forbidden(&lt, "", "b3s23"), influence(&lt, "b3o$5o$5o$5o$b3o!", "b3s23") {
        influence = influence.shift(-2, -2);
    }

    apg::pattern b2p(const apg::bitworld &bw) {
        std::vector<apg::bitworld> vbw(1, bw);
        apg::pattern pat(&lt, vbw, "b3s23");
        return pat;
    }

    void init_patterns(UniqueEnsurer &ue) {//each thread have different lt, so must have its separate pattern "constants"
        origin = b2p(ue.origin); reaction_allowed = b2p(ue.reaction_allowed); objects_forbidden = b2p(ue.objects_forbidden);

        if (ue.max_object_types_allowed>0) {
            apg::pattern blinker(&lt, "3o!", "b3s23");
                objectsSource.push_back(make_ObjectWithCost(blinker, 467));
                objectsSource.push_back(make_ObjectWithCost(blinker[1], 467));
        }
        if (ue.max_object_types_allowed>1) {
            apg::pattern block(&lt, "2o$2o!", "b3s23");
                objectsSource.push_back(make_ObjectWithCost(block, 493));
        }
        if (ue.max_object_types_allowed>2) {
            apg::pattern beehive(&lt, "b2o$o2bo$b2o!", "b3s23");
                objectsSource.push_back(make_ObjectWithCost(beehive, 615));
                objectsSource.push_back(make_ObjectWithCost(beehive.transform("rcw", 0, 0), 615));
        }
        if (ue.max_object_types_allowed>3) {
            apg::pattern pond(&lt, "b2o$o2bo$o2bo$b2o!", "b3s23");
                objectsSource.push_back(make_ObjectWithCost(pond, 685));
        }
        if (ue.max_object_types_allowed>4) {
            apg::pattern loaf(&lt, "2bo$bobo$o2bo$b2o!", "b3s23");
                objectsSource.push_back(make_ObjectWithCost(loaf, 714));
                objectsSource.push_back(make_ObjectWithCost(loaf.transform("rcw", 0, 0), 714));
                objectsSource.push_back(make_ObjectWithCost(loaf.transform("flip", 0, 0), 714));
                objectsSource.push_back(make_ObjectWithCost(loaf.transform("rccw", 0, 0), 714));
        }
        if (ue.max_object_types_allowed>5) {
            apg::pattern tub(&lt, "bo$obo$bo!", "b3s23");
                objectsSource.push_back(make_ObjectWithCost(tub, 726));
        }
        if (ue.max_object_types_allowed>6) {
            apg::pattern boat(&lt, "2o$obo$bo!", "b3s23");
                objectsSource.push_back(make_ObjectWithCost(boat, 743));
                objectsSource.push_back(make_ObjectWithCost(boat.transform("rcw", 0, 0), 743));
                objectsSource.push_back(make_ObjectWithCost(boat.transform("flip", 0, 0), 743));
                objectsSource.push_back(make_ObjectWithCost(boat.transform("rccw", 0, 0), 743));
        }
        if (ue.max_object_types_allowed>7) {
            apg::pattern ship(&lt, "2o$obo$b2o!", "b3s23");
                objectsSource.push_back(make_ObjectWithCost(ship, 767));
                objectsSource.push_back(make_ObjectWithCost(ship.transform("rcw", 0, 0), 767));
        }
        if (ue.max_object_types_allowed>8) {
            apg::pattern barge(&lt, "bo$obo$bobo$2bo!", "b3s23");
                objectsSource.push_back(make_ObjectWithCost(barge, 816));
                objectsSource.push_back(make_ObjectWithCost(barge.transform("rcw", 0, 0), 816));
        }
        if (ue.max_object_types_allowed>9) {
            apg::pattern longboat(&lt, "2o$obo$bobo$2bo!", "b3s23");
                objectsSource.push_back(make_ObjectWithCost(longboat, 868));
                objectsSource.push_back(make_ObjectWithCost(longboat.transform("rcw", 0, 0), 868));
                objectsSource.push_back(make_ObjectWithCost(longboat.transform("flip", 0, 0), 868));
                objectsSource.push_back(make_ObjectWithCost(longboat.transform("rccw", 0, 0), 868));
        }
        if (ue.max_object_types_allowed>10) {
            apg::pattern mango(&lt, "b2o$o2bo$bo2bo$2b2o!", "b3s23");
                objectsSource.push_back(make_ObjectWithCost(mango, 1010));
                objectsSource.push_back(make_ObjectWithCost(mango.transform("rcw", 0, 0), 1010));
        }
        if (ue.max_object_types_allowed>11) {
            apg::pattern shiptie(&lt, "2o$obo$b2o$3b2o$3bobo$4b2o!", "b3s23");
                objectsSource.push_back(make_ObjectWithCost(shiptie, 1050));
                objectsSource.push_back(make_ObjectWithCost(shiptie.transform("rcw", 0, 0), 1050));
        }
        if (ue.max_object_types_allowed>12) {
            apg::pattern eater1(&lt, "2o$obo$2bo$2b2o!", "b3s23");
                objectsSource.push_back(make_ObjectWithCost(eater1, 1120));
                objectsSource.push_back(make_ObjectWithCost(eater1.transform("rcw", 0, 0), 1120));
                objectsSource.push_back(make_ObjectWithCost(eater1.transform("flip", 0, 0), 1120));
                objectsSource.push_back(make_ObjectWithCost(eater1.transform("rccw", 0, 0), 1120));
                objectsSource.push_back(make_ObjectWithCost(eater1.transform("swap_xy", 0, 0), 1120));
                objectsSource.push_back(make_ObjectWithCost(eater1.transform("flip_x", 0, 0), 1120));
                objectsSource.push_back(make_ObjectWithCost(eater1.transform("flip_y", 0, 0), 1120));
                objectsSource.push_back(make_ObjectWithCost(eater1.transform("swap_xy_flip", 0, 0), 1120));
        }
        if (ue.max_object_types_allowed>13) {
            apg::pattern aircraftcarrier(&lt, "2o$o$2bo$b2o!", "b3s23");
                objectsSource.push_back(make_ObjectWithCost(aircraftcarrier, 1370));
                objectsSource.push_back(make_ObjectWithCost(aircraftcarrier.transform("rcw", 0, 0), 1370));
                objectsSource.push_back(make_ObjectWithCost(aircraftcarrier.transform("flip_x", 0, 0), 1370));
                objectsSource.push_back(make_ObjectWithCost(aircraftcarrier.transform("swap_xy", 0, 0), 1370));
        }
        if (ue.max_object_types_allowed>14) {
            apg::pattern elevener(&lt, "2o$obo$2bo$2b3o$5bo$4b2o!", "b3s23");
                objectsSource.push_back(make_ObjectWithCost(elevener, 1430));
                objectsSource.push_back(make_ObjectWithCost(elevener.transform("rcw", 0, 0), 1430));
                objectsSource.push_back(make_ObjectWithCost(elevener.transform("flip", 0, 0), 1430));
                objectsSource.push_back(make_ObjectWithCost(elevener.transform("rccw", 0, 0), 1430));
        }
        if (ue.max_object_types_allowed>15) {
            apg::pattern intergralsign(&lt, "2o$obo$2bo$2bobo$3b2o!", "b3s23");
                objectsSource.push_back(make_ObjectWithCost(intergralsign, 1450));
                objectsSource.push_back(make_ObjectWithCost(intergralsign.transform("rcw", 0, 0), 1450));
                objectsSource.push_back(make_ObjectWithCost(intergralsign.transform("flip_x", 0, 0), 1450));
                objectsSource.push_back(make_ObjectWithCost(intergralsign.transform("swap_xy", 0, 0), 1450));
        }
        if (ue.max_object_types_allowed>16) {
            apg::pattern shillelagh(&lt, "2o$obo$2bo$bo$b2o!", "b3s23");
                objectsSource.push_back(make_ObjectWithCost(shillelagh, 1470));
                objectsSource.push_back(make_ObjectWithCost(shillelagh.transform("rcw", 0, 0), 1470));
                objectsSource.push_back(make_ObjectWithCost(shillelagh.transform("flip", 0, 0), 1470));
                objectsSource.push_back(make_ObjectWithCost(shillelagh.transform("rccw", 0, 0), 1470));
                objectsSource.push_back(make_ObjectWithCost(shillelagh.transform("swap_xy", 0, 0), 1470));
                objectsSource.push_back(make_ObjectWithCost(shillelagh.transform("flip_x", 0, 0), 1470));
                objectsSource.push_back(make_ObjectWithCost(shillelagh.transform("flip_y", 0, 0), 1470));
                objectsSource.push_back(make_ObjectWithCost(shillelagh.transform("swap_xy_flip", 0, 0), 1470));
        }
        if (ue.max_object_types_allowed>17) {
            apg::pattern transblockonlonghook(&lt, "2ob2o$2obo$3bo$3bobo$4b2o!", "b3s23");
                objectsSource.push_back(make_ObjectWithCost(transblockonlonghook, 1480));
                objectsSource.push_back(make_ObjectWithCost(transblockonlonghook.transform("rcw", 0, 0), 1480));
                objectsSource.push_back(make_ObjectWithCost(transblockonlonghook.transform("flip", 0, 0), 1480));
                objectsSource.push_back(make_ObjectWithCost(transblockonlonghook.transform("rccw", 0, 0), 1480));
                objectsSource.push_back(make_ObjectWithCost(transblockonlonghook.transform("swap_xy", 0, 0), 1480));
                objectsSource.push_back(make_ObjectWithCost(transblockonlonghook.transform("flip_x", 0, 0), 1480));
                objectsSource.push_back(make_ObjectWithCost(transblockonlonghook.transform("flip_y", 0, 0), 1480));
                objectsSource.push_back(make_ObjectWithCost(transblockonlonghook.transform("swap_xy_flip", 0, 0), 1480));
        }
    }

    void early_and_late(apg::pattern x, uint32_t early_gens, apg::pattern &early_pattern, apg::pattern &early_env, uint32_t late_gens, apg::pattern &late_env) {

        apg::pattern y = x;
        apg::pattern z = x;

        for (uint32_t i = 0; i < late_gens; i++) {
            if (i==early_gens) {
                early_pattern=z;
                early_env=y;
                y=z;
            }
            z = z[1];
            y += z;
        }

        late_env = y;
    }

    // parsing lifehistory pattern (near (0,0)!) with history of the circuit.
    // It can contain state 4 extension marking additional reaction_allowed
    // otherwise a default neighbourhood of the pattern is chosen (not implemented yet)
    // states 5 and 1 have the same behaviour
    // state 3 outside bounding box of states 1,5 denotes forbidden directed input glider lines
    // states 1,2,5 close neighbourhood prevents added objects
    UniqueEnsurer parse_pattern(const std::vector<apg::bitworld> &vbw) {
        for(uint32_t l=0;l<vbw.size();l++) {// debug
            save_pattern(b2p(vbw[l]), "Input_plane_"+std::to_string(l)+".mc");
        }
        apg::pattern ori_added_objects_patt = b2p(vbw[2]); // state 5
        apg::pattern origin_patt = (b2p(vbw[0]) - b2p(vbw[3])); //state 1
        apg::pattern input_glider_lines_forbidden_patt = b2p(vbw[0]) - origin_patt - ori_added_objects_patt; // state 3
        apg::pattern reaction_allowed_patt = b2p(vbw[1])-input_glider_lines_forbidden_patt; //states 1,2,4,5
        apg::pattern state4 = b2p(vbw[3])-b2p(vbw[0]); // state 4
        apg::pattern objects_forbidden_patt = (reaction_allowed_patt - state4).convolve(influence);// historyenvelope (1,2,5) close neighbourhood
        if (state4.empty()) {// states 4 not present let us fill them
            //todo reaction_allowed_patt = approx_convex_hull(reaction_allowed_patt);
            apg::pattern halo2 = influence;
            halo2 = halo2.convolve(halo2).convolve(halo2); // radius 6
            halo2 = halo2.convolve(halo2); // radius 12
            halo2 = halo2.convolve(halo2); // radius 24
            reaction_allowed_patt = reaction_allowed_patt.convolve(halo2);
        }
        save_pattern(origin_patt, "gSoD_origin.mc");
        save_pattern(ori_added_objects_patt, "ori_added_objects.mc");
        save_pattern(input_glider_lines_forbidden_patt, "gSoD_input_glider_lines_forbidden.mc");
        save_pattern(reaction_allowed_patt, "gSoD_reaction_allowed.mc");
        save_pattern(objects_forbidden_patt, "gSoD_objects_forbidden.mc");

        UniqueEnsurer ue(origin_patt.flatlayer(0),reaction_allowed_patt.flatlayer(0),objects_forbidden_patt.flatlayer(0),input_glider_lines_forbidden_patt.flatlayer(0),ori_added_objects_patt.flatlayer(0));
        ue.solutions.clear();
        ue.origin_period = origin_patt.ascertain_period();
        if (ue.origin_period == 0) {
            std::cerr << "Pattern does not have short enough period!" << std::endl;
        }
        return ue;
    }

    double spanning_tree_cost_calc(apg::pattern envelope) {
        apg::pattern remainder = envelope;
        std::vector<apg::coords64> ccentres;
        int64_t bbox[4] = {0};

        while(remainder.nonempty()) {
            apg::pattern next = remainder.onecell(), curr=influence;
            while(!(curr==next)) {
                curr = next; next = curr.convolve(influence) & remainder;
            }
            curr.getrect(bbox);
            ccentres.emplace_back(bbox[0] * 2 + bbox[2], bbox[1] * 2 + bbox[3]);
            remainder -= curr;
        }

        std::vector<apg::edge64> all_edges = apg::spanning_graph(ccentres);
        /*if (power_thousands <= 0) {
            std::cerr << "Table of STT values depending on power_thousands starting on 0:" << std::endl;
            for(power_thousands = 0; power_thousands<1000; ++power_thousands) {
                double pt_stlength = apg::spanning_graph_to_tree(ccentres, 0, identity , all_edges); // actually constant 2.0 scale does not matter
                std::cerr << pt_stlength << std::endl;
            }
        }
        power_thousands = 500; // seems 1000 so identity function would work best ...
        */

        double stlength = apg::spanning_graph_to_tree(ccentres, 0, apg::pow5over8, all_edges); // actually constant 2.0 scale does not matter
        // if one wants to follow simeks program std::sqrt should be replaced by its 5/4 power (or better calculate the tree and then recompute the cost using this evaluation).
        // I would first experiment with the simple euklid distance mst.
        return stlength;
    }

    uint32_t get_stable_generation(apg::pattern starting_pattern, uint32_t max_gen, apg::pattern &out_gliders, apg::pattern &stabilised, UniqueEnsurer &ue) {
        uint32_t base_step = (ue.origin_period>8) ? ue.origin_period : 8;
        apg::pattern curr = starting_pattern, curr_base = curr & reaction_allowed;
        bool started = false;
        for(uint32_t gen=0;gen<max_gen;gen+=base_step) {
            apg::pattern next = curr[base_step], next_base = next & reaction_allowed;
            if (next_base == curr_base) {
                if (started) {
                    stabilised = curr_base;
                    //if (!((next & curr) - next_base).empty()) {return 0;} ... it could miss a case with gliders intersecting, but it would avoid stable trash of correct size
                    return gen;
                }
            } else {
                started = true;
            }
            out_gliders = (next - next_base).convolve(influence) & next; // part of the gliders could still be in allowed
            uint32_t gliders_pop = out_gliders.totalPopulation();
            if (((gliders_pop % 5)!=0) || (gliders_pop / 5 > ue.max_output_gliders_allowed)) {
                return 0; // forbidden
                // too many gliders or gliders not formed cleanly before leaving allowed
                // there could still be unprobable case of objects of propper size created out of allowed with propper sizes in each checkpoint on the way to fake extra solutions
                // I would risk it (could be made safer later)
            }
            curr = next; curr_base = next_base;
        }
        return 0;
    }

    uint32_t get_locally_stable_generation(apg::pattern starting_pattern, uint32_t stable_gen, apg::pattern stabilised, UniqueEnsurer &ue) {
        uint32_t base_step = (ue.origin_period>8) ? ue.origin_period : 8;
        apg::pattern stable_influenced = stabilised.convolve(influence);
        apg::pattern curr = starting_pattern, curr_stabble = curr;
        uint32_t last_locally_non_stable_gen = 999999999;
        for(uint32_t gen=0;gen < stable_gen;gen+=base_step) {
            curr = curr[base_step];
            curr_stabble = curr & stable_influenced;
            if (curr_stabble !=  stabilised) {
                last_locally_non_stable_gen = gen;
            }
        }
        if (last_locally_non_stable_gen == 999999999) {
            return 0;
        }
        return last_locally_non_stable_gen + base_step;
    }

    bool try_ps(ProblemState &ps, BeamSearchContainer &bsc, UniqueEnsurer &ue, uint32_t max_gen, double spanning_tree_cost_to_beat) {
        apg::pattern out_gliders = origin, stabilised = origin; //to be rewritten
        apg::pattern start = origin + b2p(ps.added_objects);
        ps.stable_generation = get_stable_generation(start, max_gen, out_gliders, stabilised, ue);
        if (ps.stable_generation==0) {//does not stabilise
            //std::cerr << "x";
            return false;
        }
        ps.num_output_gliders = out_gliders.totalPopulation() / 5;
        bool is_solution = (stabilised.totalPopulation()==0);

        ps.locally_stable_generation = get_locally_stable_generation(start, ps.stable_generation, stabilised, ue);

        ps.early_generation = ps.locally_stable_generation > ue.origin_period + 256 ? ps.locally_stable_generation - 256 : ue.origin_period;
        ps.early_hash = start[ps.early_generation].flatlayer(0).hash();
        //todo diameter? is this "lab" dependent? if so, flatlayer should be passed to ue and ue should have its lab to read hashes from
        if (ue.leq_cost_update(ps.early_hash, ps.added_objects_cost)) {
            // was overtaken (possibly by a different thread)
            // at least it saves the spanning_tree_cost calculation when bsc would filter it anyways
            // calculating spanning_tree_cost only when new early_hash appears could be a small improvement (reusing the old result)
            //std::cerr << "l";
            return false;
        }
        apg::pattern dummy1=origin, dummy2=origin, stable_envelope=origin; //to be redefined
        early_and_late(stabilised, 0, dummy1, dummy2, (ue.origin_period>2) ? ue.origin_period : 2, stable_envelope);
        //ps.output_gliders = out_gliders.flatlayer(0);
        ps.spanning_tree_cost = spanning_tree_cost_calc(stable_envelope);

        uint32_t output_glider_cost = ps.num_output_gliders * objectsSource[0].object_cost; //cost to "safely" reduce to the clean state
        int independent_work_estimate = 1200 + ps.total_cost(ue.pessimism) - ps.added_objects_cost; // expected cost for restarted problem
        int solution_clarity = is_solution ? (output_glider_cost) : -(output_glider_cost+independent_work_estimate);

        if (ue.best_cost_update(ps.added_objects_cost, solution_clarity) || is_solution) {
            ue.save_solution(ps, start, solution_clarity);
            if (is_solution && (ps.num_output_gliders==0)) {
                //std::cerr << std::endl << "\033[32;1mSolution " << ps.added_objects_cost << "/" << output_glider_cost << ":" << solution_clarity << "\033[0m";
                // trying to add another not cheaper object cannot improve the solution
                return true;
            }
        }

        if (ps.spanning_tree_cost < spanning_tree_cost_to_beat) {// single cluster remaining is so close to solution to try it, but an incentive to have the envelope closer could be required
            bsc.try_insert(ps, ue.pessimism);
        }
        //std::cerr << "i";
        return false;
    }

    bool try_extend_ps(const ProblemState &old_ps, apg::pattern added_object, uint32_t added_object_cost, BeamSearchContainer &bsc, UniqueEnsurer &ue) {
        ProblemState ps;
        ps.added_objects = (b2p(old_ps.added_objects) + added_object).flatlayer(0);
        ps.added_objects_cost = old_ps.added_objects_cost + added_object_cost;
        return try_ps(ps, bsc, ue, old_ps.stable_generation + ue.max_extra_gens, old_ps.spanning_tree_cost + 30);
    }

    /*  at most max_branching times try adding as cheapest so far not tried (in the routine execution) object to the pattern.
        the object must touch envelope of last 256 generations till reaction stabilisation,
        the does not influence earlier parts of reaction
        the object is restricted to be fully in allowed area
        if the object cost added object_cost_total would exceed current best object cost of a solution ... break
        //... always test if reaction p stabilises inside the allowed region (up to at most given number of issued gliders) test uniqueness (or improved objects cost) if so evaluate stt cost and insert to bsc
        // unique ensurer should be incorporated to bsc ... when the pattern matches with smaller cost, the bigger object cost should be replaced (and stt cost reused)
        */
    void generate_successors(const ProblemState &ps, BeamSearchContainer &bsc, UniqueEnsurer &ue) { // UniqueEnsurer &ue on added objects? The method of adding would hardly reach the same configuration by different paths

        uint32_t best_solution_cost = ue.get_best_solution_cost();
        if (best_solution_cost <= ps.added_objects_cost) {
            return;
        }

        apg::pattern early_pattern = origin, early_env = origin, late_env = origin; // to be redefined
        uint32_t late_generation = ps.early_generation + 512 < ps.stable_generation ? ps.early_generation + 512 : ps.stable_generation;

        early_and_late(origin+b2p(ps.added_objects), ps.early_generation, early_pattern, early_env, late_generation, late_env);

        uint32_t maximal_allowed_object_cost = best_solution_cost - ps.added_objects_cost;
        uint32_t last_cost = 0;
        apg::pattern current_objects_forbidden = objects_forbidden + early_env.convolve(influence);
        apg::pattern objects_allowed = reaction_allowed - current_objects_forbidden;
        apg::pattern to_touch = late_env.convolve(influence) & objects_allowed;
        /*  probably not worth ... the cycle would be fast anyways
            if (to_touch.empty()) {
                return;
            }
        */

        uint32_t branching_remaining = ue.max_branching;
        for(ObjectWithCost &owc : objectsSource) {
            last_cost = owc.object_cost;
            if ((branching_remaining == 0) || (owc.object_cost > maximal_allowed_object_cost)) {
                    break;
            }
            // find bitmap of lucorners where the object could be inserted
            apg::pattern lucorners_to_touch = to_touch.convolve(owc.object_env_flipped);
            apg::pattern collisions = lucorners_to_touch.convolve(owc.object_env) - objects_allowed;
            apg::pattern lucorners = lucorners_to_touch - collisions.convolve(owc.object_env_flipped);

            apg::bitworld lu_bw = lucorners.flatlayer(0); int64_t bbox[4] = {0};
            while (lu_bw.population()) {
                apg::bitworld lu_cell = lu_bw.get1cell();lu_bw -= lu_cell;lu_cell.getbbox(bbox);
                // if a search edge leads to no glider solution set branching_remaining to 1 to finish the loop early (there will not be cheaper no glider solution and cheap glider solutions are search sideeffects, not the goal)
                if (try_extend_ps(ps, owc.object.shift(bbox[0], bbox[1]), owc.object_cost, bsc, ue)) {
                    branching_remaining=1;
                }
                branching_remaining--;
                if (branching_remaining==0) {
                    break;
                }
            }
        }
        std::cerr << "(" << last_cost << ":" << branching_remaining << ")";
    }
};

void run_worker(std::atomic<uint32_t> *ctr, uint32_t n_tasks, SODThread *thread, BeamSearchContainer *bsc, UniqueEnsurer *ue, const ProblemState *problems) {

    for (;;) {
        // dequeue subtask:
        uint32_t idx = (uint32_t) ((*ctr)++);
        if (idx >= n_tasks) { break; }
        std::cerr << "[" << idx << "]";
        thread->generate_successors(problems[idx], *bsc, *ue);
    }
}

void run1iter(std::vector<SODThread> &gsodv, uint32_t beamwidth, UniqueEnsurer &ue, std::vector<ProblemState> &problems, uint32_t depth) {

    uint32_t n_tasks = problems.size();
    std::atomic<uint32_t> ctr{0ull};
    int parallelism = gsodv.size();

    std::cerr << "Processing " << n_tasks << " tasks on " << parallelism << " threads..." << std::endl;

    std::vector<BeamSearchContainer> bscv(parallelism);
    std::vector<std::thread> workers;

    for (int i = 0; i < parallelism; i++) {
        bscv[i].maxsize = beamwidth;
        workers.emplace_back(run_worker, &ctr, n_tasks, &(gsodv[i]), &(bscv[i]), &ue, &(problems[0]));
    }

    for (auto&& w : workers) { w.join(); }

    BeamSearchContainer bsc; bsc.maxsize = beamwidth;
    for (int i = 0; i < parallelism; i++) { bsc.add(bscv[i],ue.pessimism); }

    problems.clear();
    int beamindex=0; bool first=true, saveit=(depth % BIG_SAVE_FREQ)==0;
    for (auto it = bsc.pmpq.begin(); it != bsc.pmpq.end(); ++it) {
        ProblemState ps = bsc.contents[it->second];
        if (saveit || first) {
            first = false;
            ue.save_progress(gsodv[0].origin+gsodv[0].b2p(ps.added_objects), depth, ++beamindex, ps.spanning_tree_cost, ps.added_objects_cost);
        }
        problems.push_back(ps);
    }

    if (bsc.pmpq.empty()) {
        std::cerr << "search complete." << std::endl;
    }

    if (ue.best_solution_cost==999999999) {//solution not found yet
        if (saveit) {
            beamindex = beamwidth; // current beam could be shorter
            for(uint32_t d=2; d <= depth; d++) {
                if (((depth+1-d) % BIG_SAVE_FREQ)==0) {
                    int treshold=1+(beamwidth/d);
                    //std::cerr << "removal beamindex " << beamindex << " treshold " << treshold << std::endl;
                    while (beamindex>treshold) {
                        std::string delfilename = "SoD_Progress_"+std::to_string(depth+1-d)+"_"+std::to_string(beamindex--)+".mc";
                        std::remove(delfilename.c_str());
                    };
                    if (beamindex==1) {
                        break;
                    }
                }
            }
        }
    }
}

void glider_worker(std::atomic<uint32_t> *ctr, SODThread *thread, BeamSearchContainer *bsc, UniqueEnsurer *ue, double start_spanning_tree_cost_plus) {

    // starting glider cost 1200 corresponds to expected cost for creating/redirecting the input glider

    int64_t bbox[4] = {0, 0, 0, 0}; thread->origin.getrect(bbox);
    uint32_t n_lanes_onedir = (bbox[2]+bbox[3]+8);
    uint32_t n_tasks = 2 * n_lanes_onedir;
    std::vector<bool> n_allowed(n_tasks);std::vector<bool> s_allowed(n_tasks);
    for(uint32_t i=0;i<n_tasks;i++) {
        n_allowed[i]=true;s_allowed[i]=true;
    }
    apg::bitworld forbidden_glider_lanes = ue->input_glider_lines_forbidden; //x,y projections may not intersect bbox projections (otherwise ignored)
    int64_t cell_bbox[4] = {0, 0, 0, 0};
    while (forbidden_glider_lanes.population()) {
        apg::bitworld cell = forbidden_glider_lanes.get1cell();forbidden_glider_lanes -= cell;
        cell.getbbox(cell_bbox);
        if (cell_bbox[1]<bbox[1]) {
            if (cell_bbox[0]<bbox[0]) {
                int lane = cell_bbox[1]-cell_bbox[0];
                if (lane >= (bbox[1]-bbox[0]-bbox[2]-4)) {
                    uint32_t l_idx = lane - (bbox[1]-bbox[0]-bbox[2]-4);
                    if (l_idx < n_lanes_onedir) {
                        s_allowed[n_lanes_onedir+l_idx] = false;
                    }
                }
            }
            else if ((cell_bbox[0] > bbox[0]+bbox[2])) {
                int lane = cell_bbox[1]+cell_bbox[0];
                if (lane >= (bbox[1]+bbox[0]-7)) {
                    uint32_t l_idx = lane - (bbox[1]+bbox[0]-7);
                    if (l_idx < n_lanes_onedir) {
                        s_allowed[l_idx] = false;
                    }
                }
            }
        }
        else if ((cell_bbox[1] > bbox[1]+bbox[3])) {
            if (cell_bbox[0]<bbox[0]) {
                int lane = cell_bbox[1]+cell_bbox[0];
                if (lane >= (bbox[1]+bbox[0]-7)) {
                    uint32_t l_idx = lane - (bbox[1]+bbox[0]-7);
                    if (l_idx < n_lanes_onedir) {
                        n_allowed[l_idx] = false;
                    }
                }
            }
            else if ((cell_bbox[0] > bbox[0]+bbox[2])) {
                int lane = cell_bbox[1]-cell_bbox[0];
                if (lane >= (bbox[1]-bbox[0]-bbox[2]-4)) {
                    uint32_t l_idx = lane - (bbox[1]-bbox[0]-bbox[2]-4);
                    if (l_idx < n_lanes_onedir) {
                        n_allowed[n_lanes_onedir+l_idx] = false;
                    }
                }
            }
        }
    }
    apg::pattern glidernw(&(thread->lt), "2o$obo$o!", "b3s23");
    apg::pattern gliderse(&(thread->lt), "bo$2bo$3o!", "b3s23");
    apg::pattern gliderne(&(thread->lt), "b2o$obo$2bo!", "b3s23");
    apg::pattern glidersw(&(thread->lt), "bo$o$3o!", "b3s23");
    for (;;) {
        // dequeue subtask:
        uint32_t idx = (uint32_t) ((*ctr)++);
        int lane, x, y;
        if (idx >= n_tasks) { break; } //xpy lines and ymx lines ... for each of them gliders in both directions with origin period different stat phases should be tried ... each "unique stable pattern is added to bsc" with a glider in ps
        if (idx >= n_lanes_onedir) {//y-x lanes
            lane = idx - n_lanes_onedir + bbox[1]-bbox[0]-bbox[2]-4;
            if (n_allowed[idx]) {
                for(uint32_t ph=0;ph<ue->origin_period;ph++) {
                    ProblemState psnw;
                    psnw.added_objects_cost=1200;
                    if (lane > bbox[1]+bbox[3]-bbox[0]-bbox[2]) {
                        y = bbox[1] + bbox[3] + ue->origin_period/4 + 23;
                        x = y - lane;
                    } else {
                        x = bbox[0] + bbox[2] + ue->origin_period/4 + 23;
                        y = x + lane;
                    }
                    psnw.added_objects = (glidernw[ph].shift(x,y)).flatlayer(0);
                    thread->try_ps(psnw, *bsc, *ue, ue->max_extra_gens, start_spanning_tree_cost_plus);
                }
            }
            if (s_allowed[idx]) {
                for(uint32_t ph=0;ph<ue->origin_period;ph++) {
                    ProblemState psse;
                    psse.added_objects_cost=1200;
                    if (lane > bbox[1]-bbox[0]) {
                        x = bbox[0] - ue->origin_period/4 - 26;
                        y = x + lane;
                    } else {
                        y = bbox[1] - ue->origin_period/4 - 26;
                        x = y - lane;
                    }
                    psse.added_objects = (gliderse[ph].shift(x,y)).flatlayer(0);
                    thread->try_ps(psse, *bsc, *ue, ue->max_extra_gens, start_spanning_tree_cost_plus);
                }
            }
        }
        else {//x+y lanes
            lane = idx + bbox[1]+bbox[0]-7;
            if (n_allowed[idx]) {
                for(uint32_t ph=0;ph<ue->origin_period;ph++) {
                    ProblemState psne;
                    psne.added_objects_cost=1200;
                    if (lane > bbox[1]+bbox[3]+bbox[0]) {
                        y = bbox[1] + bbox[3] + ue->origin_period/4 + 23;
                        x = lane - y;
                    } else {
                        x = bbox[0] - ue->origin_period/4 - 26;
                        y = lane - x;
                    }
                    psne.added_objects = (gliderne[ph].shift(x,y)).flatlayer(0);
                    thread->try_ps(psne, *bsc, *ue, ue->max_extra_gens, start_spanning_tree_cost_plus);
                }
            }
            if (s_allowed[idx]) {
                for(uint32_t ph=0;ph<ue->origin_period;ph++) {
                    ProblemState pssw;
                    pssw.added_objects_cost=1200;
                    if (lane > bbox[0]+bbox[2]+bbox[1]) {
                        x = bbox[0]+bbox[2] + ue->origin_period/4 + 23;
                        y = lane - x;
                    } else {
                        y = bbox[1] - ue->origin_period/4 - 26;
                        x = lane - y;
                    }
                    pssw.added_objects = (glidersw[ph].shift(x,y)).flatlayer(0);
                    thread->try_ps(pssw, *bsc, *ue, ue->max_extra_gens, start_spanning_tree_cost_plus);
                }
            }
        }
    }
}

void start_gsod(std::vector<SODThread> &gsodv, uint32_t beamwidth, UniqueEnsurer &ue, std::vector<ProblemState> &problems) {
    std::cerr << "Starting gliders search..." << std::flush;

    int64_t bbox[4] = {0, 0, 0, 0}; gsodv[0].origin.getrect(bbox);
    uint32_t n_tasks = 2 * (bbox[2]+bbox[3]+8);

    std::atomic<uint32_t> ctr{0ull};
    int parallelism = gsodv.size();

    apg::pattern dummy1=gsodv[0].origin, dummy2=gsodv[0].origin, stable_envelope=gsodv[0].origin; //to be redefined
    gsodv[0].early_and_late(gsodv[0].origin, 0, dummy1, dummy2, (ue.origin_period>2) ? ue.origin_period : 2, stable_envelope);
    double start_spanning_tree_cost_to_beat = 0.1 + gsodv[0].spanning_tree_cost_calc(stable_envelope);

    std::cerr << "Running " << n_tasks << " lane tasks on " << parallelism << " threads... spanning tree cost to beat " <<  start_spanning_tree_cost_to_beat << std::endl;

    std::vector<BeamSearchContainer> bscv(parallelism);
    std::vector<std::thread> workers;


    for (int i = 0; i < parallelism; i++) {
        bscv[i].maxsize = beamwidth;
        workers.emplace_back(glider_worker, &ctr, &(gsodv[i]), &(bscv[i]), &ue, start_spanning_tree_cost_to_beat);
    }

    std::cerr << "Waiting for workers" << std::endl;

    BeamSearchContainer bsc; bsc.maxsize = beamwidth;
    for (auto&& w : workers) { w.join(); }

    std::cerr << "Workers finished" << std::endl;

    for (int i = 0; i < parallelism; i++) {
        bsc.add(bscv[i], ue.pessimism);
    }
    problems.clear();
    int beamindex=0;
    for (auto it = bsc.pmpq.begin(); it != bsc.pmpq.end(); ++it) {
        ProblemState ps = bsc.contents[it->second];
        ue.save_progress(gsodv[0].origin+gsodv[0].b2p(ps.added_objects), 0, ++beamindex, ps.spanning_tree_cost, 0);
        problems.push_back(ps);
    }
}

void run_main(  std::string infile,
                uint32_t beamwidth,
                uint32_t parallelism,
                uint32_t max_output_gliders_allowed,
                uint32_t max_extra_gens,
                uint32_t max_branching,
                uint32_t max_object_types_allowed,
                uint32_t pessimism) {

    std::vector<SODThread> gsodv(parallelism);
    std::vector<ProblemState> problems(1);

    uint32_t depth=0;
    std::cerr << "Loading problem..." << std::flush;
    std::vector<apg::bitworld> vbw = load_file<4>(infile);
    std::cerr << "Parsing problem..." << std::flush;
    UniqueEnsurer ue=gsodv[0].parse_pattern(vbw);
    ue.max_output_gliders_allowed = max_output_gliders_allowed;
    ue.max_extra_gens = max_extra_gens;
    ue.max_branching = max_branching;
    ue.max_object_types_allowed = max_object_types_allowed;
    ue.pessimism = pessimism;

    std::cerr << "Starting treads..." << std::flush;
    for (uint32_t i = 0; i < parallelism; i++) {
        gsodv[i].init_patterns(ue);
        //gsodv[i].power_thousands = 500; //(i==0) ? -1 : 700;
    }

    if (ue.ori_added_objects.population()) {
        ProblemState ps;
        ps.added_objects_cost=0;//
        ps.added_objects = ue.ori_added_objects;
        BeamSearchContainer bsc; bsc.maxsize = beamwidth;
        gsodv[0].try_ps(ps, bsc, ue, ue.max_extra_gens, 999999999);
        problems.clear();
        for (auto it = bsc.pmpq.begin(); it != bsc.pmpq.end(); ++it) {
            ProblemState ps = bsc.contents[it->second];
            problems.push_back(ps);
        }
    } else {
        std::cerr << "Looking for starting gliders ... " << std::flush;
        start_gsod(gsodv, beamwidth, ue, problems);
    }
    std::cerr << "Starting adding objects ... " << std::flush;
    while (!(problems.empty())) {
        std::cerr << " (" << depth << ")" ;
        run1iter(gsodv, beamwidth, ue, problems, ++depth);
    }
}

} // namespace gsod

int main(int argc, char* argv[]) {

    bool incorrectly_called=false;
    std::ostringstream err;
    if ((argc<4)||(argc>9)) {incorrectly_called = true; err << "Number of arguments = " << argc << "!";}
    for (int i = 2; i < argc; i++) {
        if (argv[i][0] == '-') { incorrectly_called = true; err << "Negative argument " << i << "!"; break;}
    }
    if (!incorrectly_called) {
        std::string filename = argv[1];
        std::ifstream f(filename.c_str());
        if (!f.good()) {
            incorrectly_called = true;
            err << "File "<< filename << " does no exist!";
        }
    }
    if (incorrectly_called) {
        std::cerr << err.str() << std::endl << "Correct call: gSoD <pattern file> <beamwidth> <parallelism> [<max_output_gliders_allowed>(200)] [<max_extra_gens>(1024)] [<max_branching>(2048)] [<max_object_types_allowed>10)] [<pessimism>70]" << std::endl;
        return 1;
    }

    uint32_t beamwidth = std::stoi(argv[2]);
    uint32_t parallelism = std::stoi(argv[3]);
    uint32_t max_output_gliders_allowed = 200;
    if (argc > 4) {
        max_output_gliders_allowed = std::stoi(argv[4]);
    }
    uint32_t max_extra_gens = 1024;
    if (argc > 5) {
        max_extra_gens = std::stoi(argv[5]);
    }
    uint32_t max_branching = 2048;
    if (argc > 6) {
        max_branching = std::stoi(argv[6]);
    }
    uint32_t max_object_types_allowed = 10;
    if (argc > 7) {
        max_object_types_allowed = std::stoi(argv[7]);
    }
    uint32_t pessimism = 70;
    if (argc > 8) {
        pessimism = std::stoi(argv[8]);
    }

    gsod::run_main(argv[1], beamwidth, parallelism, max_output_gliders_allowed, max_extra_gens, max_branching, max_object_types_allowed, pessimism);

    return 0;
}
