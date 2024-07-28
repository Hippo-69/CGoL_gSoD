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

#define BIG_SAVE_FREQ 10

namespace gsod {

template<int N>
std::vector<apg::bitworld> load_file(std::string filename) {
    std::vector<apg::bitworld> vbw;
    apg::lifetree<uint32_t, N> lt(100);
    apg::pattern pat = apg::pattern(&lt, filename);
    for (int i = 0; i < N; i++) { vbw.push_back(pat.flatlayer(i)); }
    return vbw;
}

struct SODThread {

    apg::lifetree<uint32_t, 1> lt;

    apg::pattern origin; // life pattern to destroy
    apg::pattern reaction_allowed; // includes both reaction and possible objects placement
    apg::pattern objects_forbidden; // objects cannot be placed here as it would react with pattern history
    //I should allow to restrict some input glider lanes+directions.
    apg::pattern influence; //handy constant
    std::vector<ObjectWithCost> objectsSource; //all p at most 2

    SODThread() : lt(1000), origin(&lt, "", "b3s23"), reaction_allowed(&lt, "", "b3s23"), objects_forbidden(&lt, "", "b3s23"), influence(&lt, "b3o$5o$5o$5o$b3o!", "b3s23") {
        influence = influence.shift(-2, -2);
    }

    apg::pattern b2p(const apg::bitworld &bw) {
        std::vector<apg::bitworld> vbw(1, bw);
        apg::pattern pat(&lt, vbw, "b3s23");
        return pat;
    }

    void init_patterns(apg::bitworld _origin, apg::bitworld _reaction_allowed, apg::bitworld _objects_forbidden) {//each thread have different lt, so must have its separate pattern "constants"
        origin = b2p(_origin); reaction_allowed = b2p(_reaction_allowed); objects_forbidden = b2p(_objects_forbidden);

        apg::pattern blinker(&lt, "3o!", "b3s23");
            objectsSource.push_back({blinker, 467});
            objectsSource.push_back({blinker[1], 467});
        apg::pattern block(&lt, "2o$2o!", "b3s23");
            objectsSource.push_back({block, 493});
        apg::pattern beehive(&lt, "b2o$o2bo$b2o!", "b3s23");
            objectsSource.push_back({beehive, 615});
            objectsSource.push_back({beehive.transform("rcw", 0, 0), 615});
        apg::pattern pond(&lt, "b2o$o2bo$o2bo$b2o!", "b3s23");
            objectsSource.push_back({pond, 685});
        apg::pattern loaf(&lt, "2bo$bobo$o2bo$b2o!", "b3s23");
            objectsSource.push_back({loaf, 714});
            objectsSource.push_back({loaf.transform("rcw", 0, 0), 714});
            objectsSource.push_back({loaf.transform("flip", 0, 0), 714});
            objectsSource.push_back({loaf.transform("rccw", 0, 0), 714});
        apg::pattern tub(&lt, "bo$obo$bo!", "b3s23");
            objectsSource.push_back({tub, 726});
        apg::pattern boat(&lt, "2o$obo$bo!", "b3s23");
            objectsSource.push_back({boat, 743});
            objectsSource.push_back({boat.transform("rcw", 0, 0), 743});
            objectsSource.push_back({boat.transform("flip", 0, 0), 743});
            objectsSource.push_back({boat.transform("rccw", 0, 0), 743});
        apg::pattern ship(&lt, "2o$obo$b2o!", "b3s23");
            objectsSource.push_back({ship, 767});
            objectsSource.push_back({ship.transform("rcw", 0, 0), 767});
        apg::pattern barge(&lt, "bo$obo$bobo$2bo!", "b3s23");
            objectsSource.push_back({barge, 816});
            objectsSource.push_back({barge.transform("rcw", 0, 0), 816});
        apg::pattern longboat(&lt, "2o$obo$bobo$2bo!", "b3s23");
            objectsSource.push_back({longboat, 868});
            objectsSource.push_back({longboat.transform("rcw", 0, 0), 868});
            objectsSource.push_back({longboat.transform("flip", 0, 0), 868});
            objectsSource.push_back({longboat.transform("rccw", 0, 0), 868});
        apg::pattern mango(&lt, "b2o$o2bo$bo2bo$2b2o!", "b3s23");
            objectsSource.push_back({mango, 1010});
            objectsSource.push_back({mango.transform("rcw", 0, 0), 1010});
        apg::pattern shiptie(&lt, "2o$obo$b2o$3b2o$3bobo$4b2o!", "b3s23");
            objectsSource.push_back({shiptie, 1050});
            objectsSource.push_back({shiptie.transform("rcw", 0, 0), 1050});
        apg::pattern eater1(&lt, "2o$obo$2bo$2b2o!", "b3s23");
            objectsSource.push_back({eater1, 1120});
            objectsSource.push_back({eater1.transform("rcw", 0, 0), 1120});
            objectsSource.push_back({eater1.transform("flip", 0, 0), 1120});
            objectsSource.push_back({eater1.transform("rccw", 0, 0), 1120});
            objectsSource.push_back({eater1.transform("swap_xy", 0, 0), 1120});
            objectsSource.push_back({eater1.transform("flip_x", 0, 0), 1120});
            objectsSource.push_back({eater1.transform("flip_y", 0, 0), 1120});
            objectsSource.push_back({eater1.transform("swap_xy_flip", 0, 0), 1120});
        /*
            apg::pattern aircraftcarrier(&lt, "2o$o$2bo$b2o!", "b3s23");
                objectsSource.push_back({aircraftcarrier, 1370});
                objectsSource.push_back({aircraftcarrier.transform("rcw", 0, 0), 1370});
                objectsSource.push_back({aircraftcarrier.transform("flip_x", 0, 0), 1370});
                objectsSource.push_back({aircraftcarrier.transform("swap_xy", 0, 0), 1370});
            apg::pattern elevener(&lt, "2o$obo$2bo$2b3o$5bo$4b2o!", "b3s23");
                objectsSource.push_back({elevener, 1430});
                objectsSource.push_back({elevener.transform("rcw", 0, 0), 1430});
                objectsSource.push_back({elevener.transform("flip", 0, 0), 1430});
                objectsSource.push_back({elevener.transform("rccw", 0, 0), 1430});
            apg::pattern intergralsign(&lt, "2o$obo$2bo$2bobo$3b2o!", "b3s23");
                objectsSource.push_back({intergralsign, 1450});
                objectsSource.push_back({intergralsign.transform("rcw", 0, 0), 1450});
                objectsSource.push_back({intergralsign.transform("flip_x", 0, 0), 1450});
                objectsSource.push_back({intergralsign.transform("swap_xy", 0, 0), 1450});
            objectsSource.push_back(apg::pattern shillelagh(&lt, "2o$obo$2bo$bo$b2o!", "b3s23");
                objectsSource.push_back({{shillelagh, 1470})
                objectsSource.push_back({shillelagh.transform("rcw", 0, 0), 1470});
                objectsSource.push_back({shillelagh.transform("flip", 0, 0), 1470});
                objectsSource.push_back({shillelagh.transform("rccw", 0, 0), 1470});
                objectsSource.push_back({shillelagh.transform("swap_xy", 0, 0), 1470});
                objectsSource.push_back({shillelagh.transform("flip_x", 0, 0), 1470});
                objectsSource.push_back({shillelagh.transform("flip_y", 0, 0), 1470});
                objectsSource.push_back({shillelagh.transform("swap_xy_flip", 0, 0), 1470});
            apg::pattern transblockonlonghook(&lt, "2ob2o$2obo$3bo$3bobo$4b2o!", "b3s23");
                objectsSource.push_back({transblockonlonghook, 1480});
                objectsSource.push_back({transblockonlonghook.transform("rcw", 0, 0), 1480});
                objectsSource.push_back({transblockonlonghook.transform("flip", 0, 0), 1480});
                objectsSource.push_back({transblockonlonghook.transform("rccw", 0, 0), 1480});
                objectsSource.push_back({transblockonlonghook.transform("swap_xy", 0, 0), 1480});
                objectsSource.push_back({transblockonlonghook.transform("flip_x", 0, 0), 1480});
                objectsSource.push_back({transblockonlonghook.transform("flip_y", 0, 0), 1480});
                objectsSource.push_back({transblockonlonghook.transform("swap_xy_flip", 0, 0), 1480});
        */
    }

    void early_and_late(apg::pattern x, uint64_t early_gens, apg::pattern &early_pattern, apg::pattern &early_env, uint64_t late_gens, apg::pattern &late_env) {

        apg::pattern y = x;
        apg::pattern z = x;

        for (uint64_t i = 0; i < late_gens; i++) {
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

    // parsing lifehistory pattern with history of the circuit. It can contain state 4 extension marking reaction_allowed
    // otherwise a default neighbourhood of the pattern is chosen (not implemented yet)
    UniqueEnsurer parse_pattern(const std::vector<apg::bitworld> &vbw) {
        for(uint32_t l=0;l<vbw.size();l++) {// debug
            std::ofstream out("Input_plane_"+std::to_string(l)+".mc");
            b2p(vbw[l]).write_macrocell(out);
            out.close();
        }
        apg::pattern origin_patt = b2p(vbw[0]);
        apg::pattern reaction_allowed_patt = b2p(vbw[0])+b2p(vbw[1])+b2p(vbw[2]);//to be debugged
        if (b2p(vbw[1]).empty()) {
            //todo reaction_allowed_patt = approx_convex_hull(reaction_allowed_patt);
            apg::pattern halo2 = influence;
            halo2 = halo2.convolve(halo2).convolve(halo2); // radius 6
            halo2 = halo2.convolve(halo2); // radius 12
            reaction_allowed_patt = reaction_allowed_patt.convolve(halo2);
        }
        apg::pattern objects_forbidden_patt = (b2p(vbw[0])+b2p(vbw[2])).convolve(influence);// to be debugged

        UniqueEnsurer ue(vbw[0],reaction_allowed_patt.flatlayer(0),objects_forbidden_patt.flatlayer(0));
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
        double stlength = apg::spanning_graph_to_tree(ccentres, 0, &std::sqrt, all_edges); // actually constant 2.0 scale does not matter
        // if one wants to follow simeks program std::sqrt should be replaced by its 5/4 power (or better calculate the tree and then recompute the cost using this evaluation).
        // I would first experiment with the simple euklid distance mst.
        return stlength;
    }

    uint64_t get_stable_generation(apg::pattern starting_pattern, uint64_t max_gen, apg::pattern &out_gliders, apg::pattern &stabilised, UniqueEnsurer &ue) {
        uint64_t base_step = (ue.origin_period>8) ? ue.origin_period : 8;
        apg::pattern curr = starting_pattern, curr_base = curr & reaction_allowed;
        for(uint64_t gen=0;gen<max_gen;gen+=base_step) {
            apg::pattern next = curr[base_step], next_base = next & reaction_allowed;
            if (next_base == curr_base) {
                stabilised = curr_base;
                //if (!((next & curr) - next_base).empty()) {return 0;} ... it could miss a case with gliders intersecting, but it would avoid stable trash of correct size
                return gen;
            }
            out_gliders = (next - next_base).convolve(influence) & next; // part of the gliders could still be in allowed
            uint64_t gliders_pop = out_gliders.totalPopulation();
            if (((gliders_pop % 5)!=0) || (gliders_pop / 5 > ue.max_gliders_allowed)) {
                return 0; // forbidden
                // too many gliders or gliders not formed cleanly before leaving allowed
                // there could still be unprobable case of objects of propper size created out of allowed with propper sizes in each checkpoint on the way to fake extra solutions
                // I would risk it (could be made safer later)
            }
            curr = next; curr_base = next_base;
        }
        return 0;
    }

    bool try_ps(ProblemState &ps, BeamSearchContainer &bsc, UniqueEnsurer &ue, uint64_t max_gen) {
        apg::pattern out_gliders = origin, stabilised = origin; //to be rewritten
        apg::pattern start = origin + b2p(ps.added_objects);
        ps.stable_generation = get_stable_generation(start, max_gen, out_gliders, stabilised, ue);
        if (ps.stable_generation==0) {//does not stabilise
            return false;
        }
        ps.num_output_gliders = out_gliders.totalPopulation() / 5;
        bool is_solution = (stabilised.totalPopulation()==0);
        ps.early_generation = ps.stable_generation > ue.origin_period + 256 ? ps.stable_generation - 256 : ue.origin_period;
        ps.early_hash = (origin+b2p(ps.added_objects))[ps.early_generation].totalHash(1000);//todo diameter? is this "lab" dependent? if so, flatlayer should be passed to ue and ue should have its lab to read hashes from
        if (ue.leq_cost_update(ps.early_hash, ps.added_objects_cost, is_solution && (ps.num_output_gliders==0))) {
            // was overtaken (possibly by a different thread)
            // at least it saves the spanning_tree_cost calculation when bsc would filter it anyways
            // calculating spanning_tree_cost only when new early_hash appears could be a small improvement (reusing the old result)
            return false;
        }
        if (is_solution) {
            ue.save_solution(ps, start);
            if (ps.num_output_gliders==0) {
                return true;
            }
        }
        apg::pattern dummy1=origin, dummy2=origin, stable_envelope=origin; //to be redefined
        early_and_late(stabilised, 0, dummy1, dummy2, (ue.origin_period>2) ? ue.origin_period : 2, stable_envelope);
        //ps.output_gliders = out_gliders.flatlayer(0);
        ps.spanning_tree_cost = spanning_tree_cost_calc(stable_envelope);
        bsc.try_insert(ps);
        return false;
    }

    bool try_extend_ps(const ProblemState &old_ps, apg::pattern added_object, uint32_t added_object_cost, BeamSearchContainer &bsc, UniqueEnsurer &ue) {
        ProblemState ps;
        ps.added_objects = (b2p(old_ps.added_objects) + added_object).flatlayer(0);
        ps.added_objects_cost = old_ps.added_objects_cost + added_object_cost;
        return try_ps(ps, bsc, ue, old_ps.stable_generation + ue.max_extra_gens);
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

        early_and_late(origin+b2p(ps.added_objects), ps.early_generation, early_pattern, early_env, ps.stable_generation, late_env);

        uint32_t maximal_allowed_object_cost = best_solution_cost - ps.added_objects_cost;

        apg::pattern current_objects_forbidden = objects_forbidden + early_env.convolve(influence);
        apg::pattern objects_allowed = reaction_allowed - current_objects_forbidden;
        apg::pattern to_touch = late_env.convolve(influence) & objects_allowed;

        uint64_t branching_remaining = ue.max_branching;
        for(ObjectWithCost &owc : objectsSource) {
            if ((branching_remaining == 0) || (owc.object_cost > maximal_allowed_object_cost)) {
                    break;
            }
            // find bitmap of lucorners where the object could be inserted
            apg::pattern object_to_add_env = owc.object + owc.object[1], object_to_add_flip = owc.object.transform("flip", 0, 0);

            apg::pattern lucorners_to_touch = to_touch.convolve(object_to_add_flip);
            apg::pattern collisions = lucorners_to_touch.convolve(object_to_add_env) - objects_allowed;
            apg::pattern lucorners = lucorners_to_touch - collisions.convolve(object_to_add_flip);

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
    }
};

void run_worker(std::atomic<uint64_t> *ctr, uint64_t n_tasks, SODThread *thread, BeamSearchContainer *bsc, UniqueEnsurer *ue, const ProblemState *problems) {

    for (;;) {
        // dequeue subtask:
        uint64_t idx = (uint64_t) ((*ctr)++);
        if (idx >= n_tasks) { break; }
        thread->generate_successors(problems[idx], *bsc, *ue);
    }
}

void run1iter(std::vector<SODThread> &gsodv, uint32_t beamwidth, UniqueEnsurer &ue, std::vector<ProblemState> &problems, uint32_t depth) {

    uint64_t n_tasks = problems.size();
    std::atomic<uint64_t> ctr{0ull};
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
    for (int i = 0; i < parallelism; i++) { bsc += bscv[i]; }

    problems.clear();
    int beamindex=0; bool first=true, saveit=(depth % BIG_SAVE_FREQ)==0;
    for (auto it = bsc.pmpq.begin(); it != bsc.pmpq.end(); ++it) {
        ProblemState ps = bsc.contents[it->second];
        if (saveit || first) {
            first = false;
            ue.save_progress(gsodv[0].origin+gsodv[0].b2p(ps.added_objects), depth, ++beamindex);
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

void glider_worker(std::atomic<uint64_t> *ctr, SODThread *thread, BeamSearchContainer *bsc, UniqueEnsurer *ue) {

    int64_t bbox[4] = {0, 0, 0, 0}; thread->origin.getrect(bbox);
    uint64_t n_lanes_onedir = (bbox[2]+bbox[3]+8);
    uint64_t n_tasks = 2 * n_lanes_onedir;
    for (;;) {
        // dequeue subtask:
        uint64_t idx = (uint64_t) ((*ctr)++);
        int lane, x, y;
        if (idx >= n_tasks) { break; } //xpy lines and ymx lines ... for each of them gliders in both directions with origin period different stat phases should be tried ... each "unique stable pattern is added to bsc" with a glider in ps
        if (idx >= n_lanes_onedir) {//y-x lanes
            lane = idx - n_lanes_onedir + bbox[1]-bbox[0]-bbox[2]-4;
            apg::pattern glidernw(&(thread->lt), "2o$obo$o!", "b3s23");
            apg::pattern gliderse(&(thread->lt), "bo$2bo$3o!", "b3s23");
            for(uint64_t ph=0;ph<ue->origin_period;ph++) {
                ProblemState psnw;
                psnw.added_objects_cost=0;
                if (lane > bbox[1]+bbox[3]-bbox[0]-bbox[2]) {
                    y = bbox[1] + bbox[3] + ue->origin_period/4 + 3;
                    x = y - lane;
                } else {
                    x = bbox[0] + bbox[2] + ue->origin_period/4 + 3;
                    y = x + lane;
                }
                psnw.added_objects = (glidernw[ph].shift(x,y)).flatlayer(0);
                thread->try_ps(psnw, *bsc, *ue, ue->max_extra_gens);
                ProblemState psse;
                psse.added_objects_cost=0;
                if (lane > bbox[1]-bbox[0]) {
                    x = bbox[0] - ue->origin_period/4 - 6;
                    y = x + lane;
                } else {
                    y = bbox[1] - ue->origin_period/4 - 6;
                    x = y - lane;
                }
                psse.added_objects = (gliderse[ph].shift(x,y)).flatlayer(0);
                thread->try_ps(psse, *bsc, *ue, ue->max_extra_gens);
            }
        }
        else {//x+y lanes
            lane = idx + bbox[1]+bbox[0]-7;
            apg::pattern gliderne(&(thread->lt), "b2o$obo$2bo!", "b3s23");
            apg::pattern glidersw(&(thread->lt), "bo$o$3o!", "b3s23");
            for(uint64_t ph=0;ph<ue->origin_period;ph++) {
                ProblemState psne;
                psne.added_objects_cost=0;
                if (lane > bbox[1]+bbox[3]+bbox[0]) {
                    y = bbox[1] + bbox[3] + ue->origin_period/4 + 3;
                    x = lane - y;
                } else {
                    x = bbox[0] - ue->origin_period/4 - 6;
                    y = lane - x;
                }
                psne.added_objects = (gliderne[ph].shift(x,y)).flatlayer(0);
                thread->try_ps(psne, *bsc, *ue, ue->max_extra_gens);
                ProblemState pssw;
                pssw.added_objects_cost=0;
                if (lane > bbox[0]+bbox[2]+bbox[1]) {
                    x = bbox[0]+bbox[2] + ue->origin_period/4 + 3;
                    y = lane - x;
                } else {
                    y = bbox[1] - ue->origin_period/4 - 6;
                    x = lane - y;
                }
                pssw.added_objects = (glidersw[ph].shift(x,y)).flatlayer(0);
                thread->try_ps(pssw, *bsc, *ue, ue->max_extra_gens);
            }
        }
    }
}

void start_gsod(std::vector<SODThread> &gsodv, uint32_t beamwidth, UniqueEnsurer &ue, std::vector<ProblemState> &problems) {
    std::cerr << "Starting gliders search..." << std::flush;

    int64_t bbox[4] = {0, 0, 0, 0}; gsodv[0].origin.getrect(bbox);
    uint64_t n_tasks = 2 * (bbox[2]+bbox[3]+8);

    std::atomic<uint64_t> ctr{0ull};
    int parallelism = gsodv.size();

    std::cerr << "Running " << n_tasks << " lane tasks on " << parallelism << " threads..." << std::endl;

    std::vector<BeamSearchContainer> bscv(parallelism);
    std::vector<std::thread> workers;

    for (int i = 0; i < parallelism; i++) {
        bscv[i].maxsize = beamwidth;
        workers.emplace_back(glider_worker, &ctr, &(gsodv[i]), &(bscv[i]), &ue);
    }

    std::cerr << "Waiting for workers" << std::endl;

    BeamSearchContainer bsc; bsc.maxsize = beamwidth;
    for (auto&& w : workers) { w.join(); }

    std::cerr << "Workers finished" << std::endl;

    for (int i = 0; i < parallelism; i++) {
        bsc += bscv[i];
    }
    problems.clear();
    for (auto it = bsc.pmpq.begin(); it != bsc.pmpq.end(); ++it) {
        problems.push_back(bsc.contents[it->second]);
    }
}

void run_main(  std::string infile,
                uint32_t beamwidth,
                uint32_t parallelism,
                uint32_t max_gliders_allowed,
                uint32_t max_extra_gens,
                uint32_t max_branching) {

    std::vector<SODThread> gsodv(parallelism);
    std::vector<ProblemState> problems(1);

    uint32_t depth=0;
    std::cerr << "Loading problem..." << std::flush;
    std::vector<apg::bitworld> vbw = load_file<4>(infile);
    UniqueEnsurer ue=gsodv[0].parse_pattern(vbw);
    ue.max_gliders_allowed = max_gliders_allowed;
    ue.max_extra_gens = max_extra_gens;
    ue.max_branching = max_branching;

    for (uint32_t i = 0; i < parallelism; i++) {
        gsodv[i].init_patterns(ue.origin, ue.reaction_allowed, ue.objects_forbidden);
    }

    start_gsod(gsodv, beamwidth, ue, problems);
    while (!(problems.empty())) {
        std::cerr << " (" << depth << ")" ;
        run1iter(gsodv, beamwidth, ue, problems, ++depth);
    }
}

} // namespace gsod


int main(int argc, char* argv[]) {

    bool incorrectly_called;
    if ((argc<4)||(argc>7)) {incorrectly_called = true;}
    for (int i = 2; i < argc; i++) {
        if (argv[i][0] == '-') { incorrectly_called = true; }
    }
    if (incorrectly_called) {
        std::cerr << " correct call gSoD <pattern file> <beamwidth> <parallelism> [<max_gliders_allowed>(2)] [<max_extra_gens>(1024)] [<max_branching>(16384)] " << std::endl;
        return 1;
    }

    uint32_t beamwidth = std::stoi(argv[2]);
    uint32_t parallelism = std::stoi(argv[3]);
    uint32_t max_gliders_allowed = 2;
    if (argc > 4) {
        max_gliders_allowed = std::stoi(argv[4]);
    }
    uint32_t max_extra_gens = 1024;
    if (argc > 5) {
        max_extra_gens = std::stoi(argv[5]);
    }
    uint32_t max_branching = 16384;
    if (argc > 6) {
        max_branching = std::stoi(argv[6]);
    }

    gsod::run_main(argv[1], beamwidth, parallelism, max_gliders_allowed, max_extra_gens, max_branching);

    return 0;
}
