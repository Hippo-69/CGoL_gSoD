#pragma once
#include <string>
#include <vector>
#include "iterators_b3s23.h"
namespace apg {

    std::string get_zoi(std::string rule) {
        (void) rule;
        return "99";
    }

    std::string get_syms(std::string rule) {
        (void) rule;
        return "99";
    }

    int rule2int(std::string rule) {
        if (rule == "b3s23") { return 0; }
        return -1;
    }

    std::vector<std::string> get_all_rules() {
        std::vector<std::string> x;
        x.emplace_back("b3s23");
        return x;
    }

    int uli_get_family(int rule) {
        switch (rule) {
            case 0 : return 0;
        }
        return 0;
    }

    uint64_t uli_valid_mantissa(int rule) {
        switch (rule) {
            case 0 : return 511;
        }
        return 3;
    }

    int iterate_var_leaf(int rule, int n, uint64_t * inleaves, uint64_t * outleaf) {
        switch(rule) {
            case 0 :
                return b3s23::iterate_var_leaf(n, inleaves, outleaf);
        }
        return -1;
    }

    int iterate_var_32_28(int rule, uint32_t* d, uint32_t * diffs) {
        switch(rule) {
            case 0 :
                return b3s23::iterate_var_32_28(d, diffs);
        }
        return -1;
    }



#ifdef __AVX512F__

    int iterate_var_48_28(int rule, uint32_t* d, uint32_t * diffs) {
        switch(rule) {
            case 0 :
                return b3s23::iterate_var_48_28(d, diffs);
        }
        return -1;
    }



#endif

    int iterate_var_leaf(int rule, int n, uint64_t * inleaves, uint64_t * hleaves, uint64_t * outleaf) {
        switch(rule) {
            case 0 :
                return b3s23::iterate_var_leaf(n, inleaves, hleaves, outleaf);
        }
        return -1;
    }

    int iterate_var_32_28(int rule, uint32_t* d, uint32_t* h, uint32_t * diffs) {
        switch(rule) {
            case 0 :
                return b3s23::iterate_var_32_28(d, h, diffs);
        }
        return -1;
    }



#ifdef __AVX512F__

    int iterate_var_48_28(int rule, uint32_t* d, uint32_t* h, uint32_t * diffs) {
        switch(rule) {
            case 0 :
                return b3s23::iterate_var_48_28(d, h, diffs);
        }
        return -1;
    }



#endif

    int iterate_var_leaf(int rule, int n, uint64_t * inleaves, uint64_t * hleaves, uint64_t * jleaves, uint64_t * outleaf) {
        switch(rule) {
            case 0 :
                return b3s23::iterate_var_leaf(n, inleaves, hleaves, jleaves, outleaf);
        }
        return -1;
    }

    int iterate_var_32_28(int rule, uint32_t* d, uint32_t* h, uint32_t* j, uint32_t * diffs) {
        switch(rule) {
            case 0 :
                return b3s23::iterate_var_32_28(d, h, j, diffs);
        }
        return -1;
    }



#ifdef __AVX512F__

    int iterate_var_48_28(int rule, uint32_t* d, uint32_t* h, uint32_t* j, uint32_t * diffs) {
        switch(rule) {
            case 0 :
                return b3s23::iterate_var_48_28(d, h, j, diffs);
        }
        return -1;
    }



#endif

    void iterate_var_grid(int rule, uint8_t *ingrid, uint8_t *outgrid) {
        switch(rule) {
        }
        (void) ingrid; (void) outgrid;
    }

}
