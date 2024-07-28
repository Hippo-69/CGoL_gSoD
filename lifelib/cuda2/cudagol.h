#pragma once

#include "basics.h"
#include "../classic/bpattern.h"


_DI_ uint4 load_bpattern(const apg::bpattern *ptr) {
    const uint4* u4ptr = ((const uint4*) ptr);
    uint4 result = u4ptr[threadIdx.x & 31];
    return result;
}

_DI_ void save_bpattern(apg::bpattern *ptr, uint4 slice) {
    uint4* u4ptr = ((uint4*) ptr);
    u4ptr[threadIdx.x & 31] = slice;
}

_DI_ uint4 load_hash(const unsigned int *hash_ptr) {

    unsigned int h = 0;

    int laneId = threadIdx.x & 31;

    if ((laneId >= 12) && (laneId < 20)) {
        h = hash_ptr[laneId - 12];
    }

    // elegantly perform necessary endianness change:
    uint4 result;
    result.x = (h & 0x0000ff00u) << 16;
    result.y = (h & 0x000000ffu);
    result.z = (h & 0xff000000u);
    result.w = (h & 0x00ff0000u) >> 16;

    return result;
}

_DI_ void hreflect_even(uint4 &b) {

    b.x = (b.x >> 8) | (b.y << 24);
    b.z = (b.z >> 8) | (b.w << 24);
    b.y = __brev(b.x);
    b.w = __brev(b.z);

}

_DI_ void hreflect_odd(uint4 &b) {

    hreflect_even(b);
    b.x <<= 1;
    b.z <<= 1;

}

_DI_ void vreflect_even(uint4 &b) {

    int laneId = threadIdx.x & 31;
    int src = (laneId < 16) ? (27 - laneId) : (laneId - 4);

    b.x = hh::shuffle_32(b.x, src);
    b.y = hh::shuffle_32(b.y, src);
    b.z = hh::shuffle_32(b.z, src);
    b.w = hh::shuffle_32(b.w, src);

    if (laneId < 16) {
        uint32_t c = b.x; b.x = b.z; b.z = c;
        c = b.y; b.y = b.w; b.w = c;
    }
}

_DI_ void vreflect_odd(uint4 &b) {

    int laneId = threadIdx.x & 31;
    int src1 = (laneId < 16) ? (28 - laneId) : (laneId - 4);
    int src2 = (laneId < 16) ? (27 - laneId) : (laneId - 4);

    b.x = hh::shuffle_32(b.x, src1);
    b.y = hh::shuffle_32(b.y, src1);
    b.z = hh::shuffle_32(b.z, src2);
    b.w = hh::shuffle_32(b.w, src2);

}

_DI_ void restrict_C2_4(uint4 &b) {

    int laneId = threadIdx.x & 31;
    if (laneId < 16) {
        b.x = 0; b.z = 0;
    } else {
        b.y = 0; b.w = 0;
    }
}

_DI_ void restrict_C2_2(uint4 &b) {

    int laneId = threadIdx.x & 31;
    if (laneId < 16) {
        b.x = 0; b.z = 0;
    } else {
        b.y &= 1;
        b.w &= 1;
    }
}

_DI_ void restrict_C2_1(uint4 &b) {

    int laneId = threadIdx.x & 31;
    if (laneId < 16) {
        b.x = 0; b.z = 0;
    } else {
        if (laneId > 16) { b.y &= 1; }
        b.w &= 1;
    }
}

/**
 * Advance a 64x64 torus by one generation using a single warp.
 *
 * Each thread holds a 64x2 rectangle of the torus in a uint4
 * datatype (representing a 2x2 array of 32x1 rectangles).
 */
_DI_ uint4 advance_torus(uint4 old) {

    // Precompute the lane indices of the threads representing
    // the 64x2 rectangles immediately above and below this one.
    int upperthread = (threadIdx.x + 31) & 31;
    int lowerthread = (threadIdx.x +  1) & 31;

    // 8 boolean ops + 8 funnel shifts:
    uint4 xor2; uint4 sh;
    {
        uint4 al; uint4 ar;
        al.x = (old.x << 1) | (old.y >> 31);
        al.y = (old.y << 1) | (old.x >> 31);
        al.z = (old.z << 1) | (old.w >> 31);
        al.w = (old.w << 1) | (old.z >> 31);
        ar.x = (old.x >> 1) | (old.y << 31);
        ar.y = (old.y >> 1) | (old.x << 31);
        ar.z = (old.z >> 1) | (old.w << 31);
        ar.w = (old.w >> 1) | (old.z << 31);

        xor2.x = al.x ^ ar.x;
        xor2.y = al.y ^ ar.y;
        xor2.z = al.z ^ ar.z;
        xor2.w = al.w ^ ar.w;
        sh.x = al.x & ar.x;
        sh.y = al.y & ar.y;
        sh.z = al.z & ar.z;
        sh.w = al.w & ar.w;
    }

    // 24 boolean ops + 4 shuffles:
    uint4 xor3; uint4 xor3prime; uint4 sv; uint4 pt8;
    xor3.x = xor2.x ^ old.x;
    xor3.y = xor2.y ^ old.y;
    xor3.z = xor2.z ^ old.z;
    xor3.w = xor2.w ^ old.w;
    xor3prime.x = hh::shuffle_32(xor3.x, lowerthread);
    xor3prime.y = hh::shuffle_32(xor3.y, lowerthread);
    xor3prime.z = hh::shuffle_32(xor3.z, upperthread);
    xor3prime.w = hh::shuffle_32(xor3.w, upperthread);
    {
        uint4 uda; uint4 udx;
        udx.x = xor3.z ^ xor3prime.z;
        udx.y = xor3.w ^ xor3prime.w;
        udx.z = xor3.x ^ xor3prime.x;
        udx.w = xor3.y ^ xor3prime.y;
        uda.x = xor3.z & xor3prime.z;
        uda.y = xor3.w & xor3prime.w;
        uda.z = xor3.x & xor3prime.x;
        uda.w = xor3.y & xor3prime.y;

        pt8.x = xor2.x ^ udx.x;
        pt8.y = xor2.y ^ udx.y;
        pt8.z = xor2.z ^ udx.z;
        pt8.w = xor2.w ^ udx.w;
        sv.x = (xor2.x & udx.x) | uda.x;
        sv.y = (xor2.y & udx.y) | uda.y;
        sv.z = (xor2.z & udx.z) | uda.z;
        sv.w = (xor2.w & udx.w) | uda.w;
    }

    // 32 boolean ops + 4 shuffles
    uint4 maj3; uint4 maj3prime; uint4 pt4; uint4 xc3;
    maj3.x = sh.x | (old.x & xor2.x);
    maj3.y = sh.y | (old.y & xor2.y);
    maj3.z = sh.z | (old.z & xor2.z);
    maj3.w = sh.w | (old.w & xor2.w);
    maj3prime.x = hh::shuffle_32(maj3.x, lowerthread);
    maj3prime.y = hh::shuffle_32(maj3.y, lowerthread);
    maj3prime.z = hh::shuffle_32(maj3.z, upperthread);
    maj3prime.w = hh::shuffle_32(maj3.w, upperthread);

    pt4.x = (maj3prime.z ^ maj3.z) ^ (sh.x ^ sv.x);
    pt4.y = (maj3prime.w ^ maj3.w) ^ (sh.y ^ sv.y);
    pt4.z = (maj3prime.x ^ maj3.x) ^ (sh.z ^ sv.z);
    pt4.w = (maj3prime.y ^ maj3.y) ^ (sh.w ^ sv.w);
    xc3.x = (maj3prime.z | maj3.z) ^ (sh.x | sv.x);
    xc3.y = (maj3prime.w | maj3.w) ^ (sh.y | sv.y);
    xc3.z = (maj3prime.x | maj3.x) ^ (sh.z | sv.z);
    xc3.w = (maj3prime.y | maj3.y) ^ (sh.w | sv.w);

    // 12 boolean ops
    uint4 newstate;
    newstate.x = pt4.x & xc3.x & (pt8.x | old.x);
    newstate.y = pt4.y & xc3.y & (pt8.y | old.y);
    newstate.z = pt4.z & xc3.z & (pt8.z | old.z);
    newstate.w = pt4.w & xc3.w & (pt8.w | old.w);

    #ifdef PEDESTRIAN_LIFE
    // correction for eighth birth transition
    newstate.x |= (maj3prime.z & maj3.z & xor3.z & xor3prime.z & sh.x &~ old.x);
    newstate.y |= (maj3prime.w & maj3.w & xor3.w & xor3prime.w & sh.y &~ old.y);
    newstate.z |= (maj3prime.x & maj3.x & xor3.x & xor3prime.x & sh.z &~ old.z);
    newstate.w |= (maj3prime.y & maj3.y & xor3.y & xor3prime.y & sh.w &~ old.w);
    #endif

    return newstate;

}


_DI_ void rotate_torus_inplace(uint4 &t, int rh, int rv) {

    if (rv & 63) {
        // translate vertically:
        uint4 d;
        d.x = (rv & 1) ? t.z : t.x;
        d.y = (rv & 1) ? t.w : t.y;
        d.z = (rv & 1) ? t.x : t.z;
        d.w = (rv & 1) ? t.y : t.w;
        int upperthread = (((rv)     >> 1) + threadIdx.x) & 31;
        int lowerthread = (((rv + 1) >> 1) + threadIdx.x) & 31;
        t.x = hh::shuffle_32(d.x, upperthread);
        t.y = hh::shuffle_32(d.y, upperthread);
        t.z = hh::shuffle_32(d.z, lowerthread);
        t.w = hh::shuffle_32(d.w, lowerthread);
    }

    if (rh & 63) {
        // translate horizontally:
        uint4 d;
        d.x = (rh & 32) ? t.y : t.x;
        d.y = (rh & 32) ? t.x : t.y;
        d.z = (rh & 32) ? t.w : t.z;
        d.w = (rh & 32) ? t.z : t.w;
        int sa = rh & 31;
        t.x = (d.x >> sa) | (d.y << (32 - sa));
        t.y = (d.y >> sa) | (d.x << (32 - sa));
        t.z = (d.z >> sa) | (d.w << (32 - sa));
        t.w = (d.w >> sa) | (d.z << (32 - sa));
    }
}


/**
 * Determine whether the pattern is too large to fit in a single tile.
 */
_DI_ bool escaped_bounding_box(uint4 &b) {

    unsigned int xz = b.x | b.z;
    unsigned int yw = b.y | b.w;

    int laneId = threadIdx.x & 31;
    bool hasMask = (laneId >= 3) && (laneId < 29);
    unsigned int xzmask = (hasMask) ? 0x0000003fu : 0xffffffffu;
    unsigned int ywmask = (hasMask) ? 0xfc000000u : 0xffffffffu;

    if (hh::ballot_32((xzmask & xz) | (ywmask & yw)) == 0) {
        // currently inside bounding box:
        return false;
    }

    // We're outside the 52x52 region, but it might be possible to
    // recentre the pattern such that it fits again. We examine the
    // horizontal and vertical extent:

    unsigned int active_threads = hh::ballot_32(xz | yw);

    int vt;
    if ((active_threads & 0xe0000007u) == 0) {
        vt = 0;
    } else if ((active_threads & 0xf0000003u) == 0) {
        vt = -2;
    } else if ((active_threads & 0xc000000fu) == 0) {
        vt = 2;
    } else if ((active_threads & 0xf8000001u) == 0) {
        vt = -4;
    } else if ((active_threads & 0x8000001fu) == 0) {
        vt = 4;
    } else if ((active_threads & 0xfc000000u) == 0) {
        vt = -6;
    } else if ((active_threads & 0x0000003fu) == 0) {
        vt = 6;
    } else {
        // pattern is too tall:
        return true;
    }

    xz |= hh::shuffle_xor_32(xz, 1);
    yw |= hh::shuffle_xor_32(yw, 1);
    xz |= hh::shuffle_xor_32(xz, 2);
    yw |= hh::shuffle_xor_32(yw, 2);
    xz |= hh::shuffle_xor_32(xz, 4);
    yw |= hh::shuffle_xor_32(yw, 4);
    xz |= hh::shuffle_xor_32(xz, 8);
    yw |= hh::shuffle_xor_32(yw, 8);
    xz |= hh::shuffle_xor_32(xz, 16);
    yw |= hh::shuffle_xor_32(yw, 16);

    int ht;
    if (xz == 0) {
        ht = 16;
    } else if (yw == 0) {
        ht = -16;
    } else {
        int rightborder = __clz(yw);
        int  leftborder = __clz(__brev(xz));

        if (leftborder + rightborder < 12) {
            // pattern is too wide:
            return true;
        }

        ht = (leftborder - rightborder) >> 1;
    }

    // recentre pattern:
    rotate_torus_inplace(b, ht, vt);
    return false;
}


/**
 * Advances a single tile and checks whether it has settled down.
 * If it settles down, return 0; otherwise, return the current
 * generation count.
 */
_DI_ int advance_initial(uint4 &a) {

    uint4 b = a;

    // Run for 30 generations without any checks. It is impossible for
    // the pattern to escape the 52x52 bounding box in this time period,
    // and only 1.6% of soups stabilise this quickly so we do not give
    // up much by skipping the stabilisation checks:
    #ifdef C1_SYMMETRY
    #define initial_gens 30
    #else
    #define initial_gens 10
    #endif

    #pragma unroll 2
    for (int i = 0; i < initial_gens; i++) {
        b = advance_torus(b);
    }

    int age = initial_gens - 6;

    // Run until generation 600, checking for bounding box escapement
    // and 2-periodicity. This loop should fit comfortably in the GPU
    // icache. About 36% of soups stabilise in this loop.
    do {

        age += 6;
        a = b;

        b = advance_torus(b);
        b = advance_torus(b);

        if (hh::ballot_32((a.x ^ b.x) | (a.y ^ b.y) | (a.z ^ b.z) | (a.w ^ b.w)) == 0) {
            // periodic with period 2:
            return 0;
        }

        b = advance_torus(b);
        b = advance_torus(b);
        b = advance_torus(b);
        b = advance_torus(b);

        if (escaped_bounding_box(b)) {
            break;
        }
    } while (age < 600);

    // Conveniently, a and b differ by 6 generations, so we can check
    // for soups which stabilise with the production of a pulsar.
    if (hh::ballot_32((a.x ^ b.x) | (a.y ^ b.y) | (a.z ^ b.z) | (a.w ^ b.w)) == 0) {
        // periodic with period 6:
        return 0;
    }

    // ================================================================
    // |                  CHECK FOR ESCAPING GLIDERS                  |
    // ================================================================

    // Take intersection of two snapshots separated by 8 generations:
    uint4 c;
    b = advance_torus(b);
    b = advance_torus(b);
    c.x = a.x & b.x;
    c.y = a.y & b.y;
    c.z = a.z & b.z;
    c.w = a.w & b.w;

    {
        uint4 d;
        d = advance_torus(c);
        d = advance_torus(d);

        if (hh::ballot_32((c.x ^ d.x) | (c.y ^ d.y) | (c.z ^ d.z) | (c.w ^ d.w)) != 0) {
            // core does not have period 2:
            return age;
        }
    }

    {
        int totpop = __popc(c.x ^ a.x) + __popc(c.y ^ a.y) + __popc(c.z ^ a.z) + __popc(c.w ^ a.w);
        totpop     = totpop << 16;
        totpop    += __popc(c.x ^ b.x) + __popc(c.y ^ b.y) + __popc(c.z ^ b.z) + __popc(c.w ^ b.w);

        totpop += hh::shuffle_xor_32(totpop, 1);
        totpop += hh::shuffle_xor_32(totpop, 2);
        totpop += hh::shuffle_xor_32(totpop, 4);
        totpop += hh::shuffle_xor_32(totpop, 8);
        totpop += hh::shuffle_xor_32(totpop, 16);

        if (totpop != 0x00050005u) {
            // diff has wrong population to be a glider:
            return age;
        }
    }

    {
        // Advance by another 6 generations to ensure glider has fully escaped:
        for (int i = 0; i < 4; i++) {
            b = advance_torus(b);
        }

        int laneId = threadIdx.x & 31;
        bool hasMask = (laneId >= 3) && (laneId < 29);
        unsigned int xzmask = (hasMask) ? 0xffffffc0u : 0;
        unsigned int ywmask = (hasMask) ? 0x03ffffffu : 0;

        int totpop = 0;

        for (int i = 0; i < 4; i++) {
            b = advance_torus(b);
            b = advance_torus(b);
            if (hh::ballot_32(((b.x & xzmask) ^ c.x) | ((b.y & ywmask) ^ c.y) | ((b.z & xzmask) ^ c.z) | ((b.w & ywmask) ^ c.w))) { return false; }
            totpop     = totpop << 8;
            totpop    += __popc(c.x ^ b.x) + __popc(c.y ^ b.y) + __popc(c.z ^ b.z) + __popc(c.w ^ b.w);
        }

        totpop += hh::shuffle_xor_32(totpop, 1);
        totpop += hh::shuffle_xor_32(totpop, 2);
        totpop += hh::shuffle_xor_32(totpop, 4);
        totpop += hh::shuffle_xor_32(totpop, 8);
        totpop += hh::shuffle_xor_32(totpop, 16);

        return ((totpop == 0x05050505u) ? 0 : age);
    }
}


_DI_ int advance_tile_inplace(uint4 &a, int rh, int rv) {

    uint4 b = advance_torus(a);
    uint4 c = advance_torus(b);

    b = advance_torus(c);
    b = advance_torus(b);
    b = advance_torus(b);
    b = advance_torus(b);

    unsigned int leftdiff;
    unsigned int rightdiff;

    {
        int laneId = threadIdx.x & 31;
        bool hasMask = (laneId >= 3) && (laneId < 29);
        unsigned int xzmask = (hasMask) ? 0xffffffc0u : 0;
        unsigned int ywmask = (hasMask) ? 0x03ffffffu : 0;

        uint4 xordiff;
        xordiff.x = (a.x ^ b.x) & xzmask;
        xordiff.y = (a.y ^ b.y) & ywmask;
        xordiff.z = (a.z ^ b.z) & xzmask;
        xordiff.w = (a.w ^ b.w) & ywmask;

        // update central 52-by-52 tile:
        a.x ^= xordiff.x;
        a.y ^= xordiff.y;
        a.z ^= xordiff.z;
        a.w ^= xordiff.w;

        xordiff.x |= xordiff.z;
        xordiff.y |= xordiff.w;
        leftdiff  = (xordiff.x >> 6) | (xordiff.y << 26);
        rightdiff = (xordiff.y << 6) | (xordiff.x >> 26);
    }

    int flags = 64;
    flags |= (hh::ballot_32(rightdiff & 0xfc000000u) ? 1 : 0);
    flags |= (hh::ballot_32( leftdiff & 0x0000003fu) ? 8 : 0);
    leftdiff  = hh::ballot_32(leftdiff);
    rightdiff = hh::ballot_32(rightdiff);

    if ((leftdiff | rightdiff) == 0) { return 0; }

    flags |= ((rightdiff & 0x00000038u) ? 2  : 0);
    flags |= (( leftdiff & 0x00000038u) ? 4  : 0);
    flags |= (( leftdiff & 0x1c000000u) ? 16 : 0);
    flags |= ((rightdiff & 0x1c000000u) ? 32 : 0);

    // Escaping glider/*WSS detection:
    if ((rh & 63) || (rv & 63)) {

        unsigned int xz = b.x | b.z;
        unsigned int yw = b.y | b.w;

        if (hh::ballot_32(xz | yw)) {

            int laneId = threadIdx.x & 31;
            bool hasMask = (laneId >= 3) && (laneId < 29);
            unsigned int xzmask = (hasMask) ? 0x0000003fu : 0xffffffffu;
            unsigned int ywmask = (hasMask) ? 0xfc000000u : 0xffffffffu;

            if (hh::ballot_32((xz & xzmask) | (yw & ywmask)) == 0) {

                rotate_torus_inplace(c, rh, rv);

                if (hh::ballot_32((c.x ^ b.x) | (c.y ^ b.y) | (c.z ^ b.z) | (c.w ^ b.w)) == 0) {

                    int totpop = __popc(b.x) + __popc(b.y) + __popc(b.z) + __popc(b.w);
                    totpop += hh::shuffle_xor_32(totpop, 1);
                    totpop += hh::shuffle_xor_32(totpop, 2);
                    totpop += hh::shuffle_xor_32(totpop, 4);
                    totpop += hh::shuffle_xor_32(totpop, 8);
                    totpop += hh::shuffle_xor_32(totpop, 16);

                    if (totpop <= 17) {
                        // either a single LWSS, MWSS, or HWSS, or between
                        // 1 and 3 gliders. Remove them from the universe:
                        a.x &= xzmask;
                        a.y &= ywmask;
                        a.z &= xzmask;
                        a.w &= ywmask;
                    }
                }
            }
        }
    }

    return flags;

}


__global__ void exhaustFirstTile(uint32_cu *hashes, uint32_cu *interesting, uint4 *output) {

    int pos = (threadIdx.x + blockIdx.x * blockDim.x) >> 5;

    uint32_t *thishash = hashes + (pos << 3);

    uint4 a = load_hash(thishash);

    #ifdef HREFLECT_EVEN
    hreflect_even(a);
    #endif
    #ifdef HREFLECT_ODD
    hreflect_odd(a);
    #endif
    #ifdef VREFLECT_EVEN
    vreflect_even(a);
    #endif
    #ifdef VREFLECT_ODD
    vreflect_odd(a);
    #endif
    #ifdef RESTRICT_C2_1
    restrict_C2_1(a);
    #endif
    #ifdef RESTRICT_C2_2
    restrict_C2_2(a);
    #endif
    #ifdef RESTRICT_C2_4
    restrict_C2_4(a);
    #endif

    int gencount = advance_initial(a);

    interesting[pos] = gencount;

    if (gencount && output) {
        int laneId = threadIdx.x & 31;
        output[(pos << 12) + laneId + 32] = a;

        uint4 b; b.x = 0; b.y = 0; b.z = 0; b.w = 0;

        if (laneId == 0) { b.y = 127; }

        output[(pos << 12) + laneId] = b;
    }
}

#include "emt.h"
#define ULTIMATE_EMT 1
#include "emt.h"

