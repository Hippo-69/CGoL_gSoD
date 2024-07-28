#pragma once
#include "../mu.hpp"

namespace hh {

template<size_t E, typename T, size_t N, bool alternate=false>
_DI_ void warp_bitonic_sort(hh::vec<T, N> &x) {

    bool descending = false;
    if constexpr (alternate) {
        descending = (threadIdx.x >> E) & 1;
    }

    if constexpr (E > 0) {

        warp_bitonic_sort<E-1, T, N, true>(x);

        for (int i = 0; i < ((int) E); i++) {
            int m = 1 << (E-1-i);

            #pragma unroll
            for (size_t j = 0; j < N; j++) {
                T other = shuffle_xor_32(x[j], m);
                compare_and_swap(x[j], other, descending ^ ((bool) (threadIdx.x & m)));
            }
        }
    }

    x.sort(descending);
}

template<size_t E, typename T, size_t N>
_DI_ bool warp_memcmp_leq(const hh::vec<T, N> &x, const hh::vec<T, N> &y) {

    if constexpr (E == 0) {
        return x <= y;
    } else {
        uint32_t mask = ((uint32_t) -1);
        if constexpr (E <= 4) { mask &= ((threadIdx.x & 16) ? 0xffff0000u : 0x0000ffffu); }
        if constexpr (E <= 3) { mask &= ((threadIdx.x &  8) ? 0xff00ff00u : 0x00ff00ffu); }
        if constexpr (E <= 2) { mask &= ((threadIdx.x &  4) ? 0xf0f0f0f0u : 0x0f0f0f0fu); }
        if constexpr (E <= 1) { mask &= ((threadIdx.x &  2) ? 0xccccccccu : 0x33333333u); }
        uint32_t lt_mask = ballot_32(x < y) & mask;
        uint32_t gt_mask = ballot_32(x > y) & mask;
        return brev32(gt_mask) <= brev32(lt_mask);
    }
}

template<size_t E, typename T, size_t N>
_DI_ void warp_memcmp_min(const hh::vec<T, N> &src, hh::vec<T, N> &dst) {

    if (warp_memcmp_leq<E>(src, dst)) {
        #pragma unroll
        for (size_t i = 0; i < N; i++) {
            dst[i] = src[i];
        }
    }
}

}
