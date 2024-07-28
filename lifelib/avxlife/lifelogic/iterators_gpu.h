#define ADVANCE_TILE_64(a, b, tmp, dst) {              \
    uint64_cu al = ROTL64(a, 1);                       \
    uint64_cu ar = ROTR64(a, 1);                       \
    uint64_cu xor2 = al ^ ar;                          \
    (dst)[threadIdx.x] = xor2 ^ a;                     \
    __syncthreads();                                   \
    uint64_cu uda = (dst)[u] & (dst)[d];               \
    uint64_cu udx = (dst)[u] ^ (dst)[d];               \
    uint64_cu sv = (xor2 & udx) | uda;                 \
    uint64_cu pt8 = xor2 ^ udx;                        \
    uint64_cu sh = al & ar;                            \
    (tmp)[threadIdx.x] = sh | (a & xor2);              \
    __syncthreads();                                   \
    uint64_cu pt4 = ((tmp)[u] ^ (tmp)[d]) ^ (sh ^ sv); \
    uint64_cu xc3 = ((tmp)[u] | (tmp)[d]) ^ (sh | sv); \
    b = pt4 & xc3 & (pt8 | a);                         \
}

#define SKIP_18_GENS 1

