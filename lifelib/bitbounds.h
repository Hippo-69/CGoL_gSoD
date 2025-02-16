#pragma once
#include <stdint.h>

namespace apg {

    uint64_t hshear(uint64_t x) {

        uint64_t y = (x & 0xff00ff00ff00ff00ull) | ((x & 0x007f007f007f007full) << 1) | ((x & 0x0080008000800080ull) >> 7);
        y = (y & 0xffff0000ffff0000ull) | ((y & 0x00003f3f00003f3full) << 2) | ((y & 0x0000c0c00000c0c0ull) >> 6);
        y = (y & 0xffffffff00000000ull) | ((y & 0x000000000f0f0f0full) << 4) | ((y & 0x00000000f0f0f0f0ull) >> 4);
        return y;

    }

    uint64_t vshear(uint64_t x) {

        uint64_t y = (x & 0xaaaaaaaaaaaaaaaaull) | ((x & 0x0055555555555555ull) << 8) | ((x & 0x5500000000000000ull) >> 56);
        y = (y & 0xccccccccccccccccull) | ((y & 0x0000333333333333ull) << 16) | ((y & 0x3333000000000000ull) >> 48);
        y = (y & 0xf0f0f0f0f0f0f0f0ull) | ((y & 0x000000000f0f0f0full) << 32) | ((y & 0x0f0f0f0f00000000ull) >> 32);
        return y;

    }

    uint64_t mxor(uint64_t a, uint64_t b) {

        // Shear matrices to obtain:
        //    c_(i,j) := a_(i,i+j)
        //    d_(i,j) := b_(i+j,j)
        uint64_t c = hshear(a);
        uint64_t d = vshear(b);
        uint64_t e = c & d;

        // -O3 should unroll this loop:
        for (int k = 1; k < 8; k++) {
            c = ((c & 0x7f7f7f7f7f7f7f7full) << 1) | ((c & 0x8080808080808080ull) >> 7);
            d = (d << 8) | (d >> 56);
            e ^= (c & d);
        }

        // e_(i,j) = xor_{k=0}^{7} a_(i,k) b_(k,j)
        return e;

    }

    uint64_t mor(uint64_t a, uint64_t b) {

        // Shear matrices to obtain:
        //    c_(i,j) := a_(i,i+j)
        //    d_(i,j) := b_(i+j,j)
        uint64_t c = hshear(a);
        uint64_t d = vshear(b);
        uint64_t e = c & d;

        // -O3 should unroll this loop:
        for (int k = 1; k < 8; k++) {
            c = ((c & 0x7f7f7f7f7f7f7f7full) << 1) | ((c & 0x8080808080808080ull) >> 7);
            d = (d << 8) | (d >> 56);
            e |= (c & d);
        }

        // e_(i,j) = or_{k=0}^{7} a_(i,k) b_(k,j)
        return e;

    }

    uint64_t upset(uint64_t x) {
        uint64_t y = x;
        y |= ((y & 0x5555555555555555) << 1);
        y |= ((y & 0x3333333333333333) << 2);
        y |= ((y & 0x0f0f0f0f0f0f0f0f) << 4);
        y |= ((y & 0x00ff00ff00ff00ff) << 8);
        y |= ((y & 0x0000ffff0000ffff) << 16);
        y |= ((y & 0x00000000ffffffff) << 32);
        return y;
    }

    void uint64_convolve2(uint64_t a, uint64_t b, uint64_t *out, bool exclusive) {
        // Convolution (works best when b is sparse):
        uint64_t brem = b;
        while (brem) {
            uint64_t c = (brem & (-brem)); // extract a single bit from brem
            brem ^= c; // remove bit
            uint64_t tzc = __builtin_ctzll(c); // determine shifts
            uint64_t xs = tzc & 7;
            uint64_t ys = tzc & 56;
            uint64_t bitmask = (0x0101010101010101ull << xs) - 0x0101010101010101ull;
            uint64_t right = (a >> (8 - xs)) & bitmask;
            uint64_t left = (a << xs) & (~bitmask);
            if (exclusive) {
                out[0] ^= (left << ys);
                out[1] ^= (right << ys);
                if (ys) {
                    out[2] ^= (left >> (64 - ys));
                    out[3] ^= (right >> (64 - ys));
                }
            } else {
                out[0] |= (left << ys);
                out[1] |= (right << ys);
                if (ys) {
                    out[2] |= (left >> (64 - ys));
                    out[3] |= (right >> (64 - ys));
                }
            }
        }
    }

    void uint64_convolve(uint64_t a, uint64_t b, uint64_t *out, bool exclusive) {
        // Convolve two 8-by-8 squares to produce a 16-by-16 result:
        if (__builtin_popcountll(a) > __builtin_popcountll(b)) {
            uint64_convolve2(a, b, out, exclusive);
        } else {
            uint64_convolve2(b, a, out, exclusive);
        }
    }

    uint64_t uint64_bottom(uint64_t tile) {
        uint64_t dy = 0;
        if (tile & 0xffffffff00000000ull) {
            dy = 4;
            if (tile & 0xffff000000000000ull) {
                dy = 6;
                if (tile & 0xff00000000000000ull) {
                    dy = 7;
                }
            } else if (tile & 0x0000ff0000000000ull) {
                dy = 5;
            }
        } else if (tile & 0x00000000ffff0000ull) {
            dy = 2;
            if (tile & 0x00000000ff000000ull) {
                dy = 3;
            }
        } else if (tile & 0x000000000000ff00ull) {
            dy = 1;
        }
        return dy;
    }

    uint64_t uint64_bl(uint64_t tile) {
        uint64_t dz = 14;
        if (tile & 0xff7f3f1f0f070301ull) {
            dz = 7;
            if (tile & 0x0f07030100000000ull) {
                dz = 3;
                if (tile & 0x0301000000000000ull) {
                    dz = 1;
                    if (tile & 0x0100000000000000ull) {
                        dz = 0;
                    }
                } else if (tile & 0x0402010000000000ull) {
                    dz = 2;
                }
            } else if (tile & 0x30180c0603010000ull) {
                dz = 5;
                if (tile & 0x1008040201000000ull) {
                    dz = 4;
                }
            } else if (tile & 0x4020100804020100ull) {
                dz = 6;
            }
        } else if (tile & 0x0080c0e0f0783c1eull) {
            dz = 11;
            if (tile & 0x0080c06030180c06ull) {
                dz = 9;
                if (tile & 0x0080402010080402ull) {
                    dz = 8;
                }
            } else if (tile & 0x0000008040201008ull) {
                dz = 10;
            }
        } else if (tile & 0x000000000080c060ull) {
            dz = 13;
            if (tile & 0x0000000000804020ull) {
                dz = 12;
            }
        }
        return dz;
    }

    uint64_t uint64_tr(uint64_t tile) {
        uint64_t dz = 0;
        if (tile & 0x80c0e0f0f8fcfeffull) {
            dz = 7;
            if (tile & 0x0000000080c0e0f0ull) {
                dz = 11;
                if (tile & 0x00000000000080c0ull) {
                    dz = 13;
                    if (tile & 0x0000000000000080ull) {
                        dz = 14;
                    }
                } else if (tile & 0x0000000000804020ull) {
                    dz = 12;
                }
            } else if (tile & 0x000080c06030180cull) {
                dz = 9;
                if (tile & 0x0000008040201008ull) {
                    dz = 10;
                }
            } else if (tile & 0x0080402010080402ull) {
                dz = 8;
            }
        } else if (tile & 0x783c1e0f07030100ull) {
            dz = 3;
            if (tile & 0x6030180c06030100ull) {
                dz = 5;
                if (tile & 0x4020100804020100ull) {
                    dz = 6;
                }
            } else if (tile & 0x1008040201000000ull) {
                dz = 4;
            }
        } else if (tile & 0x0603010000000000ull) {
            dz = 1;
            if (tile & 0x0402010000000000ull) {
                dz = 2;
            }
        }
        return dz;
    }

    uint64_t uint64_tl(uint64_t tile) {
        uint64_t dz = 14;
        if (tile & 0x0103070f1f3f7fffull) {
            dz = 7;
            if (tile & 0x000000000103070full) {
                dz = 3;
                if (tile & 0x0000000000000103ull) {
                    dz = 1;
                    if (tile & 0x0000000000000001ull) {
                        dz = 0;
                    }
                } else if (tile & 0x0000000000010204ull) {
                    dz = 2;
                }
            } else if (tile & 0x00000103060c1830ull) {
                dz = 5;
                if (tile & 0x0000000102040810ull) {
                    dz = 4;
                }
            } else if (tile & 0x0001020408102040ull) {
                dz = 6;
            }
        } else if (tile & 0x1e3c78f0e0c08000ull) {
            dz = 11;
            if (tile & 0x060c183060c08000ull) {
                dz = 9;
                if (tile & 0x0204081020408000ull) {
                    dz = 8;
                }
            } else if (tile & 0x0810204080000000ull) {
                dz = 10;
            }
        } else if (tile & 0x60c0800000000000ull) {
            dz = 13;
            if (tile & 0x2040800000000000ull) {
                dz = 12;
            }
        }
        return dz;
    }

    uint64_t uint64_br(uint64_t tile) {
        uint64_t dz = 0;
        if (tile & 0xfffefcf8f0e0c080ull) {
            dz = 7;
            if (tile & 0xf0e0c08000000000ull) {
                dz = 11;
                if (tile & 0xc080000000000000ull) {
                    dz = 13;
                    if (tile & 0x8000000000000000ull) {
                        dz = 14;
                    }
                } else if (tile & 0x2040800000000000ull) {
                    dz = 12;
                }
            } else if (tile & 0x0c183060c0800000ull) {
                dz = 9;
                if (tile & 0x0810204080000000ull) {
                    dz = 10;
                }
            } else if (tile & 0x0204081020408000ull) {
                dz = 8;
            }
        } else if (tile & 0x000103070f1e3c78ull) {
            dz = 3;
            if (tile & 0x000103060c183060ull) {
                dz = 5;
                if (tile & 0x0001020408102040ull) {
                    dz = 6;
                }
            } else if (tile & 0x0000000102040810ull) {
                dz = 4;
            }
        } else if (tile & 0x0000000000010306ull) {
            dz = 1;
            if (tile & 0x0000000000010204ull) {
                dz = 2;
            }
        }
        return dz;
    }

    uint64_t uint64_top(uint64_t tile) {
        uint64_t dy = 7;
        if (tile & 0x00000000ffffffffull) {
            dy = 3;
            if (tile & 0x000000000000ffffull) {
                dy = 1;
                if (tile & 0x00000000000000ffull) {
                    dy = 0;
                }
            } else if (tile & 0x0000000000ff0000ull) {
                dy = 2;
            }
        } else if (tile & 0x0000ffff00000000ull) {
            dy = 5;
            if (tile & 0x000000ff00000000ull) {
                dy = 4;
            }
        } else if (tile & 0x00ff000000000000ull) {
            dy = 6;
        }
        return dy;
    }

    uint64_t uint64_right(uint64_t tile) {
        uint64_t dx = 0;
        if (tile & 0xf0f0f0f0f0f0f0f0ull) {
            dx = 4;
            if (tile & 0xc0c0c0c0c0c0c0c0ull) {
                dx = 6;
                if (tile & 0x8080808080808080ull) {
                    dx = 7;
                }
            } else if (tile & 0x2020202020202020ull) {
                dx = 5;
            }
        } else if (tile & 0x0c0c0c0c0c0c0c0cull) {
            dx = 2;
            if (tile & 0x0808080808080808ull) {
                dx = 3;
            }
        } else if (tile & 0x0202020202020202ull) {
            dx = 1;
        }
        return dx;
    }

    uint64_t uint64_left(uint64_t tile) {
        uint64_t dx = 7;
        if (tile & 0x0f0f0f0f0f0f0f0full) {
            dx = 3;
            if (tile & 0x0303030303030303ull) {
                dx = 1;
                if (tile & 0x0101010101010101ull) {
                    dx = 0;
                }
            } else if (tile & 0x0404040404040404ull) {
                dx = 2;
            }
        } else if (tile & 0x3030303030303030ull) {
            dx = 5;
            if (tile & 0x1010101010101010ull) {
                dx = 4;
            }
        } else if (tile & 0x4040404040404040ull) {
            dx = 6;
        }
        return dx;
    }
}
