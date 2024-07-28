#include <cpads/core.hpp>
#include "gtest/gtest.h"


TEST(Core, FloorSqrt) {

    for (uint64_t y = 0; y < 100; y++) {
        for (uint64_t x = y*y; x < (y+1)*(y+1); x++) {
            EXPECT_EQ(hh::floor_sqrt(x), y);
        }
    }

    EXPECT_EQ(hh::floor_sqrt((uint64_t) -1), 0xffffffffull);

}
