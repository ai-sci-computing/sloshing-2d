/**
 * @file test_grid.cpp
 * @brief Unit tests for the MAC staggered grid data structure.
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include "grid.hpp"
#include <cmath>

TEST_CASE("Grid construction and dimensions") {
    Grid grid(32, 16, 2.0, 1.0);

    CHECK(grid.NX == 32);
    CHECK(grid.NY == 16);
    CHECK(grid.dx == doctest::Approx(2.0 / 32));
    CHECK(grid.dy == doctest::Approx(1.0 / 16));
    CHECK(grid.Lx == doctest::Approx(2.0));
    CHECK(grid.Ly == doctest::Approx(1.0));
}

TEST_CASE("Grid array sizes") {
    Grid grid(32, 16, 2.0, 1.0);

    // u has (NX+1) * NY elements
    CHECK(grid.u.size() == 33 * 16);
    // v has NX * (NY+1) elements
    CHECK(grid.v.size() == 32 * 17);
    // Pressure, VOF, cell_type have NX * NY elements
    CHECK(grid.p.size() == 32 * 16);
    CHECK(grid.vof.size() == 32 * 16);
    CHECK(grid.cell_type.size() == 32 * 16);
}

TEST_CASE("Grid field accessors") {
    Grid grid(8, 4, 1.0, 0.5);

    // Write and read back u-velocity
    grid.U(3, 2) = 1.5;
    CHECK(grid.U(3, 2) == doctest::Approx(1.5));

    // Write and read back v-velocity
    grid.V(5, 3) = -2.0;
    CHECK(grid.V(5, 3) == doctest::Approx(-2.0));

    // Write and read back pressure
    grid.P(4, 2) = 100.0;
    CHECK(grid.P(4, 2) == doctest::Approx(100.0));

    // Write and read back VOF
    grid.VOF(6, 1) = 0.75;
    CHECK(grid.VOF(6, 1) == doctest::Approx(0.75));
}

TEST_CASE("Cell classification") {
    Grid grid(8, 4, 1.0, 0.5);

    // No obstacle flags set → classification is purely VOF-based.
    grid.VOF(3, 2) = 0.8;     // FLUID
    grid.VOF(4, 2) = 0.005;   // below threshold → EMPTY
    grid.VOF(0, 0) = 1.0;     // border cell: no longer forced SOLID
    grid.VOF(7, 3) = 0.0;     // EMPTY

    grid.classify_cells();

    CHECK(grid.Type(3, 2) == CellType::FLUID);
    CHECK(grid.Type(4, 2) == CellType::EMPTY);
    CHECK(grid.Type(0, 0) == CellType::FLUID);
    CHECK(grid.Type(7, 3) == CellType::EMPTY);

    // Internal obstacle flag forces SOLID regardless of VOF.
    grid.is_obstacle[grid.idx(2, 1)] = 1;
    grid.VOF(2, 1) = 0.9;
    grid.classify_cells();
    CHECK(grid.Type(2, 1) == CellType::SOLID);
}

TEST_CASE("Water initialization") {
    Grid grid(16, 8, 2.0, 1.0);
    grid.initialize_water(0.5); // Fill bottom half

    // No SOLID border any more: bottom-row cells are submerged.
    // Water height = 0.5 m, dy = 1/8 = 0.125 m.
    // j=0..3 fully submerged; j=4..7 empty.
    CHECK(grid.VOF(8, 0) == doctest::Approx(1.0));
    CHECK(grid.VOF(0, 2) == doctest::Approx(1.0));
    CHECK(grid.VOF(8, 3) == doctest::Approx(1.0));
    CHECK(grid.VOF(8, 4) == doctest::Approx(0.0));
    CHECK(grid.VOF(8, 6) == doctest::Approx(0.0));

    for (auto val : grid.u) CHECK(val == doctest::Approx(0.0));
    for (auto val : grid.v) CHECK(val == doctest::Approx(0.0));
}

TEST_CASE("Velocity interpolation at cell centers") {
    Grid grid(8, 4, 1.0, 0.5);

    // Set uniform u=2.0, v=0 flow
    for (int i = 0; i <= grid.NX; ++i)
        for (int j = 0; j < grid.NY; ++j)
            grid.U(i, j) = 2.0;

    // Interpolation at cell center should return ~2.0 for u
    double cx = 4.5 * grid.dx; // Center of cell (4, 2)
    double cy = 2.5 * grid.dy;
    auto [u, v] = grid.interpolate_velocity(cx, cy);
    CHECK(u == doctest::Approx(2.0).epsilon(0.01));
}
