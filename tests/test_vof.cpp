/**
 * @file test_vof.cpp
 * @brief Unit tests for VOF (Volume of Fluid) transport.
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include "grid.hpp"
#include "vof.hpp"
#include <cmath>

TEST_CASE("Total volume computation") {
    Grid grid(16, 8, 2.0, 1.0);
    grid.initialize_water(0.5);

    VOFTransport vof;
    double vol = vof.total_volume(grid);

    // With no SOLID border, 16 cols × 4 full rows × cell_area = 16*4*(2/16)*(1/8) = 1.0
    CHECK(vol == doctest::Approx(1.0).epsilon(0.02));
}

TEST_CASE("VOF advection with zero velocity preserves field") {
    Grid grid(16, 8, 2.0, 1.0);
    grid.initialize_water(0.5);
    grid.classify_cells();

    // All velocities are zero after initialize_water
    VOFTransport vof;
    double vol_before = vof.total_volume(grid);

    // Advect with zero velocity
    for (int step = 0; step < 10; ++step) {
        vof.advect(grid, 0.001);
    }

    double vol_after = vof.total_volume(grid);
    CHECK(vol_after == doctest::Approx(vol_before).epsilon(0.01));
}

TEST_CASE("VOF advection conserves mass approximately") {
    Grid grid(32, 16, 2.0, 1.0);
    grid.initialize_water(0.4);
    grid.classify_cells();

    // Uniform rightward flow in the interior; wall faces stay at zero.
    for (int i = 1; i < grid.NX; ++i)
        for (int j = 0; j < grid.NY; ++j)
            grid.U(i, j) = 0.5;

    VOFTransport vof;
    double vol_before = vof.total_volume(grid);

    // Advect several steps
    for (int step = 0; step < 20; ++step) {
        vof.advect(grid, 0.001);
    }

    double vol_after = vof.total_volume(grid);

    // Mass should be well conserved — the advect() pipeline includes mass correction,
    // so the remaining error comes only from clamping at [0,1] boundaries.
    double relative_change = std::abs(vol_after - vol_before) / vol_before;
    CHECK(relative_change < 0.1);
}

TEST_CASE("VOF values stay bounded [0,1]") {
    Grid grid(16, 8, 2.0, 1.0);
    grid.initialize_water(0.5);
    grid.classify_cells();

    for (int i = 1; i < grid.NX; ++i)
        for (int j = 0; j < grid.NY; ++j)
            grid.U(i, j) = 1.0;

    VOFTransport vof;
    for (int step = 0; step < 50; ++step) {
        vof.advect(grid, 0.002);
    }

    for (int i = 0; i < grid.NX; ++i) {
        for (int j = 0; j < grid.NY; ++j) {
            CHECK(grid.VOF(i, j) >= 0.0);
            CHECK(grid.VOF(i, j) <= 1.0);
        }
    }
}

TEST_CASE("VOF alternating sweep order") {
    Grid grid(16, 8, 2.0, 1.0);
    grid.initialize_water(0.5);
    grid.classify_cells();

    VOFTransport vof;

    // Just verify it doesn't crash with alternating sweeps
    for (int step = 0; step < 10; ++step) {
        vof.advect(grid, 0.001);
    }

    // Volume should still be reasonable
    double vol = vof.total_volume(grid);
    CHECK(vol > 0.0);
    CHECK(vol < grid.Lx * grid.Ly);
}
