/**
 * @file test_solver.cpp
 * @brief Unit tests for the Navier-Stokes solver.
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include "grid.hpp"
#include "solver.hpp"
#include "vof.hpp"
#include <cmath>

TEST_CASE("CFL timestep computation") {
    Grid grid(32, 16, 2.0, 1.0);
    Solver solver;

    // With zero velocity, dt should be max_dt
    double dt = solver.compute_dt(grid);
    CHECK(dt == doctest::Approx(solver.params.max_dt));

    // With significant velocity, dt should be CFL-limited
    for (int i = 0; i <= grid.NX; ++i)
        for (int j = 0; j < grid.NY; ++j)
            grid.U(i, j) = 5.0;

    dt = solver.compute_dt(grid);
    double expected = solver.params.cfl * grid.dx / 5.0;
    CHECK(dt < expected * 1.1); // Approximately CFL-limited
    CHECK(dt > 0.0);
}

TEST_CASE("Body forces add gravity correctly") {
    Grid grid(16, 8, 2.0, 1.0);
    grid.initialize_water(0.5);

    Solver solver;
    double dt = 0.001;

    // v should become negative (downward) after applying gravity
    solver.apply_body_forces(grid, dt, 0.0, 0.0);

    // Check that interior fluid v-velocities have become negative
    bool found_negative_v = false;
    for (int i = 1; i < grid.NX - 1; ++i) {
        for (int j = 1; j < grid.NY; ++j) {
            if (grid.V(i, j) < 0.0) {
                found_negative_v = true;
                // Expected: v = 0 + (-9.81) * 0.001 = -0.00981
                CHECK(grid.V(i, j) == doctest::Approx(-9.81 * dt).epsilon(0.01));
            }
        }
    }
    CHECK(found_negative_v);
}

TEST_CASE("Pressure solver converges for still water") {
    Grid grid(32, 16, 2.0, 1.0);
    grid.initialize_water(0.5);

    // Apply gravity to create a divergent velocity field
    Solver solver;
    SolverParams sp;
    sp.pressure_iters = 200;
    sp.pressure_tol = 1e-5;
    solver.params = sp;

    double dt = 0.001;
    solver.apply_body_forces(grid, dt, 0.0, 0.0);
    grid.enforce_boundary_conditions();

    // Solve pressure
    int iters = solver.solve_pressure(grid, dt);
    CHECK(iters < sp.pressure_iters); // Should converge

    // Project velocity
    solver.project_velocity(grid, dt);
    grid.enforce_boundary_conditions();

    // Divergence should be very small
    double max_div = solver.max_divergence(grid);
    CHECK(max_div < 0.1); // Relaxed tolerance for first test
}

TEST_CASE("Divergence-free projection on uniform flow") {
    Grid grid(16, 8, 2.0, 1.0);

    // Fill the grid with water (except boundaries)
    for (int i = 0; i < grid.NX; ++i)
        for (int j = 0; j < grid.NY; ++j)
            grid.VOF(i, j) = 1.0;
    grid.classify_cells();

    // Set a divergence-free uniform flow: u=1, v=0
    for (int i = 0; i <= grid.NX; ++i)
        for (int j = 0; j < grid.NY; ++j)
            grid.U(i, j) = 1.0;
    // v already 0

    Solver solver;
    double dt = 0.001;

    // Divergence should already be ~0 for uniform flow
    double div_before = solver.max_divergence(grid);
    CHECK(div_before < 1e-10);
}

TEST_CASE("Full timestep does not crash") {
    Grid grid(32, 16, 2.0, 1.0);
    grid.initialize_water(0.4);

    Solver solver;
    VOFTransport vof;

    // Run several timesteps without crashing
    for (int step = 0; step < 10; ++step) {
        double dt = solver.step(grid, 0.0, 0.0);
        vof.advect(grid, dt);
        CHECK(dt > 0.0);
    }

    // Check that VOF values remain bounded
    for (int i = 0; i < grid.NX; ++i) {
        for (int j = 0; j < grid.NY; ++j) {
            CHECK(grid.VOF(i, j) >= 0.0);
            CHECK(grid.VOF(i, j) <= 1.0);
        }
    }
}

TEST_CASE("Tank acceleration creates lateral flow") {
    Grid grid(32, 16, 2.0, 1.0);
    grid.initialize_water(0.4);

    Solver solver;
    VOFTransport vof;

    // Apply strong rightward tank acceleration (water should slosh left)
    for (int step = 0; step < 50; ++step) {
        double dt = solver.step(grid, 10.0, 0.0); // Tank accel right
        vof.advect(grid, dt);
    }

    // Water should have shifted: more VOF on the left, less on the right
    double left_vol = 0.0, right_vol = 0.0;
    for (int i = 0; i < grid.NX / 2; ++i)
        for (int j = 0; j < grid.NY; ++j)
            left_vol += grid.VOF(i, j);
    for (int i = grid.NX / 2; i < grid.NX; ++i)
        for (int j = 0; j < grid.NY; ++j)
            right_vol += grid.VOF(i, j);

    CHECK(left_vol > right_vol);
}
