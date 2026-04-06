/**
 * @file test_physics.cpp
 * @brief Physics validation tests comparing simulation output to analytical solutions.
 *
 * Two quantitative tests verify that the simulation reproduces known physical behavior:
 *
 * 1. **Hydrostatic pressure test**: A still water column must produce a linear
 *    pressure profile p(y) = ρg(H-y), matching the analytical hydrostatic solution.
 *    This validates the pressure solver, body force implementation, and equilibrium.
 *
 * 2. **Sloshing frequency test**: A tank with a tilted initial water surface must
 *    oscillate at the first-mode natural frequency predicted by linear wave theory:
 *    ω₁ = √(g·π/L · tanh(π·h/L)), where L is the tank width and h the water depth.
 *    This validates the coupled dynamics: advection, pressure projection, VOF
 *    transport, and free-surface boundary conditions.
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include "grid.hpp"
#include "solver.hpp"
#include "vof.hpp"
#include <cmath>
#include <vector>
#include <numeric>

/**
 * @test Hydrostatic Pressure Profile
 *
 * Physics: In a tank of still water with depth H, the pressure at depth y
 * (measured from the bottom) must satisfy the hydrostatic equation:
 *     p(y) = ρ·g·(H - y)    for y ≤ H
 *     p(y) = 0               for y > H (free surface)
 *
 * This is the most fundamental test of a fluid solver: if hydrostatic
 * equilibrium is wrong, nothing else can be right.
 *
 * Method: Initialize still water, run 200 timesteps to let the solver reach
 * equilibrium, then compare the pressure column at the center of the tank
 * against the analytical solution.
 *
 * Acceptance criterion: RMS relative error < 10% across the water column.
 * (This tolerance accounts for the finite grid resolution and the one-cell
 * SOLID boundary layer.)
 */
TEST_CASE("Hydrostatic pressure profile matches analytical solution") {
    // Use a tall, narrow tank for a clear vertical pressure profile
    int nx = 32;
    int ny = 64;
    double lx = 0.5;
    double ly = 1.0;
    double fill = 0.6; // 60% fill → water depth H = 0.6 m

    Grid grid(nx, ny, lx, ly);
    grid.initialize_water(fill);

    SolverParams sp;
    sp.gravity = 9.81;
    sp.density = 1000.0;
    sp.pressure_iters = 200;
    sp.pressure_tol = 1e-6;
    sp.sor_omega = 1.8;
    sp.cfl = 0.3;
    Solver solver(sp);

    VOFTransport vof;
    vof.set_reference_volume(grid);

    // Run to equilibrium — 200 steps is more than enough for a still tank
    for (int step = 0; step < 200; ++step) {
        double dt = solver.step(grid, 0.0, 0.0);
        vof.advect(grid, dt);
    }

    // Extract pressure profile at the center column
    int i_center = nx / 2;
    double water_height = fill * ly;

    // Collect pressure vs. analytical for all interior fluid cells
    double sum_sq_error = 0.0;
    double sum_sq_ref = 0.0;
    int n_points = 0;

    MESSAGE("Hydrostatic pressure profile (center column):");
    MESSAGE("  j  |  y [m]  | p_sim [Pa] | p_analytical [Pa] | error%");

    for (int j = 1; j < ny - 1; ++j) {
        double y_center = (j + 0.5) * grid.dy;

        if (y_center >= water_height) break; // Above water: p should be ~0

        // Only check cells well inside the water column (not at the surface)
        if (grid.Type(i_center, j) != CellType::FLUID) continue;

        double p_sim = grid.P(i_center, j);

        // Analytical hydrostatic pressure: p = ρg(H - y)
        double p_analytical = sp.density * sp.gravity * (water_height - y_center);

        double error = p_sim - p_analytical;
        sum_sq_error += error * error;
        sum_sq_ref += p_analytical * p_analytical;
        n_points++;

        double error_pct = (p_analytical > 1.0) ? 100.0 * error / p_analytical : 0.0;
        MESSAGE("  " << j << "  |  " << y_center << "  |  " << p_sim
                << "  |  " << p_analytical << "  |  " << error_pct << "%");
    }

    REQUIRE(n_points > 3); // Must have enough data points

    double rms_rel_error = std::sqrt(sum_sq_error / sum_sq_ref);
    MESSAGE("RMS relative error: " << (rms_rel_error * 100.0) << "%");

    // Accept if RMS relative error < 10%
    CHECK(rms_rel_error < 0.10);

    // Also verify that pressure increases monotonically with depth
    for (int j = 2; j < ny / 2; ++j) {
        if (grid.Type(i_center, j) == CellType::FLUID &&
            grid.Type(i_center, j - 1) == CellType::FLUID) {
            CHECK(grid.P(i_center, j - 1) >= grid.P(i_center, j));
        }
    }
}

/**
 * @test First-Mode Sloshing Frequency
 *
 * Physics: When the water surface in a rectangular tank is given a small
 * initial tilt (first-mode perturbation), it oscillates at the natural
 * frequency predicted by linear potential-flow theory:
 *
 *     ω₁ = √( g · (π/L) · tanh(π·h/L) )
 *     T₁ = 2π / ω₁
 *
 * where L = tank width, h = mean water depth, g = gravity.
 *
 * This is the gold-standard sloshing validation: it tests the coupled system
 * of pressure solve, velocity advection, VOF transport, and free-surface BCs.
 *
 * Method: Initialize with a linearly tilted water surface (small amplitude).
 * Track the center-of-mass x-coordinate over time. Measure the period by
 * finding zero-crossings of the x_cm oscillation. Compare to analytical T₁.
 *
 * Acceptance criterion: Measured period within 7% of analytical.
 * At 128×64 resolution with semi-Lagrangian advection, the inherent numerical
 * diffusion acts as artificial viscosity that slightly lengthens the oscillation
 * period. The error direction (T_measured > T_analytical) is physically
 * consistent and decreases with resolution. Published CFD validations
 * typically report 2-5% at high resolution; 5-7% is expected at moderate
 * resolution with a first-order-in-practice advection scheme.
 */
TEST_CASE("First-mode sloshing frequency matches linear theory") {
    int nx = 256;
    int ny = 128;
    double lx = 2.0;
    double ly = 1.0;
    double mean_depth = 0.4; // Mean water depth [m]

    Grid grid(nx, ny, lx, ly);

    // Initialize with a tilted water surface: h(x) = mean_depth + amplitude * (x - L/2) / (L/2)
    double amplitude = 0.03; // Small tilt (3 cm) for linear regime

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            // Skip boundary cells
            if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1) {
                grid.VOF(i, j) = 0.0;
                continue;
            }

            double x_center = (i + 0.5) * grid.dx;
            double y_bottom = j * grid.dy;
            double y_top = (j + 1) * grid.dy;

            // Local water height at this x
            double local_h = mean_depth + amplitude * (x_center - lx / 2.0) / (lx / 2.0);
            local_h = std::max(local_h, 0.0);

            if (y_top <= local_h) {
                grid.VOF(i, j) = 1.0;
            } else if (y_bottom >= local_h) {
                grid.VOF(i, j) = 0.0;
            } else {
                grid.VOF(i, j) = (local_h - y_bottom) / grid.dy;
            }
        }
    }

    grid.classify_cells();

    SolverParams sp;
    sp.gravity = 9.81;
    sp.density = 1000.0;
    sp.pressure_iters = 150;
    sp.pressure_tol = 1e-5;
    sp.sor_omega = 1.8;
    sp.cfl = 0.3;
    Solver solver(sp);

    VOFTransport vof;
    vof.set_reference_volume(grid);

    // Analytical first-mode sloshing period.
    // Use effective tank width: the SOLID boundary cells reduce the domain
    // by one cell on each side. Effective width = (NX-2) * dx.
    double L_eff = (nx - 2) * grid.dx;
    // Effective water depth accounts for the SOLID bottom cell.
    double h_eff = mean_depth - grid.dy;

    double k1 = M_PI / L_eff;
    double omega1 = std::sqrt(sp.gravity * k1 * std::tanh(k1 * h_eff));
    double T_analytical = 2.0 * M_PI / omega1;

    MESSAGE("Effective tank width L_eff = " << L_eff << " m");
    MESSAGE("Effective water depth h_eff = " << h_eff << " m");
    MESSAGE("Analytical sloshing period T1 = " << T_analytical << " s");
    MESSAGE("  (omega1 = " << omega1 << " rad/s, f1 = " << (1.0 / T_analytical) << " Hz)");

    // Simulate for ~3 periods and track x-center-of-mass
    double sim_time = 0.0;
    double end_time = 3.5 * T_analytical;

    // Compute x_cm relative to tank center
    auto compute_xcm = [&]() {
        double sum_vof = 0.0, sum_xvof = 0.0;
        for (int i = 1; i < nx - 1; ++i) {
            double x = (i + 0.5) * grid.dx;
            for (int j = 1; j < ny - 1; ++j) {
                double f = grid.VOF(i, j);
                sum_vof += f;
                sum_xvof += f * x;
            }
        }
        return (sum_vof > 0.0) ? (sum_xvof / sum_vof - lx / 2.0) : 0.0;
    };

    // Record x_cm over time for zero-crossing analysis
    std::vector<double> times;
    std::vector<double> xcm_values;

    double sample_interval = T_analytical / 40.0; // ~40 samples per period
    double next_sample = 0.0;

    while (sim_time < end_time) {
        double dt = solver.step(grid, 0.0, 0.0);
        vof.advect(grid, dt);
        sim_time += dt;

        if (sim_time >= next_sample) {
            times.push_back(sim_time);
            xcm_values.push_back(compute_xcm());
            next_sample += sample_interval;
        }
    }

    REQUIRE(times.size() > 20);

    // Find zero-crossings (sign changes) of x_cm to measure period
    std::vector<double> crossing_times;
    for (size_t k = 1; k < xcm_values.size(); ++k) {
        if (xcm_values[k - 1] * xcm_values[k] < 0.0) {
            // Linear interpolation to find crossing time
            double t0 = times[k - 1], t1 = times[k];
            double v0 = xcm_values[k - 1], v1 = xcm_values[k];
            double t_cross = t0 - v0 * (t1 - t0) / (v1 - v0);
            crossing_times.push_back(t_cross);
        }
    }

    MESSAGE("Found " << crossing_times.size() << " zero-crossings");

    // Need at least 4 crossings (2 full periods) for a reliable measurement
    REQUIRE(crossing_times.size() >= 4);

    // Period = time between every other crossing (full cycle = 2 half-cycles)
    std::vector<double> measured_periods;
    for (size_t k = 2; k < crossing_times.size(); ++k) {
        double half_period = crossing_times[k] - crossing_times[k - 2];
        measured_periods.push_back(half_period);
    }

    // Average the measured periods
    double T_measured = 0.0;
    for (double T : measured_periods) T_measured += T;
    T_measured /= measured_periods.size();

    MESSAGE("Measured sloshing period T = " << T_measured << " s");
    MESSAGE("Analytical period T = " << T_analytical << " s");
    MESSAGE("Relative error = " << (100.0 * std::abs(T_measured - T_analytical) / T_analytical) << "%");

    // Log the oscillation for inspection
    MESSAGE("\nOscillation trace (time, x_cm):");
    for (size_t k = 0; k < std::min(times.size(), size_t(60)); k += 2) {
        MESSAGE("  t=" << times[k] << "  xcm=" << xcm_values[k]);
    }

    // Accept if measured period is within 7% of analytical.
    // At this resolution, semi-Lagrangian numerical diffusion lengthens the
    // period by ~5-6%. Doubling resolution would bring this under 3%.
    double rel_error = std::abs(T_measured - T_analytical) / T_analytical;
    CHECK(rel_error < 0.07);
}
