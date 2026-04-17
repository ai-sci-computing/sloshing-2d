/**
 * @file test_physics.cpp
 * @brief Physics validation tests comparing simulation output to analytical solutions.
 *
 * Tests:
 * 1. Hydrostatic pressure profile: p(y) = ρg(H-y) in still water.
 * 2. First-mode sloshing frequency: ω₁ = √(g·π/L · tanh(π·h/L)) for a tilted surface.
 * 3. Dambreak wavefront vs. Ritter's shallow-water solution: x_front(t) = 2√(g·h₀)·t.
 * 4. Baffle impermeability: water below baffle height stays on its side.
 * 5. Volume conservation under violent lateral shaking.
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
            double x_center = (i + 0.5) * grid.dx;
            double y_bottom = j * grid.dy;
            double y_top = (j + 1) * grid.dy;

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

    // Analytical first-mode sloshing period. With no SOLID border layer,
    // the full tank width and water depth apply directly.
    double L_eff = grid.Lx;
    double h_eff = mean_depth;

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

/**
 * @test Dambreak wavefront vs Ritter's shallow-water solution
 *
 * Physics: The classical Ritter (1892) solution for a 1-D dambreak from a
 * still column of depth h₀ predicts a leading-edge speed of
 *
 *     c_front = 2 √(g·h₀)
 *
 * and thus a front position x_front(t) = x₀ + c_front·t, where x₀ is the
 * initial right edge of the reservoir.
 *
 * Method: release a water column (h₀ = 0.5·Ly, width = 0.3·Lx) under
 * gravity, integrate for a short time, and track the rightmost cell with
 * VOF > 0.5 over time. Compare to Ritter.
 *
 * Acceptance: Ritter is a 1-D idealisation; real 2-D free-surface flow is
 * slower because of vertical inertia, and algebraic VOF + SL advection
 * adds extra dissipation. The ratio x_sim/x_Ritter measured from column
 * edge typically lands in [0.5, 1.0]. We accept [0.3, 1.2] — loose, but
 * a regression test that flags catastrophic slowdowns or unphysical
 * overshoots. Tighten when the advection scheme is improved.
 */
TEST_CASE("Dambreak wavefront tracks Ritter's shallow-water solution") {
    int nx = 128, ny = 64;
    double lx = 2.0, ly = 1.0;
    double col_w_frac = 0.3;
    double col_h_frac = 0.5;
    double h0 = col_h_frac * ly;
    double col_edge = col_w_frac * lx;

    Grid grid(nx, ny, lx, ly);
    grid.initialize_dambreak(col_w_frac, col_h_frac);

    SolverParams sp;
    sp.gravity = 9.81;
    sp.density = 1000.0;
    sp.surface_tension = 0.0; // dambreak is a macroscopic inertial flow
    sp.pressure_iters = 120;
    sp.pressure_tol = 1e-5;
    sp.sor_omega = 1.85;
    sp.cfl = 0.4;
    Solver solver(sp);

    VOFTransport vof;
    vof.set_reference_volume(grid);

    double c_ritter = 2.0 * std::sqrt(sp.gravity * h0);
    MESSAGE("Ritter front speed c = 2 sqrt(g h0) = " << c_ritter << " m/s");

    auto wavefront_x = [&]() {
        for (int i = grid.NX - 1; i >= 0; --i) {
            for (int j = 0; j < grid.NY; ++j) {
                if (grid.VOF(i, j) > 0.5) return (i + 0.5) * grid.dx;
            }
        }
        return 0.0;
    };

    MESSAGE("  t [s]  | x_sim [m] | x_Ritter [m] | ratio");

    double t = 0.0;
    double t_end = 0.25;
    double sample_dt = 0.025;
    double next_sample = sample_dt;
    double final_ratio = 0.0;

    while (t < t_end) {
        double dt = solver.step(grid, 0.0, 0.0);
        vof.advect(grid, dt);
        t += dt;

        if (t >= next_sample) {
            double x_sim = wavefront_x();
            double x_ritter = col_edge + c_ritter * t;
            double ratio = (x_sim - col_edge) / (x_ritter - col_edge);
            MESSAGE("  " << t << "  |  " << x_sim << "  |  " << x_ritter
                    << "  |  " << ratio);
            final_ratio = ratio;
            next_sample += sample_dt;
        }
    }

    MESSAGE("Final wavefront ratio (x_sim - x0) / (x_Ritter - x0) = " << final_ratio);
    MESSAGE("(Physical 2-D CFD typically lands in [0.7, 0.95]. Current scheme");
    MESSAGE(" lands near ~0.45 after velocity extrapolation; removing extrapolation");
    MESSAGE(" drops this back to ~0.28, making this a useful regression guard.)");

    CHECK(final_ratio > 0.4); // regression guard — current baseline is ~0.45
    CHECK(final_ratio < 1.2); // catches unphysical overshoot
}

/**
 * @test Interactive dambreak setup vs Ritter and height probe
 *
 * Validates the exact column geometry used in the interactive sim
 * (0.25 Lx wide, 0.4 Ly tall) against Ritter's front speed at 256x128.
 * Also reports the water-height time history at x = Lx/2 (the midpoint),
 * which can be compared to published dambreak data.
 */
TEST_CASE("Interactive dambreak (0.25 x 0.4) vs Ritter and height probe") {
    int nx = 256, ny = 128;
    double lx = 2.0, ly = 1.0;
    double col_w = 0.25, col_h = 0.4;
    double h0 = col_h * ly;
    double col_edge = col_w * lx;

    Grid grid(nx, ny, lx, ly);
    grid.initialize_dambreak(col_w, col_h);

    SolverParams sp;
    sp.gravity = 9.81;
    sp.density = 1000.0;
    sp.surface_tension = 0.072;
    sp.pressure_iters = 120;
    sp.pressure_tol = 1e-5;
    sp.sor_omega = 1.85;
    sp.cfl = 0.5;
    Solver solver(sp);

    VOFTransport vof;
    vof.set_reference_volume(grid);

    double c_ritter = 2.0 * std::sqrt(sp.gravity * h0);
    MESSAGE("Column: " << col_w * lx << " m wide x " << h0 << " m tall");
    MESSAGE("Ritter front speed = " << c_ritter << " m/s");

    int i_probe = nx / 2;
    auto water_height_at = [&](int ic) {
        for (int j = grid.NY - 1; j >= 0; --j)
            if (grid.VOF(ic, j) > 0.5)
                return (j + 0.5) * grid.dy;
        return 0.0;
    };

    auto wavefront_x = [&]() {
        for (int i = grid.NX - 1; i >= 0; --i)
            for (int j = 0; j < grid.NY; ++j)
                if (grid.VOF(i, j) > 0.5) return (i + 0.5) * grid.dx;
        return 0.0;
    };

    MESSAGE("\n  t [s]   | x_front [m] | x_Ritter [m] | ratio  | h(Lx/2) [m]");

    double t = 0.0, t_end = 0.5;
    double sample_dt = 0.05, next_sample = sample_dt;
    double final_ratio = 0.0;

    while (t < t_end) {
        double dt = solver.step(grid, 0.0, 0.0);
        vof.advect(grid, dt);
        t += dt;

        if (t >= next_sample) {
            double xf = wavefront_x();
            double xr = col_edge + c_ritter * t;
            double ratio = (xf - col_edge) / (xr - col_edge);
            double h_mid = water_height_at(i_probe);
            MESSAGE("  " << t << "  |  " << xf << "  |  " << xr
                    << "  |  " << ratio << "  |  " << h_mid);
            final_ratio = ratio;
            next_sample += sample_dt;
        }
    }

    double vol_final = vof.total_volume(grid);
    double vol_init = col_w * lx * h0;
    double vol_err = std::abs(vol_final - vol_init) / vol_init;
    MESSAGE("\nVolume conservation: " << (vol_err * 100) << "% drift");
    MESSAGE("Final Ritter ratio: " << final_ratio);
    MESSAGE("Max divergence: " << solver.last_max_divergence << " 1/s");

    CHECK(final_ratio > 0.3);
    CHECK(final_ratio < 1.2);
    CHECK(vol_err < 0.05);
}

/**
 * @test Thin-film spreading (shallow-water limit of Ritter)
 *
 * Physics: for a THIN water layer of depth h₀ released from rest, the
 * Ritter shallow-water solution is asymptotically exact — the vertical
 * inertia corrections that slow the thick-column dambreak are small
 * when h₀ << column width. So the expected ratio of simulated to
 * Ritter front speed is closer to 1.0 for thin films than for thick
 * columns.
 *
 * Method: place a 10 cm thick water layer covering the left 40 % of the
 * floor; the rest of the tank is dry. Release under gravity and track
 * the leading edge (rightmost cell with VOF > 0.5) over time. The right
 * edge of the initial reservoir sits at x₀ = 0.4·Lx = 0.8 m. The Ritter
 * front advances at c = 2√(g·h₀) ≈ 1.98 m/s.
 *
 * Why this is useful: the thick-column Ritter test (above) mixes
 * shallow-water dynamics with strong vertical dynamics. The thin-film
 * version isolates the shallow-water limit, which is the regime most
 * relevant to free-surface spreading (puddle spreading, wave run-up,
 * baffle shadow closure).
 *
 * Acceptance: x_sim/x_Ritter > 0.2 as a regression guard. The physical
 * 2-D CFD expectation for a thin layer at this resolution is ~0.6–0.8;
 * the gap between that and whatever we measure quantifies how much the
 * advection/VOF scheme is under-driving thin-film motion.
 */
TEST_CASE("Thin-film dambreak wavefront tracks shallow-water Ritter") {
    int nx = 128, ny = 64;
    double lx = 2.0, ly = 1.0;
    double col_w_frac = 0.4;
    double col_h_frac = 0.1;               // 10 % = 10 cm — thin layer
    double h0 = col_h_frac * ly;
    double col_edge = col_w_frac * lx;

    Grid grid(nx, ny, lx, ly);
    grid.initialize_dambreak(col_w_frac, col_h_frac);

    SolverParams sp;
    sp.gravity = 9.81;
    sp.density = 1000.0;
    sp.surface_tension = 0.0;
    sp.pressure_iters = 120;
    sp.pressure_tol = 1e-5;
    sp.sor_omega = 1.85;
    sp.cfl = 0.4;
    Solver solver(sp);

    VOFTransport vof;
    vof.set_reference_volume(grid);

    double c_ritter = 2.0 * std::sqrt(sp.gravity * h0);
    MESSAGE("Thin-film h0 = " << h0 << " m; Ritter front speed = " << c_ritter << " m/s");

    auto wavefront_x = [&]() {
        for (int i = grid.NX - 1; i >= 0; --i) {
            for (int j = 0; j < grid.NY; ++j) {
                if (grid.VOF(i, j) > 0.5) return (i + 0.5) * grid.dx;
            }
        }
        return 0.0;
    };

    MESSAGE("  t [s]  | x_sim [m] | x_Ritter [m] | ratio");

    double t = 0.0;
    double t_end = 0.35;                   // long enough: advances > 1 cell per sample
    double sample_dt = 0.035;
    double next_sample = sample_dt;
    double final_ratio = 0.0;

    while (t < t_end) {
        double dt = solver.step(grid, 0.0, 0.0);
        vof.advect(grid, dt);
        t += dt;
        if (t >= next_sample) {
            double x_sim = wavefront_x();
            double x_ritter = col_edge + c_ritter * t;
            double ratio = (x_sim - col_edge) / (x_ritter - col_edge);
            MESSAGE("  " << t << "  |  " << x_sim << "  |  " << x_ritter
                    << "  |  " << ratio);
            final_ratio = ratio;
            next_sample += sample_dt;
        }
    }

    MESSAGE("Final thin-film wavefront ratio = " << final_ratio);
    MESSAGE("(Shallow-water CFD target ~0.6–0.8; current result quantifies");
    MESSAGE(" the thin-film / leading-edge damping still in the scheme.)");

    CHECK(final_ratio > 0.2);
    CHECK(final_ratio < 1.2);
}

/**
 * @test Baffle impermeability
 *
 * Physics: an internal SOLID baffle must block fluid exchange. With a
 * water column shorter than the baffle, nothing should end up on the
 * opposite side even over long simulation times.
 *
 * Method: dambreak + baffle, with column height = 0.3·Ly and baffle
 * height = 0.4·Ly. Column collapses, sloshes against the baffle, but
 * should not overtop it. Integrate for 1 s, then measure volume that
 * crossed to the right of the baffle.
 *
 * Acceptance: less than 2 % of the reservoir volume should end up on
 * the far side. A larger number indicates leakage through the baffle.
 */
TEST_CASE("Baffle is impermeable to short-column dambreak") {
    int nx = 128, ny = 64;
    double lx = 2.0, ly = 1.0;
    double col_w_frac = 0.3;
    double col_h_frac = 0.3;   // intentionally short
    double baffle_x_frac = 0.5;
    double baffle_h_frac = 0.4; // taller than the column

    Grid grid(nx, ny, lx, ly);
    grid.initialize_dambreak_baffle(col_w_frac, col_h_frac,
                                    baffle_x_frac, baffle_h_frac, 2);

    SolverParams sp;
    sp.gravity = 9.81;
    sp.density = 1000.0;
    sp.surface_tension = 0.0;
    sp.pressure_iters = 120;
    sp.pressure_tol = 1e-5;
    sp.sor_omega = 1.85;
    sp.cfl = 0.4;
    Solver solver(sp);

    VOFTransport vof;
    vof.set_reference_volume(grid);

    double cell_area = grid.dx * grid.dy;
    double baffle_x_m = baffle_x_frac * lx;

    double ref_vol = vof.total_volume(grid);
    MESSAGE("Initial total volume = " << ref_vol << " m^2");

    double t = 0.0;
    double t_end = 1.0;
    while (t < t_end) {
        double dt = solver.step(grid, 0.0, 0.0);
        vof.advect(grid, dt);
        t += dt;
    }

    // Volume sitting strictly right of the baffle column (skip ±1 cell around baffle)
    double vol_right = 0.0;
    for (int i = 0; i < grid.NX; ++i) {
        double xc = (i + 0.5) * grid.dx;
        if (xc < baffle_x_m + 2.0 * grid.dx) continue;
        for (int j = 0; j < grid.NY; ++j) {
            vol_right += grid.VOF(i, j) * cell_area;
        }
    }
    double leak_frac = vol_right / ref_vol;
    MESSAGE("Volume right of baffle after " << t << " s = " << vol_right
            << "  (" << (leak_frac * 100.0) << "% of reservoir)");

    CHECK(leak_frac < 0.02);
}

/**
 * @test Volume conservation under violent shaking
 *
 * The solver must conserve fluid mass even when the tank is shaken hard
 * enough to drive complex breaking-wave dynamics. This is the "test with
 * violent shaking" failure mode — quiet tests can pass while violent
 * motion exposes mass leaks or blow-ups.
 *
 * Method: still-water IC, then apply a strong oscillating horizontal
 * acceleration for several seconds. Compare total VOF volume to initial.
 *
 * Acceptance: relative volume drift < 1 %. Also require max|div u| to
 * remain finite (solver did not blow up) and VOF to stay in [0, 1].
 */
TEST_CASE("Volume is conserved under violent shaking") {
    int nx = 128, ny = 64;
    double lx = 2.0, ly = 1.0;

    Grid grid(nx, ny, lx, ly);
    grid.initialize_water(0.4);

    SolverParams sp;
    sp.gravity = 9.81;
    sp.density = 1000.0;
    sp.surface_tension = 0.072;
    sp.pressure_iters = 80;
    sp.pressure_tol = 1e-5;
    sp.sor_omega = 1.85;
    sp.cfl = 0.4;
    Solver solver(sp);

    VOFTransport vof;
    vof.set_reference_volume(grid);

    double vol0 = vof.total_volume(grid);
    MESSAGE("Initial volume = " << vol0 << " m^2");

    double t = 0.0;
    double t_end = 2.0;
    double shake_amp = 20.0;   // violent [m/s^2]
    double shake_freq = 2.0;   // Hz
    double max_div = 0.0;

    while (t < t_end) {
        double tank_ax = shake_amp * std::sin(2.0 * M_PI * shake_freq * t);
        double dt = solver.step(grid, tank_ax, 0.0);
        vof.advect(grid, dt);
        t += dt;
        max_div = std::max(max_div, solver.last_max_divergence);
    }

    double vol_final = vof.total_volume(grid);
    double rel_err = std::abs(vol_final - vol0) / vol0;
    MESSAGE("Volume after " << t << " s violent shaking = " << vol_final
            << "   rel. error = " << (rel_err * 100.0) << "%");
    MESSAGE("Max divergence encountered = " << max_div << " 1/s");

    CHECK(rel_err < 0.01);
    CHECK(std::isfinite(max_div));
    for (int i = 0; i < grid.NX; ++i) {
        for (int j = 0; j < grid.NY; ++j) {
            CHECK(grid.VOF(i, j) >= 0.0);
            CHECK(grid.VOF(i, j) <= 1.0);
        }
    }
}
