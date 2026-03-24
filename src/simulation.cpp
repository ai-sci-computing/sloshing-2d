/**
 * @file simulation.cpp
 * @brief Implementation of the top-level simulation orchestrator.
 */

#include "simulation.hpp"

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif

#include <GLFW/glfw3.h>
#include <algorithm>
#include <cstdio>

bool Simulation::init(const SimConfig& config) {
    config_ = config;

    grid_ = new Grid(config_.nx, config_.ny, config_.lx, config_.ly);
    grid_->initialize_water(config_.fill_fraction);
    vof_transport_.set_reference_volume(*grid_);

    SolverParams sp;
    sp.viscosity = 0.01;       // Effective eddy viscosity for visual simulation
    sp.pressure_iters = 100;   // Warm-start needs fewer iterations
    sp.pressure_tol = 1e-5;   // Tight tolerance for divergence-free field
    sp.sor_omega = 1.85;       // Aggressive relaxation for faster convergence
    sp.cfl = 0.4;
    solver_ = Solver(sp);

    if (!renderer_.init(config_.win_width, config_.win_height, "Sloshing Tank Simulation")) {
        return false;
    }

    std::printf("Sloshing Tank Simulation\n");
    std::printf("========================\n");
    std::printf("Grid: %d x %d  |  Tank: %.1f x %.1f m\n",
                config_.nx, config_.ny, config_.lx, config_.ly);
    std::printf("Controls:\n");
    std::printf("  Arrow keys / WASD : Move tank\n");
    std::printf("  Space             : Pause/Resume\n");
    std::printf("  R                 : Reset\n");
    std::printf("  +/-               : Speed up/down\n");
    std::printf("  Q / Esc           : Quit\n");

    return true;
}

void Simulation::run() {
    double last_time = glfwGetTime();
    double fps_timer = 0.0;
    int frame_count = 0;
    double fps = 0.0;

    while (!renderer_.should_close()) {
        double now = glfwGetTime();
        double frame_dt = now - last_time;
        last_time = now;

        // Cap frame_dt to avoid spiral of death
        frame_dt = std::min(frame_dt, 0.05);

        // FPS counter
        fps_timer += frame_dt;
        frame_count++;
        if (fps_timer >= 1.0) {
            fps = frame_count / fps_timer;
            fps_timer = 0.0;
            frame_count = 0;
        }

        // Process input
        InputState input = renderer_.poll_input();

        if (input.quit) break;

        if (input.reset) {
            grid_->initialize_water(config_.fill_fraction);
            vof_transport_.set_reference_volume(*grid_);
            sim_time_ = 0.0;
        }

        // Debounced pause toggle
        if (input.pause && !pause_toggled_) {
            paused_ = !paused_;
            pause_toggled_ = true;
        } else if (!input.pause) {
            pause_toggled_ = false;
        }

        // Speed control
        if (input.speed_up) {
            config_.speed_multiplier = std::min(config_.speed_multiplier * 1.02, 10.0);
        }
        if (input.speed_down) {
            config_.speed_multiplier = std::max(config_.speed_multiplier * 0.98, 0.1);
        }

        // Compute tank acceleration from input
        double tank_ax = 0.0, tank_ay = 0.0;
        if (input.left)  tank_ax -= config_.tank_accel;
        if (input.right) tank_ax += config_.tank_accel;
        if (input.up)    tank_ay += config_.tank_accel;
        if (input.down)  tank_ay -= config_.tank_accel;

        // Simulation substeps
        if (!paused_) {
            double budget = frame_dt * config_.speed_multiplier;
            double accumulated = 0.0;
            int max_substeps = 20;

            for (int s = 0; s < max_substeps && accumulated < budget; ++s) {
                double dt = solver_.step(*grid_, tank_ax, tank_ay);
                vof_transport_.advect(*grid_, dt);
                sim_time_ += dt;
                accumulated += dt;
            }
        }

        // Render
        renderer_.render(*grid_, sim_time_, fps);
    }
}

void Simulation::shutdown() {
    renderer_.shutdown();
    delete grid_;
    grid_ = nullptr;
}
