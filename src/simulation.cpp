/**
 * @file simulation.cpp
 * @brief Implementation of the top-level simulation orchestrator.
 */

#include "simulation.hpp"
#include <memory>

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif

#include <GLFW/glfw3.h>
#include <algorithm>
#include <cstdio>

const char* Simulation::ic_label(InitialCondition ic) {
    switch (ic) {
        case InitialCondition::StillWater:     return "still water";
        case InitialCondition::Dambreak:       return "dambreak";
        case InitialCondition::DambreakBaffle: return "dambreak + baffle";
        case InitialCondition::TiltedSurface:  return "tilted surface";
        case InitialCondition::CentralColumn:  return "central column";
        case InitialCondition::TwoColumns:     return "two columns";
        case InitialCondition::Shallow:        return "shallow";
        case InitialCondition::DambreakBaffle34: return "dambreak + baffle 3/4";
        case InitialCondition::OpenFOAMDambreak: return "OpenFOAM damBreak";
    }
    return "unknown";
}

void Simulation::apply_initial_condition(InitialCondition ic) {
    switch (ic) {
        case InitialCondition::StillWater:
            grid_->initialize_water(config_.fill_fraction);
            break;
        case InitialCondition::Dambreak:
            grid_->initialize_dambreak();
            break;
        case InitialCondition::DambreakBaffle:
            grid_->initialize_dambreak_baffle();
            break;
        case InitialCondition::TiltedSurface:
            grid_->initialize_tilted_surface();
            break;
        case InitialCondition::CentralColumn:
            grid_->initialize_central_column();
            break;
        case InitialCondition::TwoColumns:
            grid_->initialize_two_columns();
            break;
        case InitialCondition::Shallow:
            grid_->initialize_shallow();
            break;
        case InitialCondition::DambreakBaffle34:
            grid_->initialize_dambreak_baffle(0.25, 0.4, 0.75, 0.15, 2);
            break;
        case InitialCondition::OpenFOAMDambreak:
            // Rebuild grid to match the OpenFOAM square domain (0.584 x 0.584 m)
            grid_ = std::make_unique<Grid>(config_.nx, config_.nx, 0.584, 0.584);
            // column: 0.1461 wide x 0.292 tall; baffle at x=0.292, height 0.048
            grid_->initialize_dambreak_baffle(
                0.1461 / 0.584,   // column_width_frac ≈ 0.25
                0.292  / 0.584,   // column_height_frac = 0.50
                0.292  / 0.584,   // baffle_x_frac = 0.50
                0.048  / 0.584,   // baffle_height_frac ≈ 0.082
                2);
            break;
    }
    vof_transport_.set_reference_volume(*grid_);
    sim_time_ = 0.0;
    current_ic_ = ic;
    std::printf("Initial condition: %s\n", ic_label(ic));
}

bool Simulation::init(const SimConfig& config) {
    config_ = config;

    grid_ = std::make_unique<Grid>(config_.nx, config_.ny, config_.lx, config_.ly);

    SolverParams sp;
    sp.pressure_iters = 60;
    sp.pressure_tol = 1e-4;
    sp.sor_omega = 1.85;
    sp.cfl = 0.5;
    solver_ = Solver(sp);
    tension_default_ = solver_.params.surface_tension;
    tension_on_ = true;

    apply_initial_condition(InitialCondition::StillWater);

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
    std::printf("  R                 : Reset (re-apply current IC)\n");
    std::printf("  1                 : Still water\n");
    std::printf("  2                 : Dambreak\n");
    std::printf("  3                 : Dambreak with baffle\n");
    std::printf("  4                 : Tilted surface (sloshing mode)\n");
    std::printf("  5                 : Central column\n");
    std::printf("  6                 : Two-column collision\n");
    std::printf("  7                 : Shallow fill\n");
    std::printf("  8                 : Dambreak with baffle at 3/4\n");
    std::printf("  9                 : OpenFOAM damBreak tutorial\n");
    std::printf("  T                 : Toggle surface tension on/off\n");
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
            apply_initial_condition(current_ic_);
        }

        // Debounced IC selection (number keys). Only re-apply on the rising edge.
        if (input.ic_select != 0) {
            if (input.ic_select != ic_key_held_) {
                InitialCondition ic = current_ic_;
                if      (input.ic_select == 1) ic = InitialCondition::StillWater;
                else if (input.ic_select == 2) ic = InitialCondition::Dambreak;
                else if (input.ic_select == 3) ic = InitialCondition::DambreakBaffle;
                else if (input.ic_select == 4) ic = InitialCondition::TiltedSurface;
                else if (input.ic_select == 5) ic = InitialCondition::CentralColumn;
                else if (input.ic_select == 6) ic = InitialCondition::TwoColumns;
                else if (input.ic_select == 7) ic = InitialCondition::Shallow;
                else if (input.ic_select == 8) ic = InitialCondition::DambreakBaffle34;
                else if (input.ic_select == 9) ic = InitialCondition::OpenFOAMDambreak;
                apply_initial_condition(ic);
                ic_key_held_ = input.ic_select;
            }
        } else {
            ic_key_held_ = 0;
        }

        // Debounced pause toggle
        if (input.pause && !pause_toggled_) {
            paused_ = !paused_;
            pause_toggled_ = true;
        } else if (!input.pause) {
            pause_toggled_ = false;
        }

        // Debounced surface-tension toggle
        if (input.toggle_tension && !tension_key_held_) {
            tension_on_ = !tension_on_;
            solver_.params.surface_tension = tension_on_ ? tension_default_ : 0.0;
            std::printf("Surface tension: %s\n", tension_on_ ? "on" : "off");
            tension_key_held_ = true;
        } else if (!input.toggle_tension) {
            tension_key_held_ = false;
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
            int max_substeps = 10;

            for (int s = 0; s < max_substeps && accumulated < budget; ++s) {
                double dt = solver_.step(*grid_, tank_ax, tank_ay);
                vof_transport_.advect(*grid_, dt);
                sim_time_ += dt;
                accumulated += dt;
            }
        }

        // Render
        renderer_.render(*grid_, sim_time_, fps,
                         solver_.last_max_divergence,
                         ic_label(current_ic_), tension_on_);
    }
}

void Simulation::shutdown() {
    renderer_.shutdown();
    grid_.reset();
}
