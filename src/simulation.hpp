/**
 * @file simulation.hpp
 * @brief Top-level simulation orchestrator tying together solver, VOF, and renderer.
 *
 * Manages the main loop: process input → compute simulation substeps → render.
 * Tank motion is converted from keyboard input to body forces applied in the
 * fluid's reference frame.
 */

#pragma once

#include "grid.hpp"
#include "solver.hpp"
#include "vof.hpp"
#include "renderer.hpp"
#include <memory>

/**
 * @struct SimConfig
 * @brief Configuration for the simulation.
 */
struct SimConfig {
    int nx = 256;              ///< Grid cells in x
    int ny = 128;              ///< Grid cells in y
    double lx = 2.0;           ///< Tank width [m]
    double ly = 1.0;           ///< Tank height [m]
    double fill_fraction = 0.4; ///< Initial water fill level
    double tank_accel = 15.0;  ///< Tank acceleration magnitude [m/s²]
    int win_width  = 1024;     ///< Window width [px]
    int win_height = 512;      ///< Window height [px]
    double speed_multiplier = 1.0; ///< Simulation speed scaling
};

/**
 * @class Simulation
 * @brief Orchestrates the sloshing tank simulation.
 */
class Simulation {
public:
    /**
     * @brief Initialize the simulation with the given configuration.
     * @param config Simulation parameters
     * @return true on success
     */
    bool init(const SimConfig& config);

    /**
     * @brief Run the main simulation loop until the user quits.
     */
    void run();

    /**
     * @brief Clean up resources.
     */
    void shutdown();

private:
    SimConfig config_;
    std::unique_ptr<Grid> grid_;
    Solver solver_;
    VOFTransport vof_transport_;
    Renderer renderer_;

    double sim_time_ = 0.0;     ///< Accumulated simulation time [s]
    bool paused_ = false;       ///< Simulation paused flag
    bool pause_toggled_ = false; ///< Debounce for pause key

    /**
     * @brief Process one frame: handle input, substep simulation, render.
     * @param frame_dt Real-time elapsed since last frame [s]
     */
    void frame_step(double frame_dt);
};
