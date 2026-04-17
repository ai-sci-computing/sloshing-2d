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
 * @enum InitialCondition
 * @brief Selectable initial conditions for the tank.
 */
enum class InitialCondition {
    StillWater,      ///< Flat water surface at fill_fraction
    Dambreak,        ///< Water column on the left, dry tank elsewhere
    DambreakBaffle,  ///< Dambreak with a vertical baffle in the right half
    TiltedSurface,   ///< First sloshing mode: linearly tilted surface
    CentralColumn,   ///< Water column in the centre, symmetric release
    TwoColumns,      ///< Two columns separated by a gap, colliding waves
    Shallow,         ///< Low fill for thin-film dynamics
    DambreakBaffle34, ///< Dambreak with baffle at 3/4 tank width
    OpenFOAMDambreak  ///< Exact OpenFOAM damBreak tutorial geometry
};

/**
 * @struct SimConfig
 * @brief Configuration for the simulation.
 */
struct SimConfig {
    int nx = 256;              ///< Grid cells in x
    int ny = 128;              ///< Grid cells in y
    double lx = 2.0;           ///< Tank width [m]
    double ly = 1.0;           ///< Tank height [m]
    double fill_fraction = 0.4; ///< Initial water fill level (still-water IC)
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
    int  ic_key_held_ = 0;      ///< Debounce for IC selection keys
    bool tension_key_held_ = false; ///< Debounce for T (surface tension toggle)
    bool tension_on_ = true;    ///< Surface tension currently active?
    double tension_default_ = 0.072; ///< Default σ [N/m] captured at init
    InitialCondition current_ic_ = InitialCondition::StillWater;

    /// @brief Apply the given initial condition and reset sim time / volume.
    void apply_initial_condition(InitialCondition ic);

    /// @brief Short label for the title bar.
    static const char* ic_label(InitialCondition ic);

    /**
     * @brief Process one frame: handle input, substep simulation, render.
     * @param frame_dt Real-time elapsed since last frame [s]
     */
    void frame_step(double frame_dt);
};
