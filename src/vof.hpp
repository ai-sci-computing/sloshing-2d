/**
 * @file vof.hpp
 * @brief Volume of Fluid (VOF) transport for free-surface tracking.
 *
 * Uses operator-split advection with Weymouth-Yue divergence correction,
 * Superbee TVD flux limiter, interface compression for sharpness, and
 * global mass correction.
 */

#pragma once

#include "grid.hpp"

/**
 * @class VOFTransport
 * @brief Advects the VOF field with mass-conserving operator-split TVD scheme.
 */
class VOFTransport {
public:
    /**
     * @brief Set the initial (reference) fluid volume for mass correction.
     * @param grid The simulation grid with initial water configuration
     */
    void set_reference_volume(const Grid& grid);

    /**
     * @brief Advect the VOF field by one timestep.
     *
     * Pipeline: advect → sharpen interface → enforce boundaries → correct mass.
     *
     * @param grid The simulation grid (VOF field modified in place)
     * @param dt Timestep [s]
     */
    void advect(Grid& grid, double dt);

    /**
     * @brief Compute total fluid volume (sum of VOF * cell area).
     * @param grid The simulation grid
     * @return Total fluid volume [m²] (2D area)
     */
    double total_volume(const Grid& grid) const;

private:
    int step_count_ = 0;
    double reference_volume_ = -1;

    void advect_x(Grid& grid, double dt);
    void advect_y(Grid& grid, double dt);

    /**
     * @brief Sharpen the interface via algebraic compression.
     *
     * Applies one step of the compression equation ∂f/∂t + ∇·(f(1-f)n̂) = 0,
     * which pushes diffuse interface values toward 0 or 1 along the interface
     * normal direction. This counteracts numerical diffusion without moving
     * the interface position.
     *
     * @param grid The simulation grid
     */
    void sharpen_interface(Grid& grid);

    /**
     * @brief Zero VOF in SOLID boundary cells to prevent wall artifacts.
     * @param grid The simulation grid
     */
    void enforce_vof_boundaries(Grid& grid);

    /**
     * @brief Rescale VOF to restore the reference volume.
     * @param grid The simulation grid
     */
    void correct_mass(Grid& grid);
};
